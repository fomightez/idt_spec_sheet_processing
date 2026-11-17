#!/usr/bin/env python3
"""
Build one de-duplicated, name-sorted primer list from the OCR results of scanned IDT spec sheets, and
(optionally) keep only the primers that are NEW compared with an established, 
curated "main" list, such as a shared Google Sheet.

How the pieces fit together
    1. idt_spec_ocr.py reads the scanned IDT spec-sheet PDFs and writes one CSV per PDF.
    2. THIS script reads all of those CSVs, merges duplicates, and writes the combined list. (STAGE 1;meant before you have an established active file of Chan Lab primers)
    3. If you give it an established file (as `.tsv`) with the flag `--main`, 
    it also compares the result with that list: anything the main file already has 
    is left exactly as the main file has it, and only primers the main file does 
    not have are added. A second file containing just those new primers is written so they can be pasted at the bottom of the shared sheet. (STAGE 2;meant for when you have an established active file of Chan Lab primers)

Usage:
    python idt_primer_spec_sheets_combine_with_established_list_handling_option.py /path/to/folder
    python idt_primer_spec_sheets_combine_with_established_list_handling_option.py /path/to/folder -o primers.tsv
    python idt_primer_spec_sheets_combine_with_established_list_handling_option.py /path/to/folder --print          # also echo every entry
    python idt_primer_spec_sheets_combine_with_established_list_handling_option.py /path/to/folder --include-source # add "source: file p.N"
    python idt_primer_spec_sheets_combine_with_established_list_handling_option.py /path/to/folder --main shared_sheet.tsv
    python idt_primer_spec_sheets_combine_with_established_list_handling_option.py /path/to/folder --main shared_sheet.tsv --main-order

STAGE 1: the spec-sheet CSVs (always)
    1. Finds every *.csv in the folder and all sub-folders, skipping hidden sub-folders such as
       .ipynb_checkpoints. A CSV is used only if it has the OCR script's columns (name, sequence, tm_c);
       anything else is skipped and listed.
    2. Applies NAME_OVERRIDES (in the SETTINGS block): names you have listed are replaced with the name you want,
       before anything is merged.
    3. Merges duplicates so each primer appears once. Two rows are the same primer when their SEQUENCE is
       identical (ignoring upper/lower case and spaces), whatever they were called, for example the same
       primer ordered twice under slightly different names. One row is kept: the best-verified one (OK, then
       OK_LENGTH_ONLY, then REVIEW, then NO_SEQ); ties go to the earliest file path and page. The other names
       are added to the kept entry's comments ("also listed as: ..."). Rows with no sequence can only be
       merged with rows of the same name. The source CSV files are never changed or deleted.
    4. Sorts by primer name (case-insensitive, "natural" order, so OMS-2 comes before OMS-10).
    5. Writes a tab-separated file. The first line is a header (Sequence, Name, Comments), then one primer per line:

           SEQUENCE<TAB>NAME<TAB>21 nts; Tm= 50.2 °C; comments

       The comments part is only present when there is something to say, e.g. the row still needs a
       check against the scan, or duplicate copies disagreed on the Tm.

STAGE 2: the main file (only with `--main shared_sheet.tsv`)
    Use this when you already keep a curated list (for example a shared Google Sheet with the columns
    Sequence, Name, Comments, saved with File > Download > Tab-separated values) and only want NEW primers from
    the spec sheets added to it.
      - The main file is never altered. Its rows are copied to the output exactly as they are, including any
        extra columns, and its header row (if it has one) is reused. The script refuses to overwrite it.
      - A spec-sheet entry whose sequence is already in the main file is dropped: the main file wins, including
        its name and comments. Spec-sheet entries with no sequence are dropped if their name is already in the
        main file.
      - Only sequences the main file does not have are added, with the usual comments plus these flags:
            "same bases as main-file entry ... but different modifications"  (possible duplicate, the
                                                                             modification codes are spelled
                                                                             differently)
            "NAME ALREADY IN MAIN FILE with a different sequence"
      - Two files are written: the full list (main file + new, sorted by name; use --main-order to keep the
        main file's rows in their original order with the new ones below) and a "..._new_only.tsv" containing
        just the new rows, ready to paste at the bottom of the shared sheet.
      - The console report also warns about problems inside the main file (the same sequence on two rows, or the
        same name with different sequences) and lists main-file entries whose name differs from the spec sheets.

Settings you may want to change (SETTINGS block, just below the imports)
    DEGREE_SYMBOL, NAME_OVERRIDES, FLAG_ALTERNATE_NAMES, HEADER / WRITE_HEADER, TM_CONFLICT_TOLERANCE, STATUS_RANK.

Nothing is dropped silently: every merged name is listed in the comment and in the console report, rows with a
missing sequence are kept (and flagged), and rows that share a name but have different sequences are all kept
(and reported).
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# ================================================================ SETTINGS
# Symbol written before the "C" in Tm. 
# "\u00b0" is the true degree sign "°". (At one point, the OCR also had seen the character `º`(U+00BA, "º").)
DEGREE_SYMBOL = "\u00b0"

# A CSV must have these columns to be treated as OCR-script output.
REQUIRED_COLUMNS = {"name", "sequence", "tm_c"}

# Which duplicate to keep: lower rank wins.
STATUS_RANK = {"OK": 0, "OK_LENGTH_ONLY": 1, "REVIEW": 2, "NO_SEQ": 3}

# Duplicate copies whose Tm differs by more than this (degrees C) are flagged.
TM_CONFLICT_TOLERANCE = 0.05

# When same-sequence entries are merged, add "also listed as: <other names>" to the kept entry's comments.
FLAG_ALTERNATE_NAMES = True

# First line of the output files. (If the main file has its own header row, that one is used instead.)
HEADER = ("Sequence", "Name", "Comments")
WRITE_HEADER = True

# Hard-coded name overrides: {name as it appears in the CSVs: name you want in the output}.
# Applied to every row BEFORE duplicates are merged, so it works whether the name is the one that would have
# been kept or one of the alternates. Matching ignores upper/lower case and extra spaces. The replacement can be
# any text; it does not have to be one of the names already in the data. The console report lists every
# override that was used and warns about any that matched nothing (a typo, or a name OCR read differently).
NAME_OVERRIDES = {
    '18s rRNAR': '18s rRNA R',
    'Scramble Reverse': 'Scrambled Reverse',
    'PIK3CD 3°UTR FP1 (Sequencing)': "PIK3CD 3'UTR FP1 (Sequencing)",
    'PIK3CD 3’UTR Forward Primer': "PIK3CD 3'UTR Forward Primer",
    'PIK3CD 3’UTR FP2 (Sequencing)': "PIK3CD 3'UTR FP2 (Sequencing)",
    'PIK3CD 3’UTR FP3 (Sequencing)': "PIK3CD 3'UTR FP3 (Sequencing)",
    'PIK3CD 3’UTR Reverse Primer': "PIK3CD 3'UTR Reverse Primer",
    'PIK3CD 3’UTR RP2 (Sequencing)': "PIK3CD 3'UTR RP2 (Sequencing)",
    'UL122 qgPCRR': 'UL122 qPCR R',
    'IE1_Rev]': 'IE1_Rev1',
}
# ====================================================================================

COMMENTS_FOR_STATUS = {
    "OK_LENGTH_ONLY": "modified oligo: only length was verified, check against the scan",
    "REVIEW": "REVIEW: checks failed, verify against the scan",
    "NO_SEQ": "NO SEQUENCE FOUND, check the scan",
}


def natural_key(text):
    """Sort key so that 'OMS-2' < 'OMS-10' and case is ignored."""
    return [int(t) if t.isdigit() else t.casefold() for t in re.split(r"(\d+)", text)]


def name_key(name):
    """Comparison form of a name: lower-case with whitespace collapsed."""
    return " ".join(name.casefold().split())


def original_name(entry):
    """The name as it was written in the CSV, before any NAME_OVERRIDES replacement."""
    return entry.get("original_name", entry["name"])


def display_name(entry):
    """Name for the console report, noting when an override changed it."""
    original = entry.get("original_name")
    if original and original != entry["name"]:
        return f"{entry['name']}  [override of '{original}']"
    return entry["name"]


def clean_field(value):
    """Collapse whitespace so tabs/newlines in a cell can't break the TSV."""
    return " ".join((value or "").split())


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def read_csv_rows(path):
    """Return (rows, problem). rows is a list of dicts, or None if the file isn't OCR-script output."""
    for encoding in ("utf-8-sig", "cp1252"):
        try:
            with open(path, newline="", encoding=encoding) as fh:
                reader = csv.DictReader(fh)
                columns = {c.strip() for c in (reader.fieldnames or [])}
                if not REQUIRED_COLUMNS <= columns:
                    return None, "missing columns " + ", ".join(sorted(REQUIRED_COLUMNS - columns))
                return list(reader), None
        except UnicodeDecodeError:
            continue
        except csv.Error as exc:
            return None, f"not readable as CSV ({exc})"
    return None, "could not decode file"


def bare_key(seq):
    """Sequence with /modification/ codes removed and upper-cased: the bases only."""
    return re.sub(r"/[^/]*/", "", re.sub(r"\s+", "", seq)).upper()


def row_seq_key(cells):
    """Comparison form of a main-file row's sequence cell."""
    return re.sub(r"\s+", "", cells[0]).upper() if cells else ""


def row_name(cells):
    return clean_field(cells[1]) if len(cells) > 1 else ""


def read_main_file(path):
    """Read the ground-truth TSV. Returns (header_cells or None, rows) where each row is a list of cell strings
    kept exactly as written (only trailing empty cells and blank lines are dropped)."""
    raw = Path(path).read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):                  # UTF-16 (Excel's "Unicode Text")
        text = raw.decode("utf-16")
    else:
        try:
            text = raw.decode("utf-8-sig")                            # Google Sheets / Sublime
        except UnicodeDecodeError:
            text = raw.decode("cp1252", errors="replace")
    rows = []
    for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        cells = line.split("\t")
        while cells and not cells[-1].strip():
            cells.pop()
        if cells:
            rows.append(cells)
    header = None
    if rows and len(rows[0]) > 1 and rows[0][0].strip().casefold() == "sequence" \
            and rows[0][1].strip().casefold() == "name":
        header = rows.pop(0)
    return header, rows


def make_entry(row, csv_path, root):
    name = clean_field(row.get("name"))
    seq = re.sub(r"\s+", "", row.get("sequence") or "")
    nts = to_int(row.get("bases_reported"))
    if nts is None:
        nts = to_int(row.get("bases_counted"))
    if nts is None and seq:
        nts = len(re.sub(r"/[^/]*/", "", seq))          # fall back to counting the bases ourselves
    check = clean_field(row.get("check")) or "REVIEW"   # no status recorded -> treat as needing a look
    try:
        shown = str(csv_path.relative_to(root))
    except ValueError:
        shown = str(csv_path)
    return {
        "name": name,
        "seq": seq,
        "nts": nts,
        "tm": to_float(row.get("tm_c")),
        "check": check,
        "source": shown,
        "page": to_int(row.get("page")) or 0,
    }


def rank(entry):
    return (STATUS_RANK.get(entry["check"], 4), entry["source"], entry["page"])


def format_line(entry, include_source):
    nts = f"{entry['nts']} nts" if entry["nts"] is not None else "? nts"
    tm = f"Tm= {entry['tm']:.1f} {DEGREE_SYMBOL}C" if entry["tm"] is not None else "Tm= n/a"
    parts = [nts, tm] + entry["comments"]
    if include_source:
        page = f" p.{entry['page']}" if entry["page"] else ""
        parts.append(f"source: {entry['source']}{page}")
    return f"{entry['seq']}\t{entry['name']}\t{'; '.join(parts)}"


def main():
    ap = argparse.ArgumentParser(description="Combine IDT spec-sheet OCR CSVs into one de-duplicated, name-sorted primer list, optionally keeping only primers that are new compared with a curated main file.")
    ap.add_argument("folder", help="folder to search (including all sub-folders)")
    ap.add_argument("-o", "--output", help="output TSV path (default: <folder>/combined_primers.tsv)")
    ap.add_argument("--print", dest="echo", action="store_true", help="also print every final entry")
    ap.add_argument("--include-source", action="store_true",
                    help="append 'source: <csv> p.N' to each entry's comments")
    ap.add_argument("--encoding", default="utf-8",
                    help="encoding of the TSV (default utf-8; try utf-8-sig if Excel, or other spreadsheet software, garbles symbols.)")
    ap.add_argument("--main", dest="main_file",
                    help="ground-truth TSV (Sequence, Name, Comments), e.g. the shared Google Sheet saved as .tsv. "
                         "It is never altered; only sequences it lacks are added.")
    ap.add_argument("--main-order", dest="main_order", action="store_true",
                    help="with --main: keep the main file's rows in their original order and put new entries below "
                         "(default: sort everything by name)")
    ap.add_argument("--new-only-output",
                    help="with --main: path for the file listing only the new entries "
                         "(default: <output name>_new_only.tsv)")
    args = ap.parse_args()

    root = Path(args.folder).expanduser().resolve()
    if not root.is_dir():
        sys.exit(f"Not a folder: {root}")
    out_path = Path(args.output).expanduser() if args.output else root / "combined_primers.tsv"
    new_only_path = (Path(args.new_only_output).expanduser() if args.new_only_output
                     else out_path.with_name(out_path.stem + "_new_only" + (out_path.suffix or ".tsv")))
    main_path = None
    if args.main_file:
        main_path = Path(args.main_file).expanduser().resolve()
        if not main_path.is_file():
            sys.exit(f"Main file not found: {main_path}")
        if main_path in (out_path.resolve(), new_only_path.resolve()):
            sys.exit(f"The output would overwrite the main file ({main_path}). Choose a different -o name; "
                     "the main file must never be overwritten.")

    # ---------------------------------------------------------------- 1. read every CSV
    csv_files = sorted(
        p for p in root.rglob("*.csv")
        if p.is_file() and not any(part.startswith(".") for part in p.relative_to(root).parts[:-1])
    ) # find any `.csv` files that aren't in a hidden sub-directory, like 
    # `.ipynb_checkpoints`folders
    entries, skipped = [], []
    used_files = 0
    for path in csv_files:
        rows, problem = read_csv_rows(path)
        if rows is None:
            skipped.append((path, problem))
            continue
        used_files += 1
        entries.extend(make_entry(r, path, root) for r in rows)
    if not entries:
        sys.exit(f"No OCR-script CSV rows found under {root}")

    override_lookup = {name_key(k): (k, v) for k, v in NAME_OVERRIDES.items()}
    override_hits = defaultdict(int)
    for e in entries:
        hit = override_lookup.get(name_key(e["name"]))
        if hit:
            e["original_name"] = e["name"]
            e["name"] = clean_field(hit[1])
            override_hits[hit[0]] += 1

    # ----------------------------- 2. merge duplicates (same sequence = same primer)-----------#
    groups = defaultdict(list)
    for e in entries:
        if e["seq"]:
            groups[("seq", e["seq"].upper())].append(e)
        else:                                   # no sequence to compare, so only merge identical names
            groups[("name", name_key(e["name"]))].append(e)

    final, conflicts, merged_groups = [], [], []
    n_same_name = n_renamed = 0
    for copies in groups.values():
        copies.sort(key=rank)                   # best-verified first, then earliest file path and page
        keep = dict(copies[0])
        keep["comments"] = []
        keep["sheet_names"] = sorted({original_name(c) for c in copies}, key=natural_key)
        if keep["check"] in COMMENTS_FOR_STATUS:
            keep["comments"].append(COMMENTS_FOR_STATUS[keep["check"]])

        # Other names these rows had in the CSVs (before overrides): one representative row per distinct name.
        keep_orig = name_key(original_name(copies[0]))
        distinct = {}
        for c in copies[1:]:
            k = name_key(original_name(c))
            if k != keep_orig and k not in distinct:
                distinct[k] = c
        n_renamed += len(distinct)
        n_same_name += len(copies) - 1 - len(distinct)

        # "also listed as": every name the CSVs used for this primer EXCEPT the name that is actually shown.
        # This includes a name you replaced with NAME_OVERRIDES, but never the replacement itself.
        shown_key = name_key(keep["name"])
        listed = {}
        for c in copies:
            nm = original_name(c)
            if name_key(nm) != shown_key and name_key(nm) not in listed:
                listed[name_key(nm)] = nm
        listed_names = sorted(listed.values(), key=natural_key)

        tms = sorted({c["tm"] for c in copies if c["tm"] is not None})
        if tms and tms[-1] - tms[0] > TM_CONFLICT_TOLERANCE:
            shown = "/".join(f"{t:.1f}" for t in tms)
            keep["comments"].append(f"duplicate copies disagree on Tm ({shown}), check the scans")
            conflicts.append(f"{keep['name']}: Tm {shown}")
        nts_vals = sorted({c["nts"] for c in copies if c["nts"] is not None})
        if len(nts_vals) > 1:
            keep["comments"].append(f"duplicate copies disagree on length ({'/'.join(map(str, nts_vals))})")
            conflicts.append(f"{keep['name']}: length {'/'.join(map(str, nts_vals))}")
        if distinct:
            merged_groups.append((keep, sorted(distinct.values(), key=lambda c: natural_key(original_name(c)))))
        if listed_names and FLAG_ALTERNATE_NAMES:
            keep["comments"].append("also listed as: " + " | ".join(listed_names))
        final.append(keep)

    # ---------------------------------------------------------------- 3. sort the spec-sheet entries by name
    final.sort(key=lambda e: (natural_key(e["name"]), e["seq"].upper()))

    # ---------------------------------------------------------------- 4. compare with the main file (ground truth)
    main_header, main_rows = (read_main_file(main_path) if main_path else (None, []))
    first_row_number = 2 if main_header else 1            # row numbers as they appear in the sheet

    def sheet_row(i):
        return i + first_row_number

    main_seq_rows, main_bare_rows, main_name_rows = defaultdict(list), defaultdict(list), defaultdict(list)
    for i, cells in enumerate(main_rows):
        if row_seq_key(cells):
            main_seq_rows[row_seq_key(cells)].append(i)
            main_bare_rows[bare_key(cells[0])].append(i)
        if row_name(cells):
            main_name_rows[name_key(row_name(cells))].append(i)

    new_entries = final
    already_in_main, near_matches, name_clashes = [], [], []
    if main_path:
        new_entries = []
        for e in final:
            key = e["seq"].upper()
            if key and key in main_seq_rows:                         # main file wins
                already_in_main.append((e, main_seq_rows[key]))
                continue
            if not key and name_key(e["name"]) in main_name_rows:    # no sequence read: match on the name instead
                already_in_main.append((e, main_name_rows[name_key(e["name"])]))
                continue
            if key and bare_key(key) in main_bare_rows:              # same bases, modification codes differ
                rows_ = main_bare_rows[bare_key(key)]
                names_ = " | ".join(row_name(main_rows[i]) for i in rows_)
                e["comments"].append(f"same bases as main-file entry '{names_}' but different modifications: "
                                     "possible duplicate, check")
                near_matches.append((e, rows_))
            if name_key(e["name"]) in main_name_rows:
                e["comments"].append("NAME ALREADY IN MAIN FILE with a different sequence, check")
                name_clashes.append((e, main_name_rows[name_key(e["name"])]))
            new_entries.append(e)

    # ---------------------------------------------------------------- 5. write the files
    def sort_rows(rows):
        return [line for _, _, line in sorted(rows, key=lambda r: (r[0], r[1]))]

    new_rows = [(natural_key(e["name"]), e["seq"].upper(), format_line(e, args.include_source))
                for e in new_entries]
    new_lines = sort_rows(new_rows)
    if main_path:
        main_out = [(natural_key(row_name(c)), row_seq_key(c), "\t".join(c)) for c in main_rows]
        body = ([r[2] for r in main_out] + new_lines) if args.main_order else sort_rows(main_out + new_rows)
    else:
        body = new_lines
    header_lines = []
    if main_header:
        header_lines = ["\t".join(main_header)]
    elif WRITE_HEADER:
        header_lines = ["\t".join(HEADER)]
    lines = header_lines + body

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding=args.encoding, newline="") as fh:
        fh.write("\n".join(lines) + "\n")
    if main_path:
        new_only_path.parent.mkdir(parents=True, exist_ok=True)
        with open(new_only_path, "w", encoding=args.encoding, newline="") as fh:
            fh.write("\n".join(header_lines + new_lines) + "\n")

    # ---------------------------------------------------------------- 6. report
    by_name = defaultdict(set)
    for e in final:
        by_name[name_key(e["name"])].add(e["seq"])
    same_name_diff_seq = {n: s for n, s in by_name.items() if len(s) > 1}

    status_counts = defaultdict(int)
    for e in final:
        status_counts[e["check"]] += 1

    print(f"Searched: {root}")
    print(f"CSV files used: {used_files}   skipped: {len(skipped)}")
    for path, problem in skipped:
        print(f"   skipped {path.relative_to(root)}: {problem}")
    print(f"Rows read: {len(entries)}   exact duplicates removed (same name and sequence): {n_same_name}   "
          f"renamed duplicates merged (same sequence, different name): {n_renamed}   "
          f"final entries: {len(final)}")
    print("Status of final entries: " + ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items())))
    if NAME_OVERRIDES:
        used = [f"'{k}' -> '{NAME_OVERRIDES[k]}' ({n} row{'s' if n != 1 else ''})" for k, n in override_hits.items()]
        unused = [f"'{k}'" for k in NAME_OVERRIDES if k not in override_hits]
        if used:
            print("Name overrides applied: " + "; ".join(used))
        if unused:
            print("WARNING, name overrides that matched no row: " + ", ".join(unused))
    if conflicts:
        print(f"\nDuplicate copies that DISAGREE ({len(conflicts)}), flagged in the comments:")
        for c in conflicts:
            print("   " + c)
    if merged_groups:
        print(f"\nSame sequence under different names, merged into one entry each ({len(merged_groups)}). "
              "Each line shows that entry's own sequence, its name and where it came from:")
        for keep, alts in sorted(merged_groups, key=lambda g: natural_key(g[0]["name"])):
            print()
            page = f" p.{keep['page']}" if keep["page"] else ""
            print(f"   {keep['seq']}   {display_name(keep)}   ({keep['source']}{page})   <- kept")
            for c in alts:
                page = f" p.{c['page']}" if c["page"] else ""
                print(f"   {c['seq']}   {display_name(c)}   ({c['source']}{page})   merged into the line above")
    if same_name_diff_seq:
        print(f"\nSame name, different sequences ({len(same_name_diff_seq)}), all kept:")
        for n in sorted(same_name_diff_seq, key=natural_key):
            print(f"   {n}: {len(same_name_diff_seq[n])} different sequences")

    if main_path:
        print(f"\n==== Main file: {main_path}")
        print(f"Main-file rows: {len(main_rows)}" + ("   (header row found)" if main_header else "   (no header row)"))
        dup_in_main = {k: v for k, v in main_seq_rows.items() if len(v) > 1}
        if dup_in_main:
            print(f"WARNING, same sequence on several main-file rows ({len(dup_in_main)}); left as they are:")
            for k, rows_ in dup_in_main.items():
                shown = ", ".join(f"row {sheet_row(i)} '{row_name(main_rows[i])}'" for i in rows_)
                print(f"   {k[:40]}{'...' if len(k) > 40 else ''}: {shown}")
        name_diff = {n: rows_ for n, rows_ in main_name_rows.items()
                     if len({row_seq_key(main_rows[i]) for i in rows_}) > 1}
        if name_diff:
            print(f"WARNING, same name with different sequences in the main file ({len(name_diff)}); "
                  "left as they are:")
            for n, rows_ in name_diff.items():
                print(f"   '{row_name(main_rows[rows_[0]])}': rows " + ", ".join(str(sheet_row(i)) for i in rows_))
        print(f"Spec-sheet entries already in the main file (ignored, main file wins): {len(already_in_main)}")
        renamed = []
        for e, rows_ in already_in_main:
            main_names = {name_key(row_name(main_rows[i])) for i in rows_}
            if not any(name_key(s) in main_names for s in e["sheet_names"]):
                renamed.append((e, rows_))
        if renamed:
            print(f"   of these, the main file uses a different name from the spec sheets ({len(renamed)}); "
                  "main-file name kept:")
            for e, rows_ in renamed:
                mnames = " | ".join(row_name(main_rows[i]) for i in rows_)
                print(f"      row {sheet_row(rows_[0])}: main file '{mnames}'   spec sheets: {' | '.join(e['sheet_names'])}")
        print(f"NEW entries added: {len(new_entries)}")
        if near_matches:
            print(f"   possible duplicates, same bases but different modification codes ({len(near_matches)}), flagged:")
            for e, rows_ in near_matches:
                print(f"      new '{e['name']}' {e['seq']}")
                for i in rows_:
                    print(f"      main-file row {sheet_row(i)} '{row_name(main_rows[i])}' {main_rows[i][0]}")
        if name_clashes:
            print(f"   names already in the main file with a different sequence ({len(name_clashes)}), flagged:")
            for e, rows_ in name_clashes:
                print(f"      new '{e['name']}' {e['seq']}   vs main-file row(s) "
                      + ", ".join(str(sheet_row(i)) for i in rows_))
        print(f"\nWrote {len(body)} entries ({len(main_rows)} from the main file + {len(new_entries)} new) to {out_path}")
        print(f"Wrote {len(new_entries)} new entries only to {new_only_path}")
    else:
        print(f"\nWrote {len(final)} entries to {out_path}")
    if args.echo:
        print()
        print("\n".join(lines))


if __name__ == "__main__":
    main()