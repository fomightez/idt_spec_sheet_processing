#!/usr/bin/env python3
"""
Cross-merge the comments of two lists that describe the same primers, and write new copies of both.

    1. the spec-sheet list: a TSV  (Sequence <TAB> Name <TAB> Comments), e.g. made by the combine script
    2. your own notes list: a TXT  (Sequence <TAB> Name <TAB> Comments), possibly with other text mixed in

For every primer that is in BOTH lists:
    - the notes TXT comment is APPENDED to the end of the TSV comment            (TSV -> BETTER_<tsv name>)
    - the length and Tm from the TSV comment are added to the FRONT of the TXT comment, and the other names the
      TSV knows this primer by are added to the END of it as "also listed as: ..."   (TXT -> BETTER_<txt name>)

        TSV (named "18S rRNA F"):   20 nts; Tm= 51.9 °C; also listed as: 18S rRNA F(1); <comment from the TXT>
        TXT (named "18S_rRNA_F"):   20 nts; Tm= 51.9 °C; <the original TXT comment>; also listed as: 18S rRNA F | 18S rRNA F(1)

    The TXT is treated as the main source of names, so its names are never changed. The list added to a TXT
    comment holds, in this order:
        1. the name the TSV shows for that primer, if it differs from the TXT line's name
           (can be switched off with ADD_TSV_NAME_TO_ALSO_LISTED in the SETTINGS block), then
        2. the names in the TSV comment's own "also listed as: ..." part, if it has one (the other names the same
           sequence was ordered under, written by the combine script).
    The TXT line's OWN name is always left out, because it is not an "also" for that line, and so is any name
    the TXT comment already lists. If no names are left, nothing is added. Names are compared ignoring
    upper/lower case and extra spaces.

    The "also listed as: ..." part of a TSV comment is everything from those words up to the next "; " or the end
    of the comment. The names in it are separated by " | " (a space, a vertical bar and a space), the way the
    combine script writes them.

    Example: the TSV shows the name IE1_Rev] (an OCR misread) with "also listed as: IE1_Rev1 | IE1-72 Reverse",
    and the TXT line is named IE1_Rev1. Its TXT comment gets "also listed as: IE1_Rev] | IE1-72 Reverse".

    (Note that when this was originally developed, it was intended for ones that had no `also listed as: ` 
    occurences because they were for my own primers that had no duplicates among them and I was just wanting
    to capture & add the length and Tm details the spec sheets had to the versions of my primers that would 
    kick off the lst to become the main, established list of lab primers.; however, I realized later I may want
    to use this script to reconcile cases where there are `also listed as: ` items & so I added this handling
    towards the end and I haven't adequately tested whether the TXT file will be respected as the main version
    for the name when considering alternate names listed.)

Usage:
    python merge_primerDetails_n_comments.py extracted_info.tsv my_primers.txt
    python merge_primerDetails_n_comments.py extracted_info.tsv my_primers.txt --out-dir merged/
    python merge_primerDetails_n_comments.py extracted_info.tsv my_primers.txt --tabs

How primers are matched
    By SEQUENCE (ignoring upper/lower case and spaces), not by name, because sequence is the key to the primer working. Modification codes are part of the sequence, so /56-FAM/... has to be
    spelled the same in both lists. Primers whose bases agree but whose modification codes differ are NOT
    merged; they are listed in the console report as possible matches so you can check them.

What is left alone
    The input files are never changed; new files with the prefix BETTER_ are written. Every other line of each
    file is copied exactly as it was: blank lines, notes, headers, the Upstate list of "name <TAB> 5' - SEQ - 3'"
    lines, line endings and the file's text encoding. Lines that are already merged are not merged twice, so
    running the script again on its own output changes nothing. (The "also listed as: ..." text that this script
    put at the end of a TXT comment is not copied back into the TSV, which already has those names.) BETTER_ files
    made by an earlier version of this script should be remade from your original files, not fed back in.

Lines separated by spaces
    A primer line in the TXT may be separated by runs of two or more spaces instead of tabs. These are read the
    same way and, by default, keep their spaces. With --tabs they are converted to tabs in the new TXT, so that
    the whole file really is tab-separated.

Settings you may want to change are in the SETTINGS block just below the imports.
"""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

# ================================================================ SETTINGS
OUTPUT_PREFIX = "BETTER_"      # put in front of each input file's name to make the new file's name
COMMENT_JOINER = "; "          # between the pieces of a comment, as in "21 nts; Tm= 50.2 °C; any comments"
MIN_SEQUENCE_LENGTH = 8        # a first cell must have at least this many characters to count as a sequence
ALSO_LISTED_LABEL = "also listed as:"   # the words the combine script puts in front of the other names
ALSO_SEPARATOR = " | "                  # what the combine script puts between those names
ADD_TSV_NAME_TO_ALSO_LISTED = True      # list the TSV's name for a primer among the "also listed as" names added to
                                        # its TXT comment when the TXT line is named differently (False = only the
                                        # names from the TSV comment's own "also listed as: ..." part)
# ====================================================================================

SEQUENCE_CELL = re.compile(r"(?:/[^/\s]+/|[ACGTUNRYSWKMBDHV])+")
LEADING_INFO = re.compile(r"^\s*(?:\d+|\?)\s*nts\s*;\s*Tm=\s*(?:[\d.]+\s*\S*C|n/a)\s*(?:;\s*)?")
# "also listed as: A | B" up to the next "; " (or the end of the comment)
ALSO_LISTED = re.compile(re.escape(ALSO_LISTED_LABEL) + r"(?:(?!;\s).)*", re.IGNORECASE)
# an "also listed as: ..." part that runs to the END of a note (so it can be taken off again)
TRAILING_ALSO = re.compile(r"(?:;\s*)?" + re.escape(ALSO_LISTED_LABEL) + r"(?:(?!;\s).)*$", re.IGNORECASE)
LENGTH_AND_TM = re.compile(r"(?P<length>(?:\d+|\?)\s*nts)\s*;\s*(?P<tm>Tm=\s*(?:[\d.]+\s*\S*C|n/a))")


# ---------------------------------------------------------------- small helpers
def seq_key(text):
    """Comparison form of a sequence: no spaces, upper case."""
    return re.sub(r"\s+", "", text).upper()


def bare_key(text):
    """The bases only: modification codes removed."""
    return re.sub(r"/[^/]*/", "", seq_key(text))


def shown(name, seq=""):
    """A name for the console report; rows with no name are shown by the start of their sequence."""
    return name if name else f"(no name: {seq[:18]}...)"


def name_key(text):
    return " ".join(text.casefold().split())


def looks_like_sequence(cell):
    key = seq_key(cell)
    return len(key) >= MIN_SEQUENCE_LENGTH and SEQUENCE_CELL.fullmatch(key) is not None


def strip_leading_info(comment):
    """A TXT note without the 'N nts; Tm= X °C' this script may have put at its front. That keeps the length and
    Tm from being copied back into the TSV, and lets a second run recognise notes it has already appended."""
    return LEADING_INFO.sub("", comment, count=1).strip()


def length_and_tm(comment):
    """'20 nts; Tm= 51.9 °C' taken from a TSV comment, or None if it isn't there."""
    m = LENGTH_AND_TM.search(comment)
    if not m:
        return None
    return f"{' '.join(m.group('length').split())}{COMMENT_JOINER.rstrip()} {' '.join(m.group('tm').split())}"


def also_listed_part(comment):
    """The 'also listed as: A | B' part of a TSV comment, or None if it has none."""
    m = ALSO_LISTED.search(comment)
    return m.group(0).strip() if m else None


def split_also_names(part):
    """['A', 'B'] from 'also listed as: A | B'. The split is on ' | ' as a whole, not on '|', because a name can
    itself start with a pipe (an OCR misread of an I, such as '|E1-72 Forward')."""
    part = part.lstrip("; \t")
    if part.lower().startswith(ALSO_LISTED_LABEL.lower()):
        part = part[len(ALSO_LISTED_LABEL):]
    return [n.strip() for n in part.strip().split(ALSO_SEPARATOR) if n.strip()]


def also_listed_names(comment):
    """Every name in every 'also listed as: ...' part found in a comment."""
    names = []
    for m in ALSO_LISTED.finditer(comment):
        names += split_also_names(m.group(0))
    return names


def strip_also_listed(note, known_name_keys):
    """A TXT note without an 'also listed as: ...' part at its END, when every name in that part is one the TSV
    already lists (so it is text this script added, which must not be copied back into the TSV). A part that
    holds any name the TSV does not have is your own writing and is left alone."""
    note = note.rstrip()
    m = TRAILING_ALSO.search(note)
    if m:
        names = split_also_names(m.group(0))
        if names and all(name_key(n) in known_name_keys for n in names):
            note = note[:m.start()].rstrip()
    return note.strip()


# ---------------------------------------------------------------- reading and writing, keeping encodings
def read_text(path):
    """Return (lines, encoding) where lines is a list of [text_without_line_ending, line_ending]."""
    raw = Path(path).read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        encoding = "utf-8-sig"
    elif raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        encoding = "utf-16"
    else:
        try:
            raw.decode("utf-8")
            encoding = "utf-8"
        except UnicodeDecodeError:
            encoding = "cp1252"
    parts = raw.decode(encoding).split("\n")
    lines = []
    for i, part in enumerate(parts):
        last = i == len(parts) - 1
        if last and part == "":
            break
        eol = "" if last else "\n"
        if not last and part.endswith("\r"):
            part, eol = part[:-1], "\r\n"
        lines.append([part, eol])
    return lines, encoding


def write_text(path, lines, encoding):
    text = "".join(body + eol for body, eol in lines)
    try:
        data = text.encode(encoding)
    except UnicodeEncodeError:
        print(f"   note: {Path(path).name} could not be written as {encoding}; using utf-8 instead")
        data = text.encode("utf-8")
    Path(path).write_bytes(data)


# ---------------------------------------------------------------- parsing primer lines
def split_cells(body, allow_space_separators):
    """Split a line into (cells, separators). Tabs separate cells; if allowed, so do runs of 2+ spaces
    (only the first two separators count, because the comment itself contains single spaces)."""
    if "\t" in body:
        cells = body.split("\t")
        return cells, ["\t"] * (len(cells) - 1)
    if allow_space_separators:
        parts = re.split(r"( {2,})", body, maxsplit=2)
        return parts[0::2], parts[1::2]
    return [body], []


def is_header(cells):
    return len(cells) >= 2 and cells[0].strip().casefold() == "sequence" and cells[1].strip().casefold() == "name"


def parse_primer_lines(lines, free_text):
    """Return the lines that are primers.
    free_text=True  (the notes TXT, which has other text mixed in): the first cell must look like a sequence,
                    and spaces (2 or more in a row) may separate the cells.
    free_text=False (the TSV): any tab-separated row with something in its first cell, except the header,
                    so a row whose sequence was misread by OCR (e.g. an 'I' in it) still counts as a primer."""
    rows = []
    for i, (body, _) in enumerate(lines):
        cells, seps = split_cells(body, free_text)
        if free_text:
            is_primer = len(cells) >= 2 and looks_like_sequence(cells[0])
        else:
            is_primer = len(cells) >= 2 and cells[0].strip() != "" and not is_header(cells)
        if is_primer:
            rows.append({
                "line": i, "cells": cells, "seps": seps,
                "seq": seq_key(cells[0]), "name": cells[1].strip(),
                "comment": cells[2] if len(cells) > 2 else "",
                "space_separated": "\t" not in body,
            })
    return rows


def index_by(rows, key):
    out = defaultdict(list)
    for r in rows:
        out[key(r["seq"])].append(r)
    return out


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        description="Append TXT comments to the TSV; add length/Tm from the TSV to the front of the TXT comments and "
                    "any 'also listed as: ...' part to the end, for primers that are in both lists. "
                    "Writes new BETTER_ copies.")
    ap.add_argument("tsv", help="the spec-sheet TSV: Sequence, Name, Comments")
    ap.add_argument("txt", help="your notes TXT: Sequence, Name, Comments (other text allowed)")
    ap.add_argument("--out-dir", help="folder for the new files (default: next to each input file)")
    ap.add_argument("--prefix", default=OUTPUT_PREFIX, help=f"prefix for the new file names (default {OUTPUT_PREFIX})")
    ap.add_argument("--tabs", action="store_true",
                    help="in the new TXT, convert primer lines separated by spaces into tab-separated lines")
    args = ap.parse_args()

    tsv_path, txt_path = Path(args.tsv).expanduser(), Path(args.txt).expanduser()
    for p in (tsv_path, txt_path):
        if not p.is_file():
            sys.exit(f"File not found: {p}")
    out_dir = Path(args.out_dir).expanduser() if args.out_dir else None
    tsv_out = (out_dir or tsv_path.parent) / (args.prefix + tsv_path.name)
    txt_out = (out_dir or txt_path.parent) / (args.prefix + txt_path.name)
    for src, dst in ((tsv_path, tsv_out), (txt_path, txt_out)):
        if src.resolve() == dst.resolve():
            sys.exit(f"The new file would overwrite {src}. Use a non-empty --prefix or a different --out-dir.")

    tsv_lines, tsv_enc = read_text(tsv_path)
    txt_lines, txt_enc = read_text(txt_path)
    tsv_rows = parse_primer_lines(tsv_lines, free_text=False)
    txt_rows = parse_primer_lines(txt_lines, free_text=True)
    tsv_by_seq, txt_by_seq = index_by(tsv_rows, lambda s: s), index_by(txt_rows, lambda s: s)

    # ------------------------------------------------------------ merge, primer by primer (matched on sequence)
    matched = [s for s in tsv_by_seq if s in txt_by_seq]
    name_differs, no_info, tsv_already, txt_already = [], [], 0, 0
    for seq in matched:
        t_rows, n_rows = tsv_by_seq[seq], txt_by_seq[seq]

        # every name the TSV knows this primer by: first the name it shows, then the "also listed as: ..." names
        # (taken from the FIRST such part of each comment, the one the combine script wrote; anything after it is
        # text that was appended later)
        t_alsos, t_names, t_name_keys = [], [], set()
        main_keys = {name_key(r["name"]) for r in t_rows if r["name"]}

        def add_name(n):
            if name_key(n) not in t_name_keys:
                t_name_keys.add(name_key(n))
                t_names.append(n)

        for r in t_rows:
            if r["name"]:
                add_name(r["name"])
        for r in t_rows:
            part = also_listed_part(r["comment"])
            if part and part not in t_alsos:
                t_alsos.append(part)
            for n in (split_also_names(part) if part else []):
                add_name(n)

        # TSV: append the notes to the end of the comment
        notes = [n for n in (strip_also_listed(strip_leading_info(r["comment"]), t_name_keys) for r in n_rows) if n]
        for r in t_rows:
            current = r["comment"]
            to_add = [n for n in notes if n not in current]          # never add the same note twice
            if not to_add:
                tsv_already += bool(notes)
                continue
            r["new_comment"] = (current.rstrip() + COMMENT_JOINER if current.strip() else "") + COMMENT_JOINER.join(to_add)

        # TXT: put length and Tm from the TSV at the front of the comment, and the "also listed as: ..." part at the end
        info = next((x for x in (length_and_tm(r["comment"]) for r in t_rows) if x), None)
        for r in n_rows:
            if info is None:
                no_info.append(r["name"])
            comment = r["comment"]
            if info and not comment.lstrip().startswith(info):
                comment = info + (COMMENT_JOINER + comment if comment.strip() else "")
            added_also = False
            # the other names for THIS line: not its own name, and not any the comment already lists
            already = {name_key(n) for n in also_listed_names(comment)} | {name_key(r["name"])}
            missing = [n for n in t_names
                       if name_key(n) not in already
                       and (ADD_TSV_NAME_TO_ALSO_LISTED or name_key(n) not in main_keys)]
            if missing:
                comment = ((comment.rstrip() + COMMENT_JOINER) if comment.strip() else "") \
                          + ALSO_LISTED_LABEL + " " + ALSO_SEPARATOR.join(missing)
                added_also = True
            if comment != r["comment"]:
                r["new_comment"], r["also_added"] = comment, added_also
            elif info or t_alsos:
                txt_already += 1

        if name_key(t_rows[0]["name"]) != name_key(n_rows[0]["name"]):
            name_differs.append((n_rows[0]["name"], t_rows[0]["name"]))

    # ------------------------------------------------------------ rebuild the lines
    for r in tsv_rows:
        if "new_comment" in r:
            cells = r["cells"] + [""] * (3 - len(r["cells"]))
            cells[2] = r["new_comment"]
            tsv_lines[r["line"]][0] = "\t".join(cells)

    converted = []
    for r in txt_rows:
        cells, seps = list(r["cells"]), list(r["seps"])
        if "new_comment" in r:
            if len(cells) < 3:
                cells.append("")
                seps.append(seps[-1] if seps else "\t")
            cells[2] = r["new_comment"]
        if args.tabs and r["space_separated"]:
            seps = ["\t"] * len(seps)
            converted.append(r["name"])
        if "new_comment" in r or (args.tabs and r["space_separated"]):
            txt_lines[r["line"]][0] = "".join(c + s for c, s in zip(cells, seps + [""]))

    tsv_out.parent.mkdir(parents=True, exist_ok=True)
    txt_out.parent.mkdir(parents=True, exist_ok=True)
    write_text(tsv_out, tsv_lines, tsv_enc)
    write_text(txt_out, txt_lines, txt_enc)

    # ------------------------------------------------------------ report
    n_tsv_changed = sum("new_comment" in r for r in tsv_rows)
    n_txt_changed = sum("new_comment" in r for r in txt_rows)
    only_tsv = [r for s, rs in tsv_by_seq.items() if s not in txt_by_seq for r in rs]
    only_txt = [r for s, rs in txt_by_seq.items() if s not in tsv_by_seq for r in rs]
    space_lines = [r for r in txt_rows if r["space_separated"]]

    print(f"TSV: {tsv_path}   ({len(tsv_rows)} primer rows, encoding {tsv_enc})")
    print(f"TXT: {txt_path}   ({len(txt_rows)} primer lines, {len(txt_lines)} lines in all, encoding {txt_enc})")
    print(f"\nPrimers in both lists (matched by sequence): {len(matched)}")
    print(f"   TSV comments extended with the notes:     {n_tsv_changed}" + (f"   ({tsv_already} already had them)" if tsv_already else ""))
    print(f"   TXT comments given length and Tm in front: {n_txt_changed}" + (f"   ({txt_already} already had them)" if txt_already else ""))
    print(f"   of those, given '{ALSO_LISTED_LABEL} ...' at the end:  {sum(bool(r.get('also_added')) for r in txt_rows)}")
    if no_info:
        print(f"   WARNING, no length/Tm found in the TSV comment for: {', '.join(no_info)}")

    if name_differs:
        print(f"\nSame sequence, different name ({len(name_differs)}); each file keeps its own name:")
        for n_name, t_name in name_differs:
            print(f"   TXT '{n_name}'   TSV '{t_name}'")
    if space_lines:
        what = "converted to tabs in the new TXT" if args.tabs else "left as they are (use --tabs to convert them)"
        print(f"\n{len(space_lines)} TXT primer line(s) are separated by spaces, not tabs; {what}:")
        print("   " + ", ".join(r["name"] for r in space_lines))

    tsv_by_bare = defaultdict(list)
    for r in only_tsv:
        tsv_by_bare[bare_key(r["seq"])].append(r)
    possible = [(r, tsv_by_bare[bare_key(r["seq"])]) for r in only_txt if bare_key(r["seq"]) in tsv_by_bare]
    if possible:
        print(f"\nPossible matches NOT merged, same bases but different modification codes ({len(possible)}):")
        for n_row, t_rows in possible:
            print(f"   TXT '{n_row['name']}' {n_row['cells'][0].strip()}")
            for t in t_rows:
                print(f"   TSV '{t['name']}' {t['cells'][0].strip()}")
    shown_txt = [r for r in only_txt if not any(r is p[0] for p in possible)]
    if shown_txt:
        print(f"\nIn the TXT only ({len(shown_txt)}), left unchanged: " + ", ".join(shown(r["name"], r["seq"]) for r in shown_txt))
    print(f"In the TSV only ({len(only_tsv)}), left unchanged"
          + (f": {', '.join(shown(r['name'], r['seq']) for r in only_tsv[:12])}" + (" ..." if len(only_tsv) > 12 else "") if only_tsv else ""))

    n_other_tsv = len(tsv_lines) - len(tsv_rows)
    print(f"\nOther lines copied exactly as they were: {n_other_tsv} in the TSV (header and anything without a "
          f"name), {len(txt_lines) - len(txt_rows)} in the TXT (notes, blank lines, other lists)")
    print(f"\nWrote {tsv_out}")
    print(f"Wrote {txt_out}")


if __name__ == "__main__":
    main()