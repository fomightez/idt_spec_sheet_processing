#!/usr/bin/env python3
"""
Cross-merge the comments of two lists that describe the same primers, and write new copies of both.

    1. the spec-sheet list: a TSV  (Sequence <TAB> Name <TAB> Comments), e.g. made by the combine script
    2. your own notes list: a TXT  (Sequence <TAB> Name <TAB> Comments), possibly with other text mixed in

For every primer that is in BOTH lists:
    - the notes TXT comment is APPENDED to the end of the TSV comment            (TSV -> BETTER_<tsv name>)
    - the length and Tm from the TSV comment are added to the FRONT of the TXT comment   (TXT -> BETTER_<txt name>)

        TSV:  20 nts; Tm= 51.9 °C; <comment from the TXT>
        TXT:  20 nts; Tm= 51.9 °C; <the original TXT comment>

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
    running the script again on its own output changes nothing.

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
# ====================================================================================

SEQUENCE_CELL = re.compile(r"(?:/[^/\s]+/|[ACGTUNRYSWKMBDHV])+")
LEADING_INFO = re.compile(r"^\s*(?:\d+|\?)\s*nts\s*;\s*Tm=\s*(?:[\d.]+\s*\S*C|n/a)\s*(?:;\s*)?")
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
        description="Append TXT comments to the TSV and add length/Tm from the TSV to the front of the TXT "
                    "comments, for primers that are in both lists. Writes new BETTER_ copies.")
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

        # TSV: append the notes to the end of the comment
        notes = [n for n in (strip_leading_info(r["comment"]) for r in n_rows) if n]
        for r in t_rows:
            current = r["comment"]
            to_add = [n for n in notes if n not in current]          # never add the same note twice
            if not to_add:
                tsv_already += bool(notes)
                continue
            r["new_comment"] = (current.rstrip() + COMMENT_JOINER if current.strip() else "") + COMMENT_JOINER.join(to_add)

        # TXT: put length and Tm from the TSV at the front of the comment
        info = next((x for x in (length_and_tm(r["comment"]) for r in t_rows) if x), None)
        for r in n_rows:
            if info is None:
                no_info.append(r["name"])
            elif r["comment"].lstrip().startswith(info):
                txt_already += 1
            else:
                r["new_comment"] = info + (COMMENT_JOINER + r["comment"] if r["comment"].strip() else "")

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
