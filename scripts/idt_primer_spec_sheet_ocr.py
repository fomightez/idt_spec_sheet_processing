#!/usr/bin/env python3
"""
Extract name, sequence, number of bases and Tm from scanned IDT oligo spec sheets.

Setup:
    pip install pdf2image pytesseract pillow
    # system packages:  tesseract-ocr  and  poppler-utils
    #   macOS:   brew install tesseract poppler
    #   Ubuntu:  sudo apt install tesseract-ocr poppler-utils
    #   Windows: install Tesseract (UB Mannheim build) + poppler, add both to PATH

Usage:
    python idt_primer_spec_sheet_ocr.py my_specs.pdf                 # writes my_specs.csv
    python idt_primer_spec_sheet_ocr.py my_specs.pdf --debug         # also dumps raw OCR text per page

Text clean-up settings (SETTINGS block, just below the imports)
---------------------------------------------------------------
Two replacement rounds run on the OCR text BEFORE any parsing or checks, to handle known edge cases.

Round 1 - known_text_error_substitutions
    A dictionary of exact, case-sensitive replacements applied to the whole page text.
    Add entries as you find recurring OCR errors, e.g. {'TB4AQE': 'TB40E'}.

Round 2 - disregard_dashes_in_sequence
    Ignores stray pencil/pen marks that OCR reads as a dash inside the sequence (the text between
    "5'-" and "-3'"). Accepted values:
        False                                 never remove dashes (default)
        'all'                                 every page of every file
        ('0906_', [9])                        only files whose NAME CONTAINS '0906_', only page 9
                                              (PDF page numbers, same as the CSV 'page' column)
        [('0906_', [9]), ('0907_', [1, 2])]   several rules; a page matching any rule is fixed
    Notes:
      - Dashes inside modification codes such as /56-FAM/ are left alone.
      - En and em dashes are treated as dashes too.
      - The rule also applies when a sequence is re-read from a crop during the retry step.
      - Any other value (e.g. 'All' or True) raises an error at startup.
      - Marks that touch letters and change them (T read as F, say) are NOT fixed by this setting;
        those rows still show up as REVIEW.

Audit trail: what each round changed is recorded in the CSV 'text_fixes' column, for example
'TB4AQE->TB40E' or 'removed 2 dash(es) in sequence'. The --debug text files keep the raw Tesseract
output from before either round.
"""
import argparse
import csv
import re
from pathlib import Path

import pytesseract
from pdf2image import convert_from_path, pdfinfo_from_path
from PIL import ImageOps

# ================================================================ SETTINGS
# Two text clean-up rounds run on the OCR text BEFORE any parsing or checks (see apply_text_fixes()).

# Round 1: exact string replacements for OCR errors you have seen repeatedly (case-sensitive,
# applied to the whole page text). Add entries as you find patterns.
known_text_error_substitutions = {
    'TB4AQE': 'TB40E',
    'TBAOE': 'TB40E',
}

# Round 2: ignore stray pencil/pen marks that OCR reads as a dash inside the sequence
# (the text between "5'-" and "-3'"). Dashes inside modification codes such as /56-FAM/ are left alone.
# Options:
#   False                           -> never remove dashes 
#   'all'                           -> remove them on every page of every file
#   ('0906_', [9])                  -> only on files whose NAME CONTAINS '0906_', and only on page 9
#                                      (page numbers are PDF page numbers, same as the CSV 'page' column)
#                                      (page numbers per document can be a list of each separated by comma)
#   [('0905_', [2]),('0906_', [9,18]) ]   -> several of the above rules (a page matching any rule is fixed)
#disregard_dashes_in_sequence = [('0905_', [2]),('0906_', [9]) ] 
disregard_dashes_in_sequence = 'all'
# ====================================================================================

# ---------------------------------------------------------------- OCR step
# Whole-page layout mode for the main OCR pass. psm 6 ("one uniform block") splits rows on sheets from 2015 on,
# because of the extra "For 100 uM: add X uL" text in the middle column; psm 4 ("single column, variable
# sizes") reads them correctly and handles the older sheets too.
PAGE_PSM = 4

KEYWORDS = ("SPECIFICATION", "Sequence", "Oligo", "bases", "Tm", "Properties", "Shipped")

def _score(text):
    """How many expected spec-sheet words appear? Used to pick the right orientation."""
    return sum(text.count(k) for k in KEYWORDS)

def ocr_image(img):
    """OCR one page image, trying 0/180 degrees first (scans are often upside down),
    then 90/270 only if neither looks like a spec sheet.
    Returns (text, oriented_image) so later retry steps can reuse the right-way-up image."""
    img = ImageOps.autocontrast(img.convert("L"))
    best_text, best_img, best_score = "", img, -1
    for angles in ((0, 180), (90, 270)):
        for angle in angles:
            rotated = img.rotate(angle, expand=True) if angle else img
            text = pytesseract.image_to_string(rotated, config=f"--psm {PAGE_PSM}")
            sc = _score(text)
            if sc > best_score:
                best_text, best_img, best_score = text, rotated, sc
        if best_score >= 6:   # clearly a spec sheet; don't bother with sideways
            break
    return best_text, best_img

# ----------------------------------------------------------- parsing step
# OCR often turns the prime mark into ' ’ ` or ′, and sometimes 5 -> S.
PRIME = r"[\'\"’`′´]"
SEQ_RE = re.compile(
    rf"(?<![A-Za-z0-9])[5S]\s*{PRIME}+\s*[-–—=~._:]\s*(?P<seq>.+?)\s*-\s*3\s*{PRIME}*",
    re.DOTALL,
)
# Layout seen in real IDT sheets (OCR'd text):
#   SPECIFICATION SHEET www.idtdna.com
#   Sequence - v2_Plac_T
#   25 nmole DNA Oligo, 60 bases
#   Tm (50mM NaCl)*: 73.6 °C          <- the Tm we want
#   Strongest Folding Tm: 44.4 °C     <- a different Tm; captured separately
# Note OCR sometimes misreads "50mM" as "S0mM", so the parenthetical is matched loosely.
SHEET_SPLIT_RE = re.compile(r"SPECIFICATION\s+SHEET", re.IGNORECASE)
NAME_RE = re.compile(r"Sequence\s*[-–—~=:.]+\s*(?P<name>\S[^\n]*)", re.IGNORECASE)
NAME_TAIL_RE = re.compile(r"\s+\d+(?:\.\d+)?\s*(?:nmoles?|nm|umoles?|µmoles?|μmoles?|OD)\b.*$", re.IGNORECASE)
BASES_RE = re.compile(r"(\d+)\s*bases?\b", re.IGNORECASE)
TM_RE = re.compile(r"\bTm\s*\([^)]*\)\s*\*?\s*[:;]?\s*(\d{2}(?:\.\d)?)", re.IGNORECASE)
FOLD_TM_RE = re.compile(r"Folding\s*Tm\s*[:;]?\s*(\d{2}(?:\.\d)?)", re.IGNORECASE)
MW_RE = re.compile(r"Mol\w*\s*Weight\s*[:;]?\s*(\d{1,3}(?:,\d{3})*(?:\.\d+)?)", re.IGNORECASE)

# Molecular-weight check. DNA residue (dNMP) masses; MW = sum(residues) - 61.96 for a 5'-OH/3'-OH oligo.
# On 30 of 30 clean unmodified sheets this agreed with IDT's printed value to within ~0.1 Da.
# Any single-base substitution moves the MW by >= 9 Da, so a 1 Da tolerance is safe.
RESIDUE_MASS = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.20}
MW_TERMINI = 61.96
MW_TOLERANCE = 1.0

# Extinction-coefficient check (260 nm, nearest-neighbour method). Unlike MW, this depends on base ORDER,
# so it catches transposed bases that the MW check can't.
#   e = sum(e[NN pair, i..i+1]) - sum(e[interior base])      (M^-1 cm^-1)
# The table below reproduced IDT's printed value exactly on 31 of 31 unmodified sheets
# (note TpA = 23,400). Printed values are rounded to 100, so allow 100 of slack.
EXT_NN = {"AA": 27400, "AC": 21200, "AG": 25000, "AT": 22800,
          "CA": 21200, "CC": 14600, "CG": 18000, "CT": 15200,
          "GA": 25200, "GC": 17600, "GG": 21600, "GT": 20000,
          "TA": 23400, "TC": 16200, "TG": 19000, "TT": 16800}
EXT_MONO = {"A": 15400, "C": 7400, "G": 11500, "T": 8700}
EXT_TOLERANCE = 100
OK_STATES = ("OK", "OK_LENGTH_ONLY")   # rows that need no retry
EXT_RE = re.compile(r"Coeff\w*\s*[:;]?\s*(\d{1,3}(?:,\d{3})+|\d{4,7})", re.IGNORECASE)

def clean_sequence(raw):
    """Join wrapped lines, drop spaces and stray punctuation (e.g. 'CCA.AGG'),
    uppercase the bases but leave /modification/ codes alone."""
    s = re.sub(r"[^A-Za-z0-9/+*\-]", "", raw)
    parts = re.split(r"(/[^/]*/)", s)            # keep /mod/ segments intact
    return "".join(p if p.startswith("/") else p.upper() for p in parts)

def calc_mw(seq):
    """MW from sequence, or None if it has modifications / non-ACGT letters (masses unknown)."""
    if not seq or "/" in seq or any(c not in RESIDUE_MASS for c in seq):
        return None
    return round(sum(RESIDUE_MASS[c] for c in seq) - MW_TERMINI, 1)

def calc_ext(seq):
    """Nearest-neighbour extinction coefficient, or None if modified / non-ACGT / too short."""
    if len(seq) < 2 or "/" in seq or any(c not in EXT_MONO for c in seq):
        return None
    return (sum(EXT_NN[seq[i:i + 2]] for i in range(len(seq) - 1))
            - sum(EXT_MONO[c] for c in seq[1:-1]))

def evaluate(seq, reported, mw_sheet, ext_sheet):
    """Run all consistency checks on a sequence. Returns the derived columns."""
    bare = re.sub(r"/[^/]*/", "", seq)

    mw_calc = calc_mw(seq)
    if mw_calc is None or mw_sheet is None:
        mw_check = "n/a"                      # modified oligo, or MW not readable on the sheet
    else:
        mw_check = "match" if abs(mw_calc - mw_sheet) <= MW_TOLERANCE else "MISMATCH"

    ext_calc = calc_ext(seq)
    if ext_calc is None or ext_sheet is None:
        ext_check = "n/a"
    else:
        ext_check = "match" if abs(ext_calc - ext_sheet) <= EXT_TOLERANCE else "MISMATCH"

    if not seq:
        check = "NO_SEQ"
    elif re.search(r"[^ACGTUNRYSWKMBDHV]", bare):
        check = "REVIEW"                      # unexpected characters in the sequence
    elif reported is None or reported != len(bare):
        check = "REVIEW"                      # base count missing or doesn't match
    elif "MISMATCH" in (mw_check, ext_check):
        check = "REVIEW"                      # right length, but composition/order disagrees with the sheet
    elif mw_check == "n/a" and ext_check == "n/a":
        # Typically a modified oligo (dyes/quenchers): only the base count could be checked, so a
        # same-length misread (e.g. A->G) would NOT be caught. Look these over by eye.
        check = "OK_LENGTH_ONLY"
    else:
        check = "OK"
    return {"bases_counted": len(bare) if seq else None, "mw_calc": mw_calc, "mw_check": mw_check,
            "ext_calc": ext_calc, "ext_check": ext_check, "check": check}

def parse_sheet(text):
    """Parse one spec sheet's OCR text into a dict."""
    nm = NAME_RE.search(text)
    name = NAME_TAIL_RE.sub("", nm.group(1)).strip() if nm else ""

    sm = SEQ_RE.search(text)
    seq = clean_sequence(sm.group("seq")) if sm else ""

    tm = TM_RE.search(text)
    fold = FOLD_TM_RE.search(text)
    bm = BASES_RE.search(text)
    mm = MW_RE.search(text)
    em = EXT_RE.search(text)
    bases_reported = int(bm.group(1)) if bm else None
    mw_sheet = float(mm.group(1).replace(",", "")) if mm else None
    ext_sheet = int(em.group(1).replace(",", "")) if em else None

    return {
        "name": name,
        "sequence": seq,
        "bases_reported": bases_reported,
        "tm_c": float(tm.group(1)) if tm else None,
        "folding_tm_c": float(fold.group(1)) if fold else None,
        "mw_sheet": mw_sheet,
        "ext_sheet": ext_sheet,
        **evaluate(seq, bases_reported, mw_sheet, ext_sheet),
    }

# ------------------------------------------------------- text clean-up rounds (run before any checks)
DASH_CHARS = "-–—"      # hyphen, en dash, em dash

def _dash_rules():
    """Validate the disregard_dashes_in_sequence setting and return it as a list of (tag, pages) rules,
    or the string 'all', or [] when switched off."""
    s = disregard_dashes_in_sequence
    if s is False or s is None:
        return []
    if s == 'all':
        return 'all'
    def is_rule(r):
        return (isinstance(r, tuple) and len(r) == 2 and isinstance(r[0], str)
                and isinstance(r[1], (list, tuple, set)))
    if is_rule(s):
        return [s]
    if isinstance(s, list) and s and all(is_rule(r) for r in s):
        return s
    raise ValueError("disregard_dashes_in_sequence must be False, 'all', a tuple like ('0906_', [9]), "
                     f"or a list of such tuples; got {s!r}")

def dashes_enabled(file_name, page):
    """Should stray dashes in the sequence be ignored for this file name and page number?"""
    rules = _dash_rules()
    if rules == 'all':
        return True
    return any(tag in file_name and page in pages for tag, pages in rules)

def strip_dashes_outside_mods(seq_text):
    """Remove dash characters from sequence text, except inside /modification/ codes like /56-FAM/."""
    parts = re.split(r"(/[^/]*/)", seq_text)
    return "".join(p if p.startswith("/") else re.sub(f"[{DASH_CHARS}]", "", p) for p in parts)

def remove_sequence_dashes(text):
    """Round 2: drop stray dashes inside the 5'-...-3' region of OCR text. Returns (new_text, n_removed)."""
    m = SEQ_RE.search(text)
    if not m:
        return text, 0
    old = m.group("seq")
    new = strip_dashes_outside_mods(old)
    removed = sum(old.count(c) for c in DASH_CHARS) - sum(new.count(c) for c in DASH_CHARS)
    return text[:m.start("seq")] + new + text[m.end("seq"):], removed

def apply_text_fixes(text, file_name, page):
    """Both clean-up rounds, in order. Returns (fixed_text, notes) where notes lists what changed."""
    notes = []
    for wrong, right in known_text_error_substitutions.items():          # round 1
        n = text.count(wrong)
        if n:
            text = text.replace(wrong, right)
            notes.append(f"{wrong}->{right}" + (f" x{n}" if n > 1 else ""))
    if dashes_enabled(file_name, page):                                   # round 2
        text, removed = remove_sequence_dashes(text)
        if removed:
            notes.append(f"removed {removed} dash(es) in sequence")
    return text, notes

# ------------------------------------------------------- retry step
# Tesseract's whole-page mode (--psm 6) occasionally mis-reads a character on the sequence line
# (e.g. 'CCG TC' -> 'CCGG TC', or 'ACG' -> 'AGG') even when the scan is clean. The sheet prints the
# base count and the molecular weight, so we use them as checksums: if the extracted sequence fails
# either, crop just the sequence region and re-OCR it with other settings until both agree.
START_TOKEN_RE = re.compile(r"^[5S]['\"’`′´]+[-–—=~._:]")
END_TOKEN_RE = re.compile(r"(?:^|-)3['\"’`′´°.]*$")

def _sequence_box(img):
    """Locate the pixel box around the '5'- ... -3'' region (may span several lines)."""
    d = pytesseract.image_to_data(img, config="--psm 6", output_type=pytesseract.Output.DICT)
    idx = [i for i, t in enumerate(d["text"]) if t.strip()]
    start = next((i for i in idx if START_TOKEN_RE.match(d["text"][i].strip())), None)
    if start is None:
        return None
    end = next((i for i in idx if i >= start and END_TOKEN_RE.search(d["text"][i].strip())), None)
    if end is None:
        return None
    span = [i for i in idx if start <= i <= end]
    pad = 12
    return (max(min(d["left"][i] for i in span) - pad, 0),
            max(min(d["top"][i] for i in span) - pad, 0),
            max(d["left"][i] + d["width"][i] for i in span) + pad,
            max(d["top"][i] + d["height"][i] for i in span) + pad)

def refine_sequence(img, reported, mw_sheet, ext_sheet, strip_dashes=False):
    """Return (sequence, description_of_what_worked) or None if no attempt passes every check."""
    if reported is None:
        return None
    attempts = []
    box = _sequence_box(img)
    if box:
        crop = img.crop(box)
        for scale in (2, 3, 1):
            big = ImageOps.expand(crop.resize((crop.width * scale, crop.height * scale)), border=30, fill=255)
            for psm in (7, 6):
                attempts.append((f"crop x{scale} psm{psm}", big, psm))
    for psm in (6, 11):
        attempts.append((f"full page psm{psm}", img, psm))
    for label, im, psm in attempts:
        candidate = pytesseract.image_to_string(im, config=f"--psm {psm}")
        if strip_dashes:
            candidate, _ = remove_sequence_dashes(candidate)
        m = SEQ_RE.search(candidate)
        if m:
            seq = clean_sequence(m.group("seq"))
            if evaluate(seq, reported, mw_sheet, ext_sheet)["check"] in OK_STATES:
                return seq, label
    return None

# --------------------------------------------- re-reading the printed MW / extinction coefficient
# Newer sheets have extra text ("For 100 uM: add 273 uL") in the right-hand column that sits between two
# rows of the left column. Tesseract's whole-page mode groups words into lines by height, so it can split
# the "Molecular Weight: 6,778.5" row into pieces and drop part of the number. When a printed value fails its
# check we re-read just that row from a tight crop (left column only).
MW_CROP_RE = re.compile(r"Weight\s*[:;]?\s*(\d{1,3}(?:,\d{3})*(?:\.\d+)?)", re.IGNORECASE)
FIELD_ANCHORS = {"mw": (re.compile(r"^Weight", re.IGNORECASE), MW_CROP_RE),
                 "ext": (re.compile(r"^Coeff", re.IGNORECASE), EXT_RE)}

def reread_field(img, field, data):
    """Candidate values for a printed field, re-OCR'd from a one-row crop around its label."""
    anchor_re, value_re = FIELD_ANCHORS[field]
    found = []
    for i, t in enumerate(data["text"]):
        if not anchor_re.match(t.strip()):
            continue
        top, h, left = data["top"][i], data["height"][i], data["left"][i]
        crop = img.crop((max(left - 260, 0), max(top - 10, 0), left + 560, top + h + 10))
        for scale in (2, 3):
            big = ImageOps.expand(crop.resize((crop.width * scale, crop.height * scale)), border=30, fill=255)
            m = value_re.search(pytesseract.image_to_string(big, config="--psm 7"))
            if m:
                val = m.group(1).replace(",", "")
                found.append(float(val) if field == "mw" else int(val))
    return found

def reread_printed_fields(img, parsed):
    """Try re-read MW / extinction values against the already-extracted sequence.
    Returns (updates, description) or None."""
    seq, reported = parsed["sequence"], parsed["bases_reported"]
    if not seq:
        return None
    data = pytesseract.image_to_data(img, config="--psm 6", output_type=pytesseract.Output.DICT)
    mws = [parsed["mw_sheet"]] + [v for v in reread_field(img, "mw", data) if v != parsed["mw_sheet"]]
    exts = [parsed["ext_sheet"]] + [v for v in reread_field(img, "ext", data) if v != parsed["ext_sheet"]]
    for mw in mws:
        for ext in exts:
            if (mw, ext) == (parsed["mw_sheet"], parsed["ext_sheet"]):
                continue
            ev = evaluate(seq, reported, mw, ext)
            if ev["check"] == "OK":
                changed = [n for n, new, old in (("MW", mw, parsed["mw_sheet"]), ("ext coeff", ext, parsed["ext_sheet"])) if new != old]
                return {"mw_sheet": mw, "ext_sheet": ext, **ev}, "reread " + " + ".join(changed)
    return None

def attempt_fixes(img, parsed, strip_dashes=False):
    """Order the two fixes sensibly. If the length and extinction coefficient agree but only the MW doesn't,
    the sequence is almost certainly right and the printed MW was misread, so re-read it first."""
    mw_only = (parsed["sequence"] and parsed["bases_reported"] == parsed["bases_counted"]
               and parsed["ext_check"] == "match" and parsed["mw_check"] != "match")
    steps = ("reread", "refine") if mw_only else ("refine", "reread")
    for step in steps:
        if step == "refine":
            fixed = refine_sequence(img, parsed["bases_reported"], parsed["mw_sheet"], parsed["ext_sheet"],
                                    strip_dashes=strip_dashes)
            if fixed:
                return {"sequence": fixed[0],
                        **evaluate(fixed[0], parsed["bases_reported"], parsed["mw_sheet"], parsed["ext_sheet"])}, fixed[1]
        else:
            fixed = reread_printed_fields(img, parsed)
            if fixed:
                return fixed
    return None

def split_sheets(page_text):
    """A page may hold one or several sheets; split on the 'SPECIFICATION SHEET' header."""
    chunks = [c for c in SHEET_SPLIT_RE.split(page_text) if c.strip()]
    return chunks or [page_text]

# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--debug", action="store_true", help="save raw OCR text per page")
    args = ap.parse_args()

    _dash_rules()                      # fail fast if disregard_dashes_in_sequence is mis-typed
    pdf = Path(args.pdf)
    n_pages = pdfinfo_from_path(str(pdf))["Pages"]

    rows = []
    for i in range(1, n_pages + 1):          # one page at a time keeps memory low
        page = convert_from_path(str(pdf), dpi=args.dpi, first_page=i, last_page=i)[0]
        raw_text, upright = ocr_image(page)
        if args.debug:                                   # save exactly what Tesseract produced
            (pdf.parent / f"{pdf.stem}_page{i}.txt").write_text(raw_text)
        text, text_notes = apply_text_fixes(raw_text, pdf.name, i)   # round 1 + round 2, before any checks
        strip_dashes = dashes_enabled(pdf.name, i)
        j = 0
        for chunk in split_sheets(text):
            parsed = parse_sheet(chunk)
            if not parsed["name"] and not parsed["sequence"] and parsed["tm_c"] is None:
                continue  # header fragment / blank / continuation page
            parsed["fix"] = ""
            parsed["text_fixes"] = "; ".join(text_notes)
            if parsed["check"] not in OK_STATES:
                fixed = attempt_fixes(upright, parsed, strip_dashes=strip_dashes)
                if fixed:
                    parsed.update(fixed[0], fix=fixed[1])
            j += 1
            row = {"file": pdf.name, "page": i, "sheet_on_page": j, **parsed}
            rows.append(row)
            print(f"p{i}.{j}: {row['name']} | {row['bases_reported']} bases | Tm {row['tm_c']} | "
                  f"{row['check']} {('(' + row['fix'] + ')') if row['fix'] else ''}"
                  f"{(' [text fixes: ' + row['text_fixes'] + ']') if row['text_fixes'] else ''}", flush=True)

    out = pdf.with_suffix(".csv")
    fields = ["file", "page", "sheet_on_page", "name", "sequence", "bases_reported",
              "bases_counted", "tm_c", "folding_tm_c", "mw_sheet", "mw_calc", "mw_check", "ext_sheet", "ext_calc", "ext_check", "check", "fix", "text_fixes"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {out}")

if __name__ == "__main__":
    main()
