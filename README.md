# idt_spec_sheet_processing
Python scripts and Jupyter Notebooks to extract data from IDT spec sheets for primers/oligonucleotides and combine into a useful form.

---------

The scripts and notebooks:
- `idt_primer_spec_sheet_ocr.py` - script that mines details scans of IDT spec sheets (as PDFs) using optical character recognition (OCR) 
- `idt_primer_spec_sheets_combine.py` - several mined detiled SCVs get combined to one set of mined data
- `Collecting_data_from_scanned_IDT_Spec_Sheets.ipynb` - makes use of `idt_primer_spec_sheet_ocr.py` and `idt_primer_spec_sheets_combine.py` to make a spreadsheet from a large collection of scans of IDT spec sheets for primers.
- `merge_primerDetails_n_comments.py` - comments section in two different files --> for each sequence entry shared in both files, the comments entry in the main file gets length and Tm added from the extracted data file , and the alternative names get appended in the 'comments' entry from the extracted data, too.
- `idt_unzip_and_extract_first_pages_of_download_order_spec_sheet_zips.ipynb` - the downloaded spec sheets are longer than the single sheets and so this adds handling the downloaded spec sheets and how you can process them to make a spreadsheet that can be used as the start of establishing a spreadsheet of primers

