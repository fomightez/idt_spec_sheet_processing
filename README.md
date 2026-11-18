# idt_spec_sheet_processing
Python scripts and Jupyter Notebooks to extract data from IDT spec sheets for primers/oligonucleotides and combine into a useful form.

---------

These are meant for mining scanned spec sheets you may have if you order your primers/oligonucleotides from [Integrated DNA Technology (IDT)](https://www.idtdna.com/), and then combines the mined data into a spreadsheet/text table you can use going forward as your laboratorys's primer listing. Then let's you mine more scanned later and continue to combine to keep updating as more spec sheets are found or orders placed.

The scripts and notebooks:
- `idt_primer_spec_sheet_ocrOct29-2026.py` - script that mines details scans of IDT spec sheets (as PDFs) using optical character recognition (OCR). The script doesn't care about orientation of the scanned pages and it to validate the scanned information with checks of the extinction coefficient and molecular weight and warns when not verified and things should be reviewed by user inspection. For scanning-generated text errors that keep coming up, there is a dictionary of exact, case-sensitive replacements applied to the whole page text to which you can add entries as you find recurring OCR errors. Tries to ignore stray pencil or pen marks in the area of sequences on the spec sheets.
- `idt_primer_spec_sheets_combine.py` - several mined detailed `.csv` files get combined to one set of mined data with redundant names of the same sequences and other matters handle. To faciliate maintaining & updating as you go forward, it can take an established primer list and treat that as the starting point so that only information gets added to that as the additional mined primers/ oligonucleotides data are incoporated.
- `Collecting_data_from_scanned_IDT_Spec_Sheets.ipynb` - makes use of `idt_primer_spec_sheet_ocr.py` and `idt_primer_spec_sheets_combine.py` to make a spreadsheet from a large collection of scans of IDT spec sheets for primers.
- `merge_primerDetails_n_comments.py` - comments section in two different files --> for each sequence entry shared in both files, the comments entry in the main file gets length and Tm added from the extracted data file , and the alternative names get appended in the 'comments' entry from the extracted data, too. This script is meant for when you have a list of primers you made and want to incorporate some of the information from the spec sheet scans.
- `idt_unzip_and_extract_first_pages_of_download_order_spec_sheet_zips.ipynb` - the downloaded spec sheets are longer than the single sheets and so this adds handling the downloaded spec sheets and how you can process them to make a spreadsheet that can be used as the start of establishing a spreadsheet of primers

### Caveats

Excel opens `combined_primers.tsv` wrong when double-clicked (`°` garbled with default encoding, or no columns with utf-8-sig).  
Easiest solution: open the `.tsv` file in Sublime Text (or other text editor; on VSCode is acceptable from Microsoft productes), and then copy-paste contents into Excel.

For now, the OCR looking at spec sheets won't handle duplexes. 

-------


## Looking for a free & convenient place to use the scripts and notebooks without touching your local system?

I have been using MyBinder sessions launched from [my 3Dscatter_plot_mod_playground-binder](https://github.com/fomightez/3Dscatter_plot_mod_playground-binder) repo to develop & run these scripts and notebooks after uploading them and the necessary data. Go to [there](https://github.com/fomightez/3Dscatter_plot_mod_playground-binder) and click the '`lauch binder`' to get started in a Jupyter session without logging in or signing up for anything. (MyBinder, by the Jupyter Community, provides temporary, anonymous remote sessions for  computational efforts with that aren't overly demanding.) 