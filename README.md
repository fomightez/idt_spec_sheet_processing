# idt_spec_sheet_processing
Python scripts and Jupyter Notebooks to extract data from IDT spec sheets for primers/oligonucleotides and combine into a useful form.

-------

### Quick Overview

Go from IDT spec sheets, either hardcopy or downloaded, to a useful spreadsheet/text data table.


Format for the generated final entries for the spreadsheet/text data table:
primer sequence <tab> name of the primer <tab> comments 

(That is essentially a remix of the standard FASTA format with the 'comments' being analogous to the optional description on the header line.)

The important part is that the software here is to avoid errors in the sequence and names of the primer by making it programmatic to collect and process the data from the IDT spec sheets to result in a primer spreadsheet/text data table.  The 'comments' field is the only part you'd want to generally be editing by hand.



-------


### Looking for a free & convenient place to use the scripts and notebooks without touching your local machine?

I have been using MyBinder-provided Jupyter sessions launched from [my 3Dscatter_plot_mod_playground-binder](https://github.com/fomightez/3Dscatter_plot_mod_playground-binder) repo to develop & run these scripts and notebooks after uploading them and the necessary data. Go to [there](https://github.com/fomightez/3Dscatter_plot_mod_playground-binder) and click the '`lauch binder`' to get started in a Jupyter session without logging in or signing up for anything. (MyBinder, by the Jupyter Community, provides temporary, anonymous remote sessions for  computational efforts with that aren't overly demanding.) 


---------


### HOW-to-USE GUIDE

I'd suggest running in MyBinder-provided Jupyter sessions launched from [my 3Dscatter_plot_mod_playground-binder](https://github.com/fomightez/3Dscatter_plot_mod_playground-binder). repo 

What do you need?
- A place to run Jupyter Notebooks using the Python kernel (not WASM-based)
- PDF files of IDT spec sheets (orientation doesn't matter)
- The appropriate Jupyter Notebook or notebooks to run

For the PDF files of the spec sheets, you can use either scanned from hardcopies or downloaded from IDT. Each spec sheet should be on a separate page in the PDF file. You can combine several pages into a single PDF as you scan hardcopies. You don't want the file size to get to giant as it will make them hard to email yourself. You can has as many PDF files as you want. The downloaded sheets from IDT come with additional pages that you don't want, and so you need to do some extra steps to extract and get each on a single page. You can use the steps in `idt_unzip_and_extract_first_pages_of_download_order_spec_sheet_zips.ipynb` BEFORE THE combine step to do this programmatically.




When you have Jupyter running and IDT spec sheets as single pages in whatever number of PDF files you end up with, the next step left is just running the appropriate provided Jupyter Notebooks and most things will be handled.  If you are just getting started, you'd use `Collecting_data_from_scanned_IDT_Spec_Sheets.ipynb`. Once you have an established main listing, you switch to using `Updating_primer_list_with_data_from_scanned_IDT_Spec_Sheets.ipynb` going forward to add more primers as you need to update.

The Optical Character Recognition (OCR) from the spec sheets isn't always perfect; however, the Jupyter notebook/scripts will signal those that need additional scrutiny. Look for `REVIEW: checks failed, verify against the scan` or `only length was verified` in the comments field to highlight those.

In addition to reviewing, you'll need to do any applicable editing & annotating by hand to the `.tsv` file produced and use it where you'd like, for example on Google Sheets to share with your lab members.


---------

### Details of the software components

The guidance above should be sufficient for most users. This section is meant for those more curious and/or troubleshooting.

These script and notebooks are meant for mining scanned spec sheets you may have if you order your primers/oligonucleotides from [Integrated DNA Technology (IDT)](https://www.idtdna.com/), and then combines the mined data into a spreadsheet/text table you can use going forward as your laboratorys's primer listing. Then let's you mine more scanned later and continue to combine to keep updating as more spec sheets are found or orders placed.

The scripts and notebooks:
- `idt_primer_spec_sheet_ocrOct29-2026.py` - script that mines details scans of IDT spec sheets (as PDFs) using optical character recognition (OCR). The script doesn't care about orientation of the individual scanned pages and it validates the scanned information with checks of the extinction coefficient and molecular weight and warns when not verified and things should be reviewed by user inspection. For scanning-generated text errors that keep coming up, there is a dictionary of exact, case-sensitive replacements applied to the whole page text to which you can add entries as you find recurring OCR errors. Tries to ignore stray pencil or pen marks in the area of sequences on the spec sheets. This script produces a `.csv` file with many details about the primers and how and from where they were mined. There will be one `.csv` produced per set of `.pdf` files of the spec sheets.
- `idt_primer_spec_sheets_combine.py` - any number of `.csv` files packed with the mined information from the scan step get combined to one set of mined data with redundant names of the same sequences and other matters handle. To faciliate maintaining & updating as you go forward, the script can take an established primer list and treat that as the starting point so that only information gets added to that as the additional mined primers/ oligonucleotides data are incoporated into the latest version of the main, established list of lab primers. To run the script, point it at the directory where you have stored your `.csv` files. The script produces a `.tsv` file that a de-duplicated listing of all the entries from all the `.csv` files, or any additional supplied listing, in one `.tsv` in  a 'boiled-down' form that is the spreadsheet/data table that should have everything needed from the mined data and can subsequently have the comments hand-edited to make it even better. 
- `Collecting_data_from_scanned_IDT_Spec_Sheets.ipynb` - makes use of `idt_primer_spec_sheet_ocrOct29-2026.py` and `idt_primer_spec_sheets_combine.py` to make a spreadsheet from a large collection of scans of IDT spec sheets for primers. (This Jupyter Notebook file assumes you are just starting out and making a list from several scans; if you aren't at the the initial stage and already have an established primer list, see `Updating_primer_list_with_data_from_scanned_IDT_Spec_Sheets.ipynb`.)
- `merge_primerDetails_n_comments.py` - comments section in two different files --> for each sequence entry shared in both files, the comments entry in the main file gets length and Tm added from the extracted data file , and the alternative names get appended in the 'comments' entry from the extracted data, too. This script is meant for when you have a list of primers you made and want to incorporate some of the information from the spec sheet scans.
- `idt_unzip_and_extract_first_pages_of_download_order_spec_sheet_zips.ipynb` - the downloaded spec sheets are longer than the single sheets and so this adds handling the downloaded spec sheets and how you can process them to make a spreadsheet that can be used as the start of establishing a spreadsheet of primers
- `Updating_primer_list_with_data_from_scanned_IDT_Spec_Sheets.ipynb` -  makes use of `idt_primer_spec_sheet_ocrOct29-2026.py` and `idt_primer_spec_sheets_combine.py` to add scans of more spec sheets to an established list of primers. Very much like `Collecting_data_from_scanned_IDT_Spec_Sheets.ipynb` above, except assumes you already have an established, main list of primers and want to add the data from some additional scanned spec sheets to it.

### Caveats

Excel opens `combined_primers.tsv` wrong when double-clicked (`°` garbled with default encoding, or no columns with utf-8-sig).  
Easiest solution: open the `.tsv` file in Sublime Text (or other text editor; on VSCode is acceptable from Microsoft productes), and then copy-paste contents into Excel.

For now, the OCR looking at spec sheets won't handle duplexes. 

