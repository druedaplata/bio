# Script to download genomic sequences from Genbank

## Requirements

1. Use `python3.7` or higher versions
2. Install all packages listed in `requirements.txt`
3. This could be done with `pip install -r requirements.txt`

## How to use?

1. Open `config.ini` and change the parameters according to your needs.
    - **Database**: Name of the database in genbank
    - **Search Query**: Make a search in genbank and copy the information in search details.
    - **Email**: Use an email registered in Genbank
    - **Max_Items_Retreived**: If left unchanged, all sequences found will be downloaded. Change if you need a quick result for testning purposes.
    - **Publication Years**: List of years to limit the search range, separated by commas, ej. 1998,1999,2000
    - **Publication Months**: List of months to limit the search range, separated by commas, ej. 01,02,03,04,10,12

2. Run the script `data_download.py` using python.