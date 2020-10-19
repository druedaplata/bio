# Steps to download and align genomic sequences from Genbank in Guane

## Requirements

1. Use `python3.7` or higher versions

2. If running in Guane, load the module and activate the environment.
    - `module load Bioinformatics/Bioconda/python3`
    - `source activate biopython_env`

2. If not, install all packages listed in `requirements.txt`
    - This could be done with `pip install -r requirements.txt`

## How to use?

1. Open `config.ini` and change the parameters according to your needs.
    - **Database**: Name of the database in genbank
    - **Search Query**: Make a search in genbank and copy the information in search details.
    - **Email**: Use an email registered in Genbank
    - **Max_Items_Retreived**: If left unchanged, all sequences found will be downloaded. Change if you need a quick result for testning purposes.
    - **Publication Years**: List of years to limit the search range, separated by commas, ej. 1998,1999,2000
    - **Publication Months**: List of months to limit the search range, separated by commas, ej. 01,02,03,04,10,12

2. Run the script `data_download.py` using python and wait until it finishes.

3. Use the bash script `make.sh` if you need a joint metadata and fasta file for all sequences in a single year.

4. For alignment we use [mafft](https://mafft.cbrc.jp/alignment/software/), read the documentation to understand the parameters and modify them to your needs.

5. Create a file and copy the following script to align the sequences using mafft on Guane, read the [SLURM](https://slurm.schedmd.com/sbatch.html) documentation if you need more advanced options.

```bash
#!/bin/bash -l

#SBATCH --time 14-12:00:00
#SBATCH --job-name ncbi_align
#SBATCH --partition normal
#SBATCH --error ncbi_align.err
#SBATCH --output ncbi_align.out
#SBATCH --exclusive

# Load the module
module load Bioinformatics/Bioconda/python3
source activate mafft_env

# Align all sequences in the first month, change the input and output to align other months
mafft --auto --nomemsave --thread -1 --op 9 --ep 0.2 --lop 9 --lep 0.2 --LEXP 0.45 01/sequences_clean.fasta > 01/alignment_clean.fasta

```

6. Run this file, assuming it is called *align.job* using `sbatch align.job`
