from Bio import Entrez, SeqIO
import pandas as pd
import os
import csv
from tqdm import tqdm


def fetch_using_id(accesion, db):
    """Fetch a single accesion id from database in text format

    Args:
        accesion (str): Accesion identifier for a single sequence.
        db (str): Database identifier, ej nucleotide

    Returns:
        SeqIO Record:  standard Sequence Input/Output interface for BioPython contains all data obtained from database
    """
    # Initialize parameters
    db = db
    paramEutils = {'usehistory': 'y'}

    # Generate query to Entrez search
    with Entrez.efetch(db, id=accesion, rettype='gb', retmode='text') as handle:
        # get Esearch result as dict objec
        record = SeqIO.read(handle, 'genbank')
    return record


def search_ncib(term, db='nucleotide', retmax=20):
    """Makes a search in db using the term query, and downloads all sequences found.

    Args:
        term (str): Search query produced in the website
        db (str, optional): Database used to make this search. Defaults to 'nucleotide'.
        retmax (int, optional): Max number of sequences obtained from database, max 100000. Defaults to 20.
    """

    with Entrez.esearch(db=db, retmax=retmax, term=term, idtype='acc') as handle:
        results = Entrez.read(handle)

    accesion_list = results['IdList']
    sequence_list = []

    confirm = input(f'You are about to download {len(accesion_list)} files.\n Do you want to continue? [Y/n] ')
    if confirm in ('y', 'Y'):
        for accesion in tqdm(accesion_list):
            sequence_list.append(fetch_using_id(accesion, db))
    return sequence_list


def write_sequences_to_disk(sequences, multiple_files=False, output_dir='output'):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if multiple_files:
        for seq in sequences:
            path = os.path.join(output_dir, f'{seq.id}.fasta')
            SeqIO.write(seq, path, 'fasta')
    else:
        path = os.path.join(output_dir, f'50_sequences.fasta')
        SeqIO.write(sequences, path, 'fasta')


if __name__ == "__main__":
    # Always tell who you are to NCIB
    email = 'druedaplata@protonmail.com'
    Entrez.email = email

    db = 'nucleotide'
    search_query = '(("Severe acute respiratory syndrome coronavirus 2"[Organism] OR sars cov 2[All Fields]) AND genome[All Fields] AND complete[All Fields]) AND ("2020/01/01"[PDAT] : "2020/04/31"[PDAT])'

    # get all ids to download
    
    sequences = search_ncib(search_query, db, retmax=50)
    print('downloaded')
    # write sequences to disk

    write_sequences_to_disk(sequences)


    # fetch all ids

    # save all ids to csv file