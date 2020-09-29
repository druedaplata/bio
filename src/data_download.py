from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import configparser
import pandas as pd
import numpy as np
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
    paramEutils = {'usehistory': 'y'}

    # Generate query to Entrez search
    with Entrez.efetch(db, id=accesion, rettype='gb',
                       retmode='text') as handle:
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

    with Entrez.esearch(db=db, retmax=retmax, term=term,
                        idtype='acc') as handle:
        results = Entrez.read(handle)
        accesion_list = results['IdList']
    return accesion_list


def write_sequences_to_disk(items, database, output_dir):

    confirm = input(
        f'You are about to download {len(items)} files.\n Do you want to continue? [Y/n] '
    )
    if confirm in ('y', 'Y'):

        # Fetch all individual items from database
        sequences = []
        for accesion in tqdm(items):
            sequences.append(fetch_using_id(accesion, database))

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        if not os.path.exists(f'{output_dir}/seqs'):
            os.mkdir(f'{output_dir}/seqs')

        columns = [
            'id', 'organism', 'host', 'country', 'isolate', 'isolation_source',
            'molecule_type', 'dbxref', 'topology', 'taxonomy',
            'collection_date', 'reported_date', 'n_count'
        ]
        meta = pd.DataFrame(columns=columns)

        records = []
        for seq in sequences:

            # Check for N in each sequence
            max_len = len(seq.seq)
            n_count = np.sum(np.array(seq.seq) == 'N')
            print(f'%N = {n_count/max_len}')

            # Write fasta file with id and sequence
            record = SeqRecord(seq.seq, f'{seq.id}', '', '')
            SeqIO.write(record, f'{output_dir}/seqs/{seq.id}.fasta', 'fasta')
            records.append(record)

            # Append meta data in dataframe
            _id = seq.id
            organism = seq.features[0].qualifiers.get('organism', ['NaN'])[0]
            host = seq.features[0].qualifiers.get('host', ['NaN'])[0]
            country = seq.features[0].qualifiers.get('country', ['NaN'])[0]
            isolate = seq.features[0].qualifiers.get('isolate', ['NaN'])[0]
            isolation_source = seq.features[0].qualifiers.get(
                'isolation_source', ['NaN'])[0]
            molecule_type = seq.features[0].qualifiers.get(
                'mol_type', ['NaN'])[0]
            dbxref = seq.features[0].qualifiers.get('db_xref', ['NaN'])[0]
            collection_date = seq.features[0].qualifiers.get(
                'collection_date', ['NaN'])[0]
            topology = seq.annotations['topology']
            reported_date = seq.annotations['date']
            taxonomy = '|'.join(seq.annotations['taxonomy'])

            data = pd.Series(
                (_id, organism, host, country, isolate, isolation_source,
                 molecule_type, dbxref, topology, taxonomy, collection_date,
                 reported_date, n_count),
                index=meta.columns)

            meta = meta.append(data, ignore_index='True')

        # Save all sequences in a single file
        SeqIO.write(records, f'{output_dir}/sequences.fasta', 'fasta')

        # Save dataframe to file
        meta.to_csv(f'{output_dir}/meta.csv', index=False)
        print('Finished.')

    else:
        print('No files downloaded...Closing.')


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read('config.ini')

    # Always tell who you are to NCIB
    email = config.get('DEFAULT', 'Email')
    Entrez.email = email
    database = config.get('DEFAULT', 'Database')
    search_query = config.get('DEFAULT', 'Search_Query')
    retmax = config.get('DEFAULT', 'Max_Items_Retreived')
    output_dir = config.get('DEFAULT', 'Output_Dir')

    pub_start_date = config.get('DEFAULT', 'Pub_Start_Date')
    pub_end_date = config.get('DEFAULT', 'Pub_End_Date')

    search_query += f'AND ({pub_start_date}[PDAT] : {pub_end_date}[PDAT])'

    # Search all items in query
    items = search_ncib(search_query, database, retmax=retmax)
    # Write all items in disk, if user wants.
    write_sequences_to_disk(items, database, output_dir)
