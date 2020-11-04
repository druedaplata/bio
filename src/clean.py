"""Short script to read all sequences and remove all sequences with any N's
"""

import os
from tqdm import tqdm
from Bio import SeqIO


for month in tqdm(('01','02','03','04','05','06','07','08','09','10','11','12')):
    clean_seqs = []
    path = f'{month}/sequences.fasta'
    if os.path.exists(path):
        with open(path, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if 'n' not in record.seq:
                    if 'N' not in record.seq:
                        clean_seqs.append(record)
        SeqIO.write(clean_seqs, f'{month}/sequences_clean.fasta', 'fasta')
                    
                                   
