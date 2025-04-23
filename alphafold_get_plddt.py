### taken and modified from (https://github.com/BioComputingUP/AlphaFold-disorder)
#Piovesan D, Monzon AM, Tosatto SCE.
#Intrinsic protein disorder and conditional folding in AlphaFoldDB. Protein Sci. 2022 Nov;31(11):e4466.
#PMID: 36210722 PMCID: PMC9601767.

### Description: get pLDDT score from list of pdb files
### input: a folder containing pdb files; the folder should be in the same directory as the py file or mention full path to the folder in the input
### output: tsv file with pLDDT score for all the pdb files

### to run this type following command in the cmd
### 	'python3 alphafold_get_pLDDT.py -i pdbs/ -o out.tsv'
### note: 'pdbs/' is path to "pdbs" folder containing pdb files


from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.SeqUtils import seq1
from Bio.PDB import DSSP
import numpy as np
import warnings
import pandas as pd
import argparse
import logging.config
import sys
import csv
from pathlib import Path, PurePath
import tempfile
import gzip
import shutil
import os

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def process_pdb(pdb_file, pdb_name):
    # Decompress the structure if necessary
    real_file = pdb_file
    if is_gz_file(pdb_file):
        file_ext = Path(Path(pdb_file).stem).suffix
        fd, real_file = tempfile.mkstemp(prefix="alphafold-disorder_", suffix=file_ext)
        with open(real_file, "wb") as tmp:
            with gzip.open(pdb_file) as pdbf:
                shutil.copyfileobj(pdbf, tmp)
        os.close(fd)

    # Load the structure
    file_ext = Path(real_file).suffix
    if file_ext in ['.pdb']:
        structure = PDBParser(QUIET=True).get_structure('', real_file)
    else:
        # assume mmCIF
        structure = FastMMCIFParser(QUIET=True).get_structure('', real_file)


    # Remove decompressed if necessary
    if real_file != pdb_file:
        Path(real_file).unlink()

    # Parse b-factor (pLDDT) and DSSP
    df = []
    for i, residue in enumerate(structure.get_residues()):
        lddt = residue['CA'].get_bfactor() / 100.0
        
        df.append((pdb_name, i + 1, seq1(residue.get_resname()), lddt, 1 - lddt))
    df = pd.DataFrame(df, columns=['name', 'pos', 'aa', 'lddt', 'disorder'])
    return df

def parse_args():
    parent_parser = argparse.ArgumentParser(add_help=False)

    group = parent_parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--in_struct', type=str,
                       help='A single file, folder or file listing containing (gzipped) PDB or mmCIF files (relative paths)')
    #group.add_argument('-d', '--in_dssp', type=str, help='A TSV file with RSA and pLDDT columns (checkpoint file)')

    parent_parser.add_argument('-o', '--out', type=str, required=True,
                               help='Output file. Automatically generate multiple files using this name, ignore extention')

    parent_parser.add_argument('-f', '--format', type=str, choices=['tsv'], default='tsv', help='Output format')
    
    parent_parser.add_argument('-ll', type=str, choices=['notset', 'debug', 'info', 'warning', 'error', 'critical'],
                               default='info', help='Log level')

    main_parser = argparse.ArgumentParser(parents=[parent_parser])

    return main_parser.parse_args()


def process_file(f):
    result = pd.DataFrame([])
    if f.stat().st_size > 0:  # and 'P52799' in file.stem:  # 'P13693', 'P52799', 'P0AE72', 'Q13148'
        logging.debug('Processing PDB {}'.format(f))
        result = process_pdb(f, f.stem.split('.')[0])
    else:
        logging.debug('Empty file {}'.format(f))
    return result


if __name__ == '__main__':

    # parse command line arguments
    args = parse_args()
    fout_path = Path(args.out)

    # Set logger
    logging.basicConfig(format='%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s',
                        level=logging.getLevelName(args.ll.upper()), stream=sys.stdout)
    logging.getLogger('numexpr').setLevel(logging.WARNING)  # Remove numexpr warning

    # Disable pandas warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)

    if args.in_struct:
        # Generate DSSP output from PDB files
        data = pd.DataFrame()
        p = Path(args.in_struct)
        if p.is_file():
            # input is a single struct file or file with list
            if ''.join(PurePath(p).suffixes) in ['.pdb', '.pdb.gz', '.cif', '.cif.gz']:
                # process single file as input
                processed_data = process_file(p)
                if not processed_data.empty:
                    data = data.append(processed_data)
            else:
                # process list of files as input (paths in list are relative)
                with open(p, 'r') as list_file:
                    for file in list_file:
                        real_file = Path(p.parent, Path(file.strip()))
                        processed_data = process_file(real_file)
                        if not processed_data.empty:
                            data = data.append(processed_data)
        else:
            # input is a directory
            for file in p.iterdir():
                processed_data = process_file(file)
                if not processed_data.empty:
                    data = data.append(processed_data)    
    else:
        data = None

    pred = pd.DataFrame()
    for name, pdb_data in data.groupby('name'):
        pred = pred.append(pdb_data.copy())
	
    # Write to file
    if args.format == 'tsv':
        fout_name = '{}/{}_pred.tsv'.format(fout_path.parent, fout_path.stem)
        pred.to_csv(fout_name, sep='\t', quoting=csv.QUOTE_NONE, index=False, float_format='%.3f')
        logging.info('Prediction written in {}'.format(fout_path))
