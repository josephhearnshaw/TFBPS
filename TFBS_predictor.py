#!/usr/bin/env python
#Essential imports to algortihm
import pandas as pd
import numpy as np
import math
import JASPAR_utils as jspr # cutting down on lines
# The below is only to have a 'test' dataset easily accessible
import os
import subprocess
import shutil
import wget
import errno
import gzip
import shutil
import time
import argparse # Easier command line interface


def pseudocounter_normalizer(col, seq_count):
    """
        Normalizes each nucelotide count in the column and adds psuedocounts where necessary
        Returns the respective counts
    """
    # Normalise counts
    count_A, count_T = col.count('A')/seq_count, col.count('T')/seq_count
    count_C, count_G = col.count('C')/seq_count, col.count('G')/seq_count
    # Add pseudocounts to avoid zero
    psuedo_check = lambda count: 0.8 if count == 0 else count
    return psuedo_check(count_A), psuedo_check(count_T), \
           psuedo_check(count_C), psuedo_check(count_G)


def PPM(seqs):
    """
    Count nucleotide occurences at each position within the column to create the PFM
    The PFM is then converted to a PPM
    Returns a matrix containing the PPM.
    """
    seqs_len = [len(l) for l in seqs] # List of each column length
    if max(seqs_len) != min(seqs_len):
        raise Exception("Sequence lengths don't match, can't proceed!")
    PPM_list, seq_len, seq_count = [], seqs_len[0], len(seqs_len) # Initate list & dict, & get number of sequences
    # Count the number of occurences of each nucleotide and append to the list
    for x in range(0, seq_len):
        col = "" #  Placeholder of the vertical aix in the matrix
        for i in range(0, len(seqs)):
            col += seqs[i][x] # Store the sequences in horizontal and vertical (i.e. A matrix)
        count_A, count_T, count_C, count_G = pseudocounter_normalizer(col, seq_count)
        PPM_dict=dict(A=count_A, T=count_T, C=count_C, G=count_G) # Add the counts to the dict for each nucleotide
        PPM_list.append(PPM_dict) # Add to the list to make a list of dicts
    return PPM_list


def PWM_maker(seqs):
    """
        Creates PWM from PPM
        Returns the PWM as a matrix
    """

    log2_likelihood = lambda f, bg: math.log2((f/bg)) # Calculates the base 2 log likliehood of each element
    PWM_matrix, PPM_matrix = [], PPM(seqs) # Initiate the PWM_matrix list#
    for dicts in PPM_matrix:
        col = [dicts[count] for count in 'ACGT'] # Create columns ACGT and add respective counts
        col = [log2_likelihood(element, (1/len(col))) for element in col] # Compute PWMs
        PWM_matrix.append(col) # Add the column to the PWM_matrix
    df_matrix = pd.DataFrame(PWM_matrix)
    df_matrix_transpose = df_matrix.transpose()
    df_matrix_transpose.rename(index={0:'A',1:'C', 2:'G', 3:'T'}, inplace=True)
    return df_matrix_transpose


def nuc_to_np_num(str_seq):
    """
    Convert nucleotides to integers with coding ACGTN = 01230
    Note that N is defaulted to equal A as it's often non problematic,
    at least with well referenced/assembled  genomes!
    Returns a numpy array with unit8 unsigned int of same length
    """
    np_sequence = np.zeros((len(str_seq), ), np.uint8)
    ref = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    for i, nuc in enumerate(str_seq):
        np_sequence[i] = ref[nuc]

    return np_sequence


def np_num_to_nuc(numpy_sequence):
    """
    Converts the numpy array back to strings
    """
    sequence_string = ['A' for i in range(len(numpy_sequence))]
    str_dict = {0:'A', 1:'C', 2:'G', 3:'T'}
    for i, no in enumerate(numpy_sequence):
        sequence_string[i] = str_dict[no]

    return ''.join(sequence_string)


def scan_PWM(PWM_matrix, sequences, threshold, outdir):
    """
    Scans PWM through the genome sequence(s) given
    """
    n_cols = PWM_matrix.shape[1]
    # Converted pd df PWM to np array & column indicies for the PWM from 0 to (n_cols-1)
    numpy_matrix, cols = PWM_matrix.to_numpy(), np.arange(n_cols)
    # Initate the DF
    colnames=['Score', 'Sequence', 'Start Pos', 'End Pos']
    matches=pd.DataFrame({colname:[] for colname in colnames}, columns=colnames)
    print("Scanning through your genome, now! This may take some time, so hold on.\n")
    # Iterate through the whole genome
    for i in range(len(sequences) - n_cols + 1):
        scan_window = sequences[i:(i+n_cols)] # Obtain window to scan in genome
        score = np.sum(numpy_matrix[scan_window, cols]) # Index correct score from 0 to (n_cols-1)
        if score > threshold:
            # Add the relevant information to appropriate columns
            matches.loc[len(matches)] = [score, np_num_to_nuc(scan_window),
                                   i + 1, i + n_cols]
    print(f"Finished, writing your data out!\n")
    matches.to_csv(f"{outdir}", sep="\t", index=None)



def parse_fasta(filename):
    """
        Iterate through the file, make sure the first line has a sequence,
        denoted by '>' & return this, otherwise it's not a sequence
    """
    with open(filename, 'r') as fileHandler:
        lines = fileHandler.read().splitlines()
    first = lines.pop(0)
    if first.startswith('>'):
        return ''.join(lines)
    else:
        return None

def load_sequence(seq):
    """
        User can either load sequence from FASTA file or directly as a sequence on commandline
    """
    # If it's a sequence beign used by the user on the command line then set it as so
    if set(seq.lower()) == set(['a', 'c', 'g', 't']):
        sequence = nuc_to_np_num(seq)
        return sequence
    # If the sequence isn't a sequence then it's the fasta header
    sequence = parse_fasta(seq)
    # You need to give a sequence, silly
    if sequence is None:
        Exception("Invalid fasta format")
        return
    sequence = nuc_to_np_num(sequence)
    return sequence


def download(URL, filepath):
    """
    Perform wget request to download whatever data
    """
    try:
        wget.download(URL, filepath)
    except socket.error as e:
        if e.errno == errno.ECONNREFUSED:
            raise  # Not error we are looking for
        pass  # Handle error here.

def gunZip(filepathIN, filepathOUT):
    """ Unzips .gz file extensions """
    filepathOUT = filepathIN.replace('.gz', '')
    with gzip.open(filepathIN, 'rb') as f_in:
        with open(filepathOUT, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def __test__():
    """
    Performs the CGI challenge! Apologies for the excessive use of UNIX commands, was running out of time.
    """

    homo_chr1_url="ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz"
    desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') # Get the Desktop dir
    try:
        shutil.rmtree(f"{desktop}/CGI_ai_challenge/")
    except OSError:
        pass
    os.mkdir(f"{desktop}/CGI_ai_challenge/")
    CGI_dir = f"{desktop}/CGI_ai_challenge"
    fasta_dir, out_dir = f"{CGI_dir}/Homo_sapiens.GRCh38.dna.chromosome.1.fa", f"{CGI_dir}/homo_spaiens_chr1_100000_bases.fa"
    download(homo_chr1_url, CGI_dir)
    gunZip(f"{CGI_dir}/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz",fasta_dir)
    os.remove(f"{CGI_dir}/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz")
    # Get first 100k bases
    command_seq_header = f"sed -n '1,1p' {fasta_dir} > {CGI_dir}/fasta_header.txt"
    subprocess.Popen([command_seq_header], shell=True)
    time.sleep(2) # Issue of doing something outside GIL is it won't wait until it's done, 2 seconds should be adequate.
    #Get the length of the fasta header
    with open(f"{CGI_dir}/fasta_header.txt", 'r') as file:
        text = file.read().strip().split()
        len_chars = sum(len(word) for word in text)

    no_to_cut=100000 + len_chars
    command = f'cat {fasta_dir} | head -c{no_to_cut} >> {out_dir}'
    subprocess.Popen([command], shell=True)
    command = f"sed -e '/^>/! s/[Nn]//g' {out_dir} >> {CGI_dir}/tmp.fa"# Remove N's
    subprocess.Popen([command], shell=True)

    myc_jasper_url="http://jaspar.genereg.net/api/v1/matrix/MA0147.3.pfm" # JASPAR MYC PFM
    download(myc_jasper_url, f"{desktop}/CGI_ai_challenge")
    print(f"\nData will be stored in the following directory: {desktop}/CGI_ai_challenge/ \nLoading sequence...\n")
    sequence = load_sequence(f"{CGI_dir}/tmp.fa") # Seq from chr1 - homo sapiens
    jsp_PWM = jspr.load_matrix(f"{CGI_dir}/MA0147.3.pfm") # PWM for MYC
    scan_PWM(jsp_PWM, sequence, 0.4, f"{CGI_dir}/results.txt") # Threshold varies from 0.4-0.9, went with lowest
    files_rm = [fasta_dir, f"{CGI_dir}/fasta_header.txt", f"{CGI_dir}/tmp.fa",
                out_dir, f"{CGI_dir}/tmp_file.jaspar", f"{desktop}/CGI_ai_challenge/MA0147.3.pfm"]
    [os.remove(x) for x in files_rm] # Clean up
    max_score = pd.read_csv(f"{CGI_dir}/results.txt", sep="\t", lineterminator='\n')
    high_score = max_score.loc[max_score['Score'].idxmax()]
    print(f"Highest scoring sequence had the following attributes: \n{high_score}")


def dna_seq_iterator(seq_dir):
    """
    Return sequences back in a list for further use
    """
    seqs = []
    with open(seq_dir) as f:
        lines = f.readlines()
        for line in lines:
            seq = line.strip().split("\t")
            seqs.append(seq)
    if len(seqs) > 1:
        sequences = sum(seqs, [])
        return sequences

    return seqs

def convert_to_df(matrix):
    df_matrix = pd.DataFrame(matrix)
    df_matrix_transpose = df_matrix.transpose()
    df_matrix_transpose.rename(index={0:'A',1:'C', 2:'G', 3:'T'}, inplace=True)
    return df_matrix_transpose

######## ####
### Start ###
#############


parser = argparse.ArgumentParser(
    description="This script will predict Transcription factor binding sites\n")
# Args init
parser.add_argument(
    '-dna_seqs', '--dnaseqsdir', help='Directory to sequence you wish to create a PWM for. Data must be tab delimited', required=False)
parser.add_argument(
    '-fasta', '--fastadir', help='Fasta file directory for the sequences to compare', required=False)
parser.add_argument(
    '-pfm', '--pfmdir', help='PFM directory, must be PFM format', required=False)
parser.add_argument(
    '-out', '--outdir', help='Directory for where you want your results to be stored', required=False)
parser.add_argument(
    '-testonly', '--test', help='Will look for matches in Chr1 with MYC from Jasper, no additional parameters needed. Set as 1 to perform.', required=True)
parser.add_argument(
    '-threshold', '--threshold', help='Threshold needed for scanning via PWM.', required=False)

# Initiate the arguments
args = vars(parser.parse_args())
dna_seq, fasta_dir, pfm_dir, out_dir, testBool = args['dnaseqsdir'], args['fastadir'],\
                                                 args['pfmdir'], args['outdir'], args['test']
threshold = args['threshold']

if testBool == "1":
    print("Performing test")
    __test__()
    import sys
    sys.exit()
else:
    pass

if dna_seq is not None and fasta_dir is not None and out is not None:
    dna_seq = dna_seq_iterator(dna_seq)
    matrix = PWM_maker(dna_seq)
    df_matrix = convert_to_df(matrix)
    fn = fasta_dir.split("/")[-1:]
    tmp_fn = fasta_dir.replace(fn[0], "tmp_file.fa")
    command = f"sed -e '/^>/! s/[Nn]//g' {fasta_dir} >> {tmp_fn}"# Remove N's
    sequence = load_sequence(tmp_fn)
    if threshold is None:
        threshold = 0.7
    scan_PWM(df_matrix, threshold, out)
if dna_seq is None and fasta_dir is not None and out is not None and pfm_dir is not None:
    fn = fasta_dir.split("/")[-1:]
    tmp_fn = fasta_dir.replace(fn[0], "tmp_file.fa")
    command = f"sed -e '/^>/! s/[Nn]//g' {fasta_dir} >> {tmp_fn}"# Remove N's
    subprocess.Popen([command], shell=True)
    sequence = load_sequence(tmp_fn)
    jsp_PWM = jspr.load_matrix(pfm_dir)
    if threshold is None:
        threshold = 0.7
    scan_PWM(jsp_PWM, sequence, threshold, out)
else:
    print("Please execute Python3 TFBS_predictor.py -h for the help menu to see what flags are avaiable.\n")
