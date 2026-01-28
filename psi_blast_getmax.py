from util import *
import re
import numpy as np
import sys, getopt
from Bio import SeqIO
import logging

USAGE = '''
This script runs psi-blast with specified e-value and etracts sequences from the iteration with max number of hits

python ''' + sys.argv[0] + '''

-i || --input                  -input file with protein sequence
[-o || --psiblastout]          -output file with psi-blast result  
[-n || --psiblast_max_out]     -output file with the part of psi-blast result, containing all the hits from the iteration 
								with max number of hits
[-t || --refseqs_out]		   -output file with refesqs extracted from the psiblast iteration with max numbe rof hits
[-k || --protens_out]          -output file with protien sequences corresponding to the iteration with the max numbe rof hits
[-a || --aligned_out]          -output file with the protions of protein sequences corresponding to the aligned part in psiblast run 
[-r || --iteration]            -number of psiblast iteration to run
[-m || --max_target_seq]       -maximum number of sequences resulting from psiblast run to keep 
[-e || --evalue]               -evalue threashold to save hits
[-b || --blastdb_path]         -full path to blast protein refseq database
[-l || --psiblast]             -full path to psiblast program
[-c || --blastdbcmd]           -full path to blastdbcmd program
'''

LOGGER = logging.getLogger('prepare_data_construct_tree')
LOGGER.setLevel(logging.DEBUG)

fh = logging.FileHandler('hmm_pipeline.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
LOGGER.addHandler(fh)
LOGGER.addHandler(ch)

#global variables
PSIBLAST_INPUT = "protein.fa"
PSIBLAST_OUTPUT = "resulting_file.txt"
PSIBLAST_OUTPUT_MAX_HITS = "resulting_file_max_hits.txt"
REFSEQS_OUTPUT = "refseqs_out.txt"
PSI_ROUNDS = 3
MAX_TARGET_SEQS = 60000
EVAL_TRESHOLD = 0.049
#blast database full path
BLAST_DB_PATH = "/ssd2/refseqDB/refseq_protein"
#psiblast program full path
PSIBLAST = "psiblast"
BLASTDBCMD = "blastdbcmd"
PROTEIN_SEQUENCES = "protein_sequences.fa"
ALIGNED_REGIONS_OUTPUT = "aligned_regions_out.fasta"

def initialize(argv):
	global PSIBLAST_INPUT, PSIBLAST_OUTPUT, PSIBLAST_OUTPUT_MAX_HITS, REFSEQS_OUTPUT, PROTEIN_SEQUENCES, ALIGNED_REGIONS_OUTPUT
	global PSI_ROUNDS, MAX_TARGET_SEQS, EVAL_TRESHOLD
	global BLAST_DB_PATH, PSIBLAST, BLASTDBCMD
	LOGGER.info('Initializing parameters') 
	try:
		opts, args = getopt.getopt(argv[1:],"hi:o:n:t:k:s:r:m:e:b:l:c",["input=", "psiblastout=", "psiblast_max_out=", \
		"refseqs_out=", "protens_out=", "aligned_out=", "iteration=", "max_target_seq=", "evalue=", "blastdb_path=", \
		"psiblast=", "blastdbcmd="])
		if len(opts) == 0:
			raise getopt.GetoptError("Options are required\n")
	except getopt.GetoptError as e:
		LOGGER.error("===========ERROR==========\n " + str(e) + USAGE)
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print (USAGE)
			sys.exit()
		elif opt in ("-i", "--input"):
			PSIBLAST_INPUT = str(arg).strip()
		elif opt in ("-o", "--psiblastout"):
			PSIBLAST_OUTPUT = str(arg).strip()
		elif opt in ("-n", "--psiblast_max_out"):
			PSIBLAST_OUTPUT_MAX_HITS = str(arg).strip()
		elif opt in ("-t", "--refseqs_out"):
			REFSEQS_OUTPUT = str(arg).strip()
		elif opt in ("-k", "--protens_out"):
			PROTEIN_SEQUENCES = str(arg).strip()
		elif opt in ("-s", "--aligned_out"):
			ALIGNED_REGIONS_OUTPUT = str(arg).strip()   
		elif opt in ("-r", "--iteration"):
			PSI_ROUNDS = int(arg) 
		elif opt in ("-m", "--max_target_seq"):
			MAX_TARGET_SEQS = int(arg) 
		elif opt in ("-e", "--evalue"):
			EVAL_TRESHOLD = float(arg) 
		elif opt in ("-b", "--blastdb_path"):
			BLAST_DB_PATH = str(arg).strip()
		elif opt in ("-l", "--psiblast"):
			PSIBLAST = str(arg).strip()
		elif opt in ("-c", "--blastdbcmd"):
			BLASTDBCMD = str(arg).strip()

def identify_longest_iteration():
	#identify every iteration in the psi-blast output, find the one with most hits, and extract all protein names in 
	#one file and extract the subset of psi-blast run with max number of hits to a separate file
	LOGGER.info('Identifying the longest psiblast iteration')
	psiblast() 
	with open(PSIBLAST_OUTPUT, 'r') as psi_result:
		header_match = re.findall(r"# (.*) hits found", psi_result.read())
		header_match2 = map(int, header_match)
		max_num = max(header_match2)
		max_idx = np.argmax(header_match2)+1
	with open(REFSEQS_OUTPUT, 'w') as out_file, open(PSIBLAST_OUTPUT_MAX_HITS, 'w') as max_hits_out:
		with open(PSIBLAST_OUTPUT, 'r') as psi_result:
			max_num_hit_found = False
			query_size = 0
			for line in psi_result:
				if re.match("# {0} hits found".format(max_num), line):
					max_num_hit_found = True
					continue
				if max_num_hit_found:
					out_file.write(line.split()[1] + "\n")
					max_hits_out.write(line)
					query_size += 1
					if query_size >= max_num:
						break
			
def extract_aligned_regions():
	LOGGER.info('Extracting aligned regions') 
	blastdbcmd()
	refSeqToCoords = dict()
	with open(PSIBLAST_OUTPUT_MAX_HITS, 'r') as psi_max_hits, open(ALIGNED_REGIONS_OUTPUT, 'w') as aligned_regions:
		for record in psi_max_hits:
			record = record.split("\t")
			refSeqToCoords[record[1]] = (int(record[8]), int(record[9]))
		#We don't use SeqRecord to save the sequences because the standart
		#way is much faster 
		for record in SeqIO.parse(open(PROTEIN_SEQUENCES, "r"), "fasta"):
			refseqId = record.name
			start = refSeqToCoords[refseqId][0] - 1
			end = refSeqToCoords[refseqId][1] - 1
			coveredPart = str(record.seq[start:end+1])			
			aligned_regions.write(">" + prepareNames(record.description) + "\n")
			aligned_regions.write(coveredPart+ "\n")

#In production psi-blast input file name will be provided as parameter 
def psiblast():
	#run psi blast using parameters above
	psiblast_commandline = '{0} -db {1} -outfmt 7 -query {2} -out {3} -num_iterations {4} ' \
	'-max_target_seqs {5} -evalue {6}'.format(PSIBLAST, BLAST_DB_PATH, PSIBLAST_INPUT, PSIBLAST_OUTPUT,  \
	PSI_ROUNDS, MAX_TARGET_SEQS, EVAL_TRESHOLD)
	LOGGER.info('Launching psiblast : \n' + psiblast_commandline) 
	runSubProcess(psiblast_commandline)

def blastdbcmd():
	blastdbcmd_commandline = '{0} -entry_batch {1} -db {2} -outfmt "%f" -out {3}' \
	.format(BLASTDBCMD, REFSEQS_OUTPUT, BLAST_DB_PATH, PROTEIN_SEQUENCES)
	LOGGER.info('Launching blastdbcmd : \n' + blastdbcmd_commandline) 
	runSubProcess(blastdbcmd_commandline)
	
def main(argv):
	initialize(argv)
	identify_longest_iteration()
	extract_aligned_regions()
	
if __name__ == "__main__":
   main(sys.argv)
