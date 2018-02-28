#!/usr/bin/python3
# coding: utf-8

## This is an additinal script for th MARVEL pipeline of analysis and retrieaval of Viral Extended sequences
# Developed by Deyvid Amgarten
# Creative commons

import sys
import os
import subprocess
import datetime

# Function declarations

# Usage
def usage():
    print('generate_bins_from_reads.py version 0.1\nAdditional script from the MARVEL pipeline.\n\nUSAGE: generate_bins_from_reads.py -1 reads-R1.fastq -2 reads-R2.fastq -t num_treads \n\t-12 Files containing R1,R2 Illumina reads in fastq format\n\t-t Number of threads to be used\n\t-h This message')


# Verify arguments
def verify_arg(arg_list):
    if arg_list[1] == '-1':
        return True
    elif arg_list[1] == '-h':
        return False
    else:
        print('There is something wrong, use -h for information and usage')
        return False

# Run
def run_shell(command):
    try:
        subprocess.call(command,shell=True)
    except:
        print('Could not run: %s'%command)

###
### Main code
###

# Verify arguments
args_list = sys.argv
if verify_arg(args_list):
    pass
else:
    usage()
    quit()

# Variables
file_reads1 = args_list[2]
file_reads2 = args_list[4]
threads = args_list[6]


# Assembly with Spades
print(str(datetime.datetime.now()))
print('Assemblying reads into contigs with Spades. This may take awhile, please be patient...\n')
run_spades = 'spades -1 '+file_reads1+' -2 '+file_reads2+' -t '+threads+' -o spades_assembly'
try:
    subprocess.call(run_spades, shell=True)
except:
    print('Could not run Spades with reads: %s and %s'%file_reads1,file_reads2)
print('Assembly has finished!\n')

# Mapping reads to contigs
print('Mapping reads to contigs with bowtie2\n')
run_bowtie_build = 'bowtie2-build -f spades_assembly/scaffolds.fasta contigs_index'
run_shell(run_bowtie_build)

run_bowtie = 'bowtie2 -x contigs_index -1 '+file_reads1+' -2 '+file_reads2+'-S mapping.sam --no-unal'
run_shell(run_bowtie)

run_samtools = 'samtools view -Sb mapping.sam > mapping.bam'
run_shell(run_samtools)

print('Finished mapping!\n')

# Binning
print('Generating bins with Metabat2.\n')
try:
    os.stat('bins_folder/')
except:
    os.mkdir('bins_folder/')

run_metabat = 'runMetabat -m 1500 -s 10000 -o bins_folder -t '+threads+' spades_assembly/scaffolds.fasta mapping.bam'
run_shell(run_metabat)
print('Finished binning!')

print(str(datetime.datetime.now()))
print('All done. Bins are in the folder: bins_folder/')
print('Thank you for using the MAVEL pipeline')


