#!/usr/bin/env python
# coding: utf-8

# Edit April, 9th 2019
## This is the MARVEL Pipeline for analysis and retrieaval of Viral Long sequences
## This is a python second-half version, which receives bins as input
# Developed by Deyvid Amgarten
# Creative commons

# Import required Python modules
import numpy as np
from Bio import SeqIO
import re
import sys
import os
import shutil
import subprocess
from collections import Counter
import pickle
import datetime
import sklearn
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
import argparse

# Function declarations

# Auxiliary function to kmer_frequencies
def chunks(l, n):
    for i in range(0, len(l) - (n - 1)):
        yield l[i:i + n]


# Calculate kmer frequencies
def kmer_frequency(dna, kmer):
    freq = {}
    freq.update(Counter(chunks(dna, 3)))  # update with counter of dinucleotide
    if kmer in freq:
        for each in freq:
            freq[each] = freq[each] / len(dna)
        return (freq[kmer])
    else:
        return 0


# Run prokka
def run_prokka(binn, input_folder, threads):
    # Check the fasta format
    prefix = get_prefix(binn)
    # Filehandle where the output of prokka will be saved
    # output_prokka = open(str(prefix)+'prokka.output', mode='w')
    # Full command line for prokka
    command_line = ('prokka --kingdom Viruses --centre X --compliant --gcode 11 --cpus ' + threads + ' --force --quiet --prefix prokka_results_' + str(prefix) + ' --fast --norrna --notrna --outdir ' + input_folder + 'results/prokka/' + str(prefix) + ' --cdsrnaolap --noanno ' + input_folder + str(binn)).split()
    return_code = subprocess.call(command_line, stderr=subprocess.PIPE)
    # Check with prokka run smothly
    if return_code == 1:
        print("Prokka may not be correctly installed. Please check that.")
        quit()

# Get prefix from bins
def get_prefix(binn):
    if re.search('.fasta', binn):
        prefix = re.sub('.fasta', '', binn)
    else:
        prefix = re.sub('.fa', '', binn)
    return(prefix)

# Extract features from genbank record and prepare vector of features
def extract_features(record):
    # iterate each feature
    count = 0
    sum_cds_length = 0
    strand_shift = 0
    non_coding_spacing = []
    for feature in record.features:
        # This is a modification for erroneus translations
        if feature.type == "CDS":
            if 'translation' in feature.qualifiers:
                if re.search('\w\*', str(feature.qualifiers['translation'])) is None:
                    count += 1
                    start = feature.location.start
                    end = feature.location.end
                    sum_cds_length += (end - start)
                    if count == 1:
                        strand_prev = feature.location.strand
                        end_prev = end
                    else:
                        non_coding_spacing.append(start - end_prev)
                        if strand_prev != feature.location.strand:
                            strand_shift += 1
                    end_prev = end
                    strand_prev = feature.location.strand
            else:
                print('WARNING: Prokka predicted a CDS, but there is no translation. Record ID: ',record.id)

    if len(non_coding_spacing) > 1:
        density = count / (len(record.seq) / 1000)
        #mean_gene_size = sum_cds_length / count
        sum_spacing = 0
        for i in non_coding_spacing:
            sum_spacing += i
        #mean_spacing_size = sum_spacing / (len(non_coding_spacing) - 1)
        #ATG_freq = kmer_frequency(str(record.seq), 'ATG')
        # Return mean size of non conding regions and the density of genes
        # print(record.id, [mean_spacing_size, density, mean_gene_size, strand_shift, '/', count, flag])
        return ([density, strand_shift / count])
    else:
        #warnings_handle.write("WARNING:" + record.id + " has 1 or zero appropriate CDS features.")
        print('WARNING: '+record.id+' has 1 or zero appropriate CDS features (those are important for prediction).')


###
### Main code
###

# Set arguments
# Modification to use argparse
parser = argparse.ArgumentParser(description='Predic phage draft genomes in metagenomic bins.')
parser.add_argument('-i', '--input',action="store", required=True, dest="input_folder", help='Path to a folder containing metagenomic bins in .fa or .fasta format (required!)')
parser.add_argument('-t', action="store", dest="threads", default='1', help='Number of CPU threads to be used by Prokka and hmmscan (default=1)')
parser.add_argument('-m', '--models', default='models/all_vogs_hmm_profiles_feb2018.hmm', help='HMM models file')
parser.add_argument('-o', '--output', default='results/', help='Output results directory')
args = parser.parse_args()

# Greeting message
print('\n**Welcome to the MARVEL Tool!\n')
print('** Please cite: Amgarten DE, Braga LP, Da Silva AM, Setubal JC. MARVEL, a Tool for Prediction of Bacteriophage Sequences in Metagenomic Bins. Frontiers in Genetics. 2018;9:304.')

# Verify databases
if not os.path.isfile(args.models):
    parser.error(f'--models {args.models} does not exist. Please specify correct directory, or run python download_and_set_models.py if you have not already')
    quit()

# Create Filehandle for warnings
#warnings_handle = open('marvel-warnings.txt', 'w')

# Important variables
input_folder = args.input_folder
threads = args.threads

# Fix folders path if missing '/'
if not re.search('/$', input_folder):
    input_folder = input_folder+'/'
if not re.search('/$', args.output):
    args.output = args.output+'/'    

# Take the input folder and list all multifasta (bins) contained inside it
list_bins_temp = os.listdir(input_folder)
list_bins = []
count_bins = 0
# Empty folder
if list_bins_temp == []:
    print('**Input folder is empty. Exiting...\n')
    quit()
else:
    for each_bin in list_bins_temp:
        if re.search('.fasta$', each_bin, re.IGNORECASE):
            list_bins.append(each_bin)
            count_bins += 1
        elif re.search('.fa$', each_bin, re.IGNORECASE):
            list_bins.append(each_bin)
            count_bins += 1

if count_bins == 0:
    print('**There is no valid bin inside the input folder (%s).\nBins should be in \'.fasta\' or \'.fa\' format.\nExiting...'%input_folder)
    quit()

print('**Arguments are OK. Checked the input folder and found %d bins.\n' % count_bins)
print('**'+str(datetime.datetime.now()))

# Create results folder
try:
    os.stat(args.output)
except:
    os.mkdir(args.output)

#####
# PROKKA
#####
# Running prokka for all the bins multfasta files in input folder
# Perform a check in each bin, then call the execute_prokka function individually
# It may take awhile
count_prokka = 0
print('**Prokka has started, this may take awhile. Be patient.\n')
for binn in list_bins:
    # Verify bin size
    len_bin = 0
    for record in SeqIO.parse(input_folder + binn, 'fasta'):
        len_bin += len(record.seq)
    #FIX: If a bin is too short, skip it
    if len_bin < 2000:
        print('**MARVEL has found a bin, which is too short to code proteins (<2000pb). As CDSs are an import feature for MARVEL, we will be skipping this bin: '+binn)
        continue
    run_prokka(binn, input_folder, threads)
    count_prokka += 1
    if count_prokka % 10 == 0:
        print('**Done with %d bins...' % count_prokka)
print('**Prokka tasks have finished!\n')

### 
# Checkpoint
# Remove prokka output that has empty fasta files
###
prokka_results = os.path.join(args.output,'prokka')
bin_dirs = os.listdir(prokka_results)
skipped_bins = []
for bin_dirname in bin_dirs:
    bin_path = os.path.join(prokka_results, bin_dirname)
    for f in os.listdir(bin_path):
        fpath = os.path.join(bin_path, f)
        if fpath.endswith('.faa') and os.stat(fpath).st_size == 0:
            skipped_bins.append(bin_dirname)
            print('**Skipping {} - No valid protein fasta was found'.format(bin_dirname))
            print('Removing {}'.format(bin_path))
            shutil.rmtree(bin_path)

if len(skipped_bins) == len(list_bins):
    print('**Error: prokka produced 0-length fasta files for all bins')
    print('\tPlease check your input!')
    sys.exit(1)

####
# HMM SEARCHES
####
print('**'+str(datetime.datetime.now()))
print('**Starting HMM scan, this may take awhile. Be patient.\n')
#print(str(datetime.datetime.now()))
# Create a new results folder for hmmscan output
try:
    os.stat(args.output + 'hmmscan/')
except:
    os.mkdir(args.output + 'hmmscan/')

# Call HMMscan to all bins
prop_hmms_hits = {}
count_hmm = 0
for binn in list_bins:
    # Prefix for naming results
    prefix = get_prefix(binn)
    #FIX: If a bin is too short, skip it
    len_bin = 0
    for record in SeqIO.parse(input_folder + binn, 'fasta'):
        len_bin += len(record.seq)
    if len_bin < 2000 or (prefix in skipped_bins):
        continue
    command_line_hmmscan = 'hmmscan -o ' + args.output + 'hmmscan/' + prefix + '_hmmscan.out --cpu ' + threads + ' --tblout ' + args.output + 'hmmscan/' + prefix + '_hmmscan.tbl --noali models/all_vogs_hmm_profiles_feb2018.hmm ' + args.output + 'prokka/' + prefix + '/prokka_results_' + prefix + '.faa'
    # In case hmmscan returns an error
    try:
        subprocess.call(command_line_hmmscan, shell=True)
        #True
    except:
        print('**Error calling HMMscan:', command_line_hmmscan)
        quit()
    count_hmm += 1
    # Iteration control
    if count_hmm % 10 == 0:
        print('**Done with %d bins HMM searches...' % count_hmm)
    # Parse hmmscan output files and find out which bins have less than 10% of their proteins
    # without any significant hits (10e-10)
    num_proteins_bin = 0
    with open(args.output + 'prokka/' + prefix + '/prokka_results_' + prefix + '.faa', 'r') as faa:
        for line in faa:
            if re.search('^>', line):
                num_proteins_bin += 1
    dic_matches = {}
    with open(args.output + 'hmmscan/' + prefix + '_hmmscan.tbl', 'r') as hmmscan_out:
        for line in hmmscan_out:
            match = re.search('^VOG\d+\s+-\s+(\S+)\s+-\s+(\S+)\s+.+$', line)
            if match:
                if match.group(1) not in dic_matches:
                    dic_matches[match.group(1)] = float(match.group(2))
    # Take the proportion of proteins matching the pVOGs and store in a dictionary
    i_sig = 0
    for key in dic_matches:
        if dic_matches[key] <= 1e-10:
            i_sig += 1
    prop_hmms_hits[prefix] = i_sig / num_proteins_bin

###
# Machine learning prediction of viral bins
###

# Go for each bin and extract relevant features from respective genbank files
# Final feature vector is "array_bins"
print('**'+str(datetime.datetime.now()))
print('\n**Extracting features from bins...\n')
data_bins = []
data_bins_index = []
# Temporary fix: Some GBL files are broken
# Verify this

# Iteration for bins
count_pred = 0
for bins in list_bins:
    prefix = get_prefix(bins)
    #FIX: If a bin is too short, skip it.
    len_bin = 0
    for record in SeqIO.parse(input_folder + bins, 'fasta'):
        len_bin += len(record.seq)
    if len_bin < 2000 or (prefix in skipped_bins):
        continue
    count_pred += 1
    sub_data_bins = []
    #FIX: Check whether prokka generated a .gbk or .gbf file
    if os.path.isfile(args.output + 'prokka/' + prefix + '/prokka_results_' + prefix + '.gbk'):
        file_name = args.output + 'prokka/' + prefix + '/prokka_results_' + prefix + '.gbk'
    else:
        file_name = args.output + 'prokka/' + prefix + '/prokka_results_' + prefix + '.gbf'
    # GBK File with gene prediction
    #file_name = args.output + 'prokka/' + prefix + '/prokka_results_' + prefix + '.gbk'
    for record in SeqIO.parse(file_name, "genbank"):
        # Extract_features still take the class as an argument
        # Sending only the record now
        temp_list = extract_features(record)
        if temp_list is None:
            pass
        else:
            sub_data_bins.append(temp_list)
    ###
    if not sub_data_bins:
        continue
    sub_data_bins_array = np.array(sub_data_bins)
    mean_for_bin = list(np.mean(sub_data_bins_array, axis=0))
    ####
    # Append here, the fifth feature: proportions of hmm hits
    mean_for_bin.append(prop_hmms_hits[prefix])
    data_bins.append(mean_for_bin)
    data_bins_index.append(prefix)
array_bins = np.array(data_bins, dtype=float)
print('**Extracted features from', len(array_bins), 'bins\n')

# Load RFC model from file
print('**'+str(datetime.datetime.now()))
print('**Doing the machine learning prediction...\n')
pkl_filename = "models/pickle_model_rfc_trained_bins8k_refseq_all_3features_den_stran_prophitshmm.pkl"
with open(pkl_filename, 'rb') as file:
    pickle_model = pickle.load(file)

# Predict wether bins are from phage or other organisms
y_test_predicted_pickle = pickle_model.predict(array_bins[:, ])

# Assess the probability of the classification
y_test_prob = pickle_model.predict_proba(array_bins[:, ])

# Retrieve bins predicted as phages
# 0.5 will be the threshold of probability for the class 1
bins_predicted_as_phages = []
i = 0
is_there_phages = 0
for pred in y_test_prob[:, 1]:
    if pred >= 0.7:
        if is_there_phages == 0:
            print("**Found phages in this sample!!!")
            print("**Bins predicted as phages and probabilities according to Random Forest algorithm:")
            is_there_phages = 1
        bins_predicted_as_phages.append(data_bins_index[i])
        print("***", data_bins_index[i], "-> ",round(pred*100,2),"%")
        # print(pred, data_bins_index[i], array_bins[i])
    i += 1

### Alternative way: Not using probability
# bins_predicted_as_phages = []
# i = 0
# for pred in y_test_predicted_pickle:
#    if pred == 1:
#        bins_predicted_as_phages.append(data_bins_index[i])
#    i += 1


print('**Finished Machine learning predictions!\n')
#print(str(datetime.datetime.now()))
# Just make sure to end the program closing the warnings filehandle
#warnings_handle.close()

if is_there_phages == 1:
    try:
        os.stat(args.output + 'phage_bins')
    except:
        os.mkdir(args.output + 'phage_bins')
    if re.search('.fasta',list_bins[0], re.IGNORECASE):
        file_format = '.fasta'
    else:
        file_format = '.fa'
    for bin_phage in bins_predicted_as_phages:
        subprocess.call('cp ' + input_folder + bin_phage + file_format + ' ' + args.output + 'phage_bins/' + bin_phage + '.fasta', shell=True)
    print('**Bins predicted as phages are in the folder:', args.output + 'phage_bins/\n')
else:
    print("**We did not find any phage bins in this sample.")

# Print ending messages
print('**'+str(datetime.datetime.now()))
print('**Thank you for using Marvel!\n')
