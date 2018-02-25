#!/usr/bin/python3
# coding: utf-8

## This is the MARVEL Pipeline for analysis and retrieaval of Viral Long sequences
## This is a python second-half version, which receives bins as input
# Developed by Deyvid Amgarten
# Creative commons

# Import required Python modules
import numpy as np
# Biopython SeqIO
from Bio import SeqIO
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn import metrics
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
import re
import sys
import os
import subprocess
from collections import Counter
import pickle
import datetime


# Function declarations

# Usage
def usage():
    print('MARVEL Usage: marvel_bins.py -i input_folder \n\t-i input folder with bins\n\t-h This message')


# Verify arguments
def verify_arg(arg_list):
    if arg_list[1] == '-i':
        return True
    elif arg_list[1] == '-h':
        return False
    else:
        print('There is something wrong, use -h for information and usage')
        return False


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
    prefix = re.sub('.fasta', '', binn)

    # Filehandle where the output of prokka will be saved
    # output_prokka = open(str(prefix)+'prokka.output', mode='w')
    # Full command line for prokka
    command_line = (
                'prokka --kingdom Viruses --gcode 11 --cpus ' + threads + ' --force --quiet --prefix prokka_results_' + str(
            prefix) + ' --fast --norrna --notrna --outdir ' + input_folder + 'results/prokka/' + str(
            prefix) + ' --cdsrnaolap --noanno ' + input_folder + str(binn)).split()
    return_code = subprocess.call(command_line)
    # Check with prokka run smothly
    if return_code == 1:
        print("Prokka did not end well")
        quit()


# Extract features from genbank record and prepare vector of features
def extract_features(record):
    # iterate each feature
    count = 0
    sum_cds_length = 0
    strand_shift = 0
    non_coding_spacing = []
    for feature in record.features:
        # This is a modification for erroneus translations
        if feature.type == "CDS" and re.search('\w\*', str(feature.qualifiers['translation'])) is None:
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

    if len(non_coding_spacing) > 1:
        density = count / (len(record.seq) / 1000)
        mean_gene_size = sum_cds_length / count
        sum_spacing = 0
        for i in non_coding_spacing:
            sum_spacing += i
        mean_spacing_size = sum_spacing / (len(non_coding_spacing) - 1)
        ATG_freq = kmer_frequency(str(record.seq), 'ATG')
        # Return mean size of non conding regions and the density of genes
        # print(record.id, [mean_spacing_size, density, mean_gene_size, strand_shift, '/', count, flag])
        return ([density, strand_shift / count])
    else:
        warnings_handle.write("WARNING:" + record.id + " has 1 or zero appropriate CDS features.")


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

# Create Filehandle for warnings
warnings_handle = open('marvel-warnings.txt', 'w')

# Important variables
threads = '22'

# Take the input folder and list all multifasta(bins) contained inside it
# Provionally, set the input folder
input_folder = args_list[2]
list_bins = []
count_bins = 0
for j in os.listdir(input_folder):
    if re.search('.fasta', j):
        list_bins.append(j)
        count_bins += 1

print('\nWelcome to the MARVEL pipeline!\n')
print('Arguments are OK. Checked the input folder (%s) and found %d bins.\n' % (input_folder, count_bins))
print(str(datetime.datetime.now()))

# Create The results folder
try:
    os.stat(input_folder + 'results')
except:
    os.mkdir(input_folder + 'results')

#####
# PROKKA
#####
# Running prokka for all the bins multfasta files in input folder
# Perform a check in each bin, then call the execute_prokka function individually
# It may take awhile
count = 0
print('Prokka has started, this may take awhile. Be patient.\n')
for binn in list_bins:
    # Verify bin size
    len_bin = 0
    for record in SeqIO.parse(input_folder + binn, 'fasta'):
        len_bin += len(record.seq)
        # print(len_bin)
    # Make sure that bacterial bins are not take into account
    if len_bin < 500000:
        pass
        # run_prokka(binn, input_folder, threads)
    count += 1
    if count % 10 == 0:
        print('Done with %d bins...' % count)
print('Prokka tasks have finished!\n')

####
# HMM SEARCHES
####
print('Starting HMM scan, this may take awhile. Be patient.\n')
print(str(datetime.datetime.now()))
# Create a new results folder for hmmscan call
try:
    os.stat(input_folder + 'results/hmmscan/')
except:
    os.mkdir(input_folder + 'results/hmmscan/')

# Call HMMscan all bins
prop_hmms_hits = {}
i = 0
for binn in list_bins:
    prefix = re.sub('.fasta', '', binn)
    command_line_hmmscan = 'hmmscan -o ' + input_folder + 'results/hmmscan/' + prefix + '_hmmscan.out --cpu ' + threads + ' --tblout ' + input_folder + 'results/hmmscan/' + prefix + '_hmmscan.tbl --noali models/all_vogs_hmm_profiles_feb2018.hmm ' + input_folder + 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.faa'
    # In case hmmscan returns an error
    try:
        # subprocess.call(command_line_hmmscan, shell=True)
        True
    except:
        print('Error calling HMMscan:', command_line_hmmscan)
        quit()
    i += 1
    # Iteration control
    if i % 10 == 0:
        print('Done with %d bins HMM searches...' % i)
    # Parse hmmscan output files and find out which bins have less than 10% of their proteins
    # without any significant hits (10e-10)
    num_proteins_bin = 0
    with open(input_folder + 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.faa', 'rU') as faa:
        for line in faa:
            if re.search('^>', line):
                num_proteins_bin += 1
    dic_matches = {}
    with open(input_folder + 'results/hmmscan/' + prefix + '_hmmscan.tbl', 'rU') as hmmscan_out:
        for line in hmmscan_out:
            match = re.search('^VOG\d+\s+-\s+(\S+)\s+-\s+(\S+)\s+.+$', line)
            if match:
                if match.group(1) not in dic_matches:
                    dic_matches[match.group(1)] = float(match.group(2))
    # Decision about this particular bin predicted as phage buy the RF classifier
    i_sig = 0
    for key in dic_matches:
        if dic_matches[key] <= 1e-10:
            i_sig += 1
    prop_hmms_hits[prefix] = i_sig / num_proteins_bin
    # if i_sig/num_proteins_bin > 0.10:
    #    bins_passed_hmmscan.append(bin_phage)
    #    print(i_sig/num_proteins_bin, bin_phage)

# print(prop_hmms_hits)
# quit()
###
# Machine learning prediction of viral bins
###

# Go for each bin and extract relevant features from respective genbank files
# Final feature vector is "array_bins"
print('Extracting features from bins...\n')
data_bins = []
data_bins_index = []
# Temporary fix: Some GBL files are broken

# Iteration for bins
i = 0
for bins in list_bins:
    i += 1
    prefix = re.sub('.fasta', '', bins)

    sub_data_bins = []
    # GBK File with gene prediction
    file_name = input_folder + 'results/prokka/' + prefix + '/prokka_results_' + prefix + '.gbk'
    for record in SeqIO.parse(file_name, "genbank"):
        # Extract_features still take the class as an argument
        # Sending only the record now
        temp_list = extract_features(record)
        if temp_list is None:
            pass
        else:
            sub_data_bins.append(temp_list)
    ###
    # Iteration control
    # if i%100 == 0:
    #    print('Done with', i, 'bins')
    # Empty list
    if not sub_data_bins:
        continue
    sub_data_bins_array = np.array(sub_data_bins)
    # print(sub_data_caudo_array)
    mean_for_bin = list(np.mean(sub_data_bins_array, axis=0))
    # print(mean_for_bin)
    ####
    # Append here, the fifth feature: proportions of hmm hits
    mean_for_bin.append(prop_hmms_hits[prefix])
    data_bins.append(mean_for_bin)
    data_bins_index.append(prefix)
array_bins = np.array(data_bins, dtype=float)
print('Extracted features from', len(array_bins), 'bins\n')

# Load RFC from file
print('Doing the machine learning prediction...\n')
pkl_filename = "models/pickle_model_rfc_trained_bins_refseq_until2015_3features_den_stran_prophitshmm.pkl"
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
for pred in y_test_prob[:, 1]:
    if pred >= 0.7:
        bins_predicted_as_phages.append(data_bins_index[i])
        # print(pred, data_bins_index[i], array_bins[i])
    i += 1

### Alternative way: Not using probability
# bins_predicted_as_phages = []
# i = 0
# for pred in y_test_predicted_pickle:
#    if pred == 1:
#        bins_predicted_as_phages.append(data_bins_index[i])
#    i += 1


try:
    os.stat(input_folder + 'results/phage_bins')
except:
    os.mkdir(input_folder + 'results/phage_bins')
for bin_phage in bins_predicted_as_phages:
    # temp_file_name = input_folder+'prokka_results_'+bin_phage+'/prokka_results_'+bin_phage+'.faa'
    subprocess.call(
        'cp ' + input_folder + bin_phage + '.fasta ' + input_folder + 'results/phage_bins/' + bin_phage + '.fasta',
        shell=True)

print('Finished Machine learning predictions!\n')

# Just make sure to end the program closing the warnings filehandle
warnings_handle.close()

# Prin ending messages
print('Bins predicted as phages are in the folder:', input_folder + 'results/phage_bins/\n')
print('Thank you for using Marvel!\n')

