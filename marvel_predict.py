#!/usr/bin/env python
import argparse
import re
import os
import Bio.SeqIO as SeqIO
import numpy as np
import sys
import logging
import pickle

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger('marvel')
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser('Random Forest predictor')
parser.add_argument('--min-prob', default=0.7, help='Minimum probability for acceptance [0.7]')
parser.add_argument('aa_fasta', help='Prokka protein predictions')
parser.add_argument('genbank', help='Prokka genbank')
parser.add_argument('tblout', help='Hmmsearch tblout')
args = parser.parse_args()

marvel_root = os.path.dirname(os.path.realpath(__file__))

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
                logger.warning('Prokka predicted a CDS, but there is no translation. Record ID: '.format(record.id))

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
        logger.warning('{} has 1 or zero appropriate CDS features used in prediction'.format(record.id))


num_proteins_bin = 0
with open(args.aa_fasta, 'r') as faa:
    for line in faa:
        if re.search('^>', line):
            num_proteins_bin += 1

dic_matches = {}
with open(args.tblout, 'r') as hmmsearch_out:
    for line in hmmsearch_out:
        if line.startswith('#'):
            continue
        fields = re.split(r' +', line)
        _id, _eval = fields[0], float(fields[4])
        dic_matches[_id] = _eval

# Take the proportion of proteins matching the pVOGs and store in a dictionary
i_sig = 0
for key in dic_matches:
    if dic_matches[key] <= 1e-10:
        i_sig += 1

match = re.search(r'prokka_results_(.*)\.faa', args.aa_fasta, re.IGNORECASE)
assert match, 'Amino acid fasta does not follow expected naming convention'
prefix = match.group(1)
prop_hmms_hits = {prefix: i_sig / num_proteins_bin}

data_bins = []
data_bins_index = []
# Temporary fix: Some GBL files are broken
# Verify this

# Iteration for bins
count_pred = 1
sub_data_bins = []

for record in SeqIO.parse(args.genbank, "genbank"):
    # Extract_features still take the class as an argument
    # Sending only the record now
    temp_list = extract_features(record)
    if temp_list is not None:
        sub_data_bins.append(temp_list)

if not sub_data_bins:
    print('No features found')
    sys.exit(0)

sub_data_bins_array = np.array(sub_data_bins)
mean_for_bin = list(np.mean(sub_data_bins_array, axis=0))
####
# Append here, the fifth feature: proportions of hmm hits
mean_for_bin.append(prop_hmms_hits[prefix])
data_bins.append(mean_for_bin)
data_bins_index.append(prefix)

array_bins = np.array(data_bins, dtype=float)
logger.info('Extracted {} features from {}'.format(len(sub_data_bins_array), prefix))

logger.info('Doing the machine learning prediction')

pkl_filename = os.path.join(marvel_root,
                            "models/pickle_model_rfc_trained_bins8k_refseq_all_3features_den_stran_prophitshmm.pkl")
with open(pkl_filename, 'rb') as file:
    pickle_model = pickle.load(file)

# Predict wether bins are from phage or other organisms
y_test_predicted_pickle = pickle_model.predict(array_bins[:, ])

# Assess the probability of the classification
y_test_prob = pickle_model.predict_proba(array_bins[:, ])

# report probability result
with open('predict_prob.tsv', 'w') as result_h:
    result_h.write('#BIN_ID\tPHAGE_PROB\tIS_PHAGE\n')
    result_h.write('{}\t{}\t{}\n'.format(prefix, y_test_prob[0, 1], y_test_prob[0, 1] > args.min_prob))
