#!/usr/bin/python3
# coding: utf-8

## This is the MARVEL Pipeline for analysis and retrieaval of Viral Long sequences
## This is an auxiliary script to download and set models
# Developed by Deyvid Amgarten
# Creative commons

# Greeting message
print('\nYou only need to run this script once!\n')

# Verify databases
if not os.path.isfile('models/all_vogs_hmm_profiles_feb2018.hmm'):
    print('Downloading flat file database. Do not worry, that will just take a few minutes and is executed only in the first time... \n')
    os.system('wget ftp://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz')
    print('Extracting database file...\n')
    subprocess.run('tar -xzf AllvogHMMprofiles.tar.gz', shell=True)
    subprocess.run('cat AllvogHMMprofiles/* > models/all_vogs_hmm_profiles_feb2018.hmm', shell=True)
    subprocess.run('rm -r AllvogHMMprofiles/ AllvogHMMprofiles.tar.gz', shell=True)
    print('Compressing hmm database...')
    subprocess.run('hmmpress models/all_vogs_hmm_profiles_feb2018.hmm', shell=True)
    print('Database is all set!\n')
else:
    print('HMM Database is already set.\n')

if not os.stat('example_data/'):
    print('Downloading example datasets. This will be fast\n')
    os.system('wget http://projetos.lbi.iq.usp.br/metazoo/deyvid/datasets/example_data.tar.gz')
    subprocess.run('tar -xzf example_data.tar.gz', shell=True)
    subprocess.run('rm example_data.tar.gz',shell=True)
else:
    print('Example datasets are already in your disk.\n')

print('All set. Thank you for using MARVEL.')