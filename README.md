
<p align="center"> <img src="logo_marvel.png" height="110" alt="MARVEL" /> </p>

# Metagenomic Analysis and Retrieval of Viral Elements

MARVEL is a tool for recovery of draft phage genomes from whole community shotgun metagenomic sequencing data.  

Main script:
   * **marvel_bins.py** - Machine learning prediction of phage bins
  
Auxiliary script:
   * **generate_bins_from_reads.py** - Generates metagenomic bins, given Illumina sequencing reads

## Reference and citation

A manuscript describing MARVEL was published in [Frontiers in Genetics](https://www.frontiersin.org/articles/10.3389/fgene.2018.00304/full).  
If you find MARVEL useful in your research, please cite:  
*Amgarten DE, Braga LP, Da Silva AM, Setubal JC. MARVEL, a Tool for Prediction of Bacteriophage Sequences in Metagenomic Bins. Frontiers in Genetics. 2018;9:304.*

## New Viral Groups

We are working to train models to new viral groups. Let us know if a particular group would be helpful in your research.

## Dependencies

All scripts from this project were coded in [Python 3](https://www.python.org/). So, first of all, make sure you have it installed and updated.  
MARVEL's main scrip (marvel_bins.py) requires Prokka and HMMER tools as dependencies. By installing Prokka and its dependencies, you will usually  install HMMER tools automatically.  
**It seems that prokka has changed one of its output file's extension from ".gbk" to ".gbf". If you enconter a "File not found" error, please let me know I will fix it asap (deyvid.amgarten@usp.br).**

* [Prokka](https://github.com/tseemann/prokka) - Rapid Prokaryotic genome annotation.
* [HMMER Tools](http://www.hmmer.org/) - Biosequence analysis using profile hidden Markov models

These Python libraries are required:

* [Numpy](http://www.numpy.org/), [Scipy](https://www.scipy.org/) - Efficiently handling arrays and scientific computing
* [Biopython](http://biopython.org/) - Handling biological sequences and records
* [Scikit-learn](http://scikit-learn.org/stable/) - Machine learning

To install these Python libraries, just type: 
```
pip3 install -U numpy scipy biopython scikit-learn
```

## Installing

Getting MARVEL ready to run is as simple as cloning this Github project or dowload and extract it to a directory inside your computer:

```
git clone https://github.com/LaboratorioBioinformatica/MARVEL
```

## Quick start

Inside the directory where MARVEL was extracted (or cloned), you will need to download and set the models. 
This is required only once and it is simple. Just run:
```
python3 download_and_set_models.py
```
All set!  
Now, to run MARVEL type:
```
python3 marvel_bins.py -i input_directory -t num_threads
```

Change 'input_directory' to the folder where bins are stored in fasta format and 'num_threads' to the number of CPU cores to be used. Several threads should be used to speed up prokka and hmm searches.  
Results will be stored in the 'Results' folder inside the input directory.  
Obs: You need to execute the scripts from the directory where MARVEL was extracted, i.e., MARVEL's root folder. 

## Running the example datasets

We provide a folder with example datasets containing mocking bins of RefSeq viral and bacterial genomes.  
To try these examples, run:

```
python3 marvel_bins.py -i example_data/bins_8k_refseq -t 12
```

## Additional scripts

MARVEL's main script receives metagenomic bins as input. However, we additionally provide a simple scrip which receives
metagenomic reads (Illumina sequencing) and generates bins.
[metaSpades](http://bioinf.spbau.ru/spades), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Metabat2](https://bitbucket.org/berkeleylab/metabat) are used for assembling, mapping and binning, respectively.  

In case you want to generate the bins by yourself, we recommend two special parameters in the binning process: -m 1500 -s 10000. These parameters tell metabat to generate bins with contigs of at least 1500 bp and with a minimum total size of 10 kbp, which makes more sense when one is trying to retrieve viral genomes.  
We can't stress enough that there are several tools for assembly and binning, which should be well-chosen according to
the researcher's purposes. Our intention here is to facilitate the use of our tool.  

```
python3 generate_bins_from_reads.py -1 reads-R1.fastq -2 reads-R2.fastq -t num_threads
```
## Suggested workflow for improved genome recovery

MARVEL is indicated to metagenomic studies where whole community DNA sequencing reads are available. The process of binning with Metabat2 will work better if you have several samples from a same environment (time-series samples for example), so Metabat2 can assess contigs' abundancy correlation among samples and create better bins. Here we suggest a workflow of analyses to generate draft phage genomes from time-series samples raw sequencing reads.

1. Assembly raw reads for each sample individually with [metaSpades](http://bioinf.spbau.ru/spades).

2. Join contigs in a single multi-FASTA and use dedupe from [BBMap](https://github.com/BioInfoTools/BBMap) tools to remove duplicated contigs.

3. Use resulting deduplicated contigs as reference for mapping samples' reads individually with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

3. Generate bins with metabat2 using BAM files obtained in step 3. At this point, it is important to set these specific parameters in metabat2: -m 1500 -s 10000. You should use these parameters, so metabat2 can generate bins with phage genomes characteristics.

4. Run MARVEL giving bins folder as input. Bins predicted as phage will be in the folder: results/phage_genomes/.

For improved draft genomes, we suggested the following additional steps:

6. Merge contigs with overlapping ends with [Phrap](http://www.phrap.org/phredphrapconsed.html) for each individual bin predicted by MARVEL as phage.

7. Further validate predicted bins by assessing bacterial/archaeal genes with [CheckM](https://github.com/Ecogenomics/CheckM/wiki) and predicting tail/capsid proteins with [VIRALpro](http://scratch.proteomics.ics.uci.edu/explanation.html#VIRALpro).

If you have any question, please contact us. We will be glad to help with your analyses.  
deyvid.amgarten@usp.br

## Simulated RefSeq datasets for training and testing

All the simulated datasets used for training and testing the Random Forest classifier, as well as predicted bins from composting samples were made available through this link:

[Browse and Download datasets](http://projetos.lbi.iq.usp.br/metazoo/deyvid/datasets/) 

## Author
[Deyvid Amgarten](https://sites.google.com/view/deyvid/english)  
This tool was developed as part of my PhD thesis by the [Bioinformatics Graduate Program](https://www.ime.usp.br/en/bioinformatics/graduate) from the University of Sao Paulo, Brazil.

## Setulab
Our group is interested in studying:

* Classical machine learning techniques to create predictors, which are being applied in genomics and metagenomics;
* Deep Learning techniques applied to more broad problems, as for instance prediction of complex biological attributes in viruses and bacteria.

Please visit [Setulab's page](http://lbi.usp.br/learning/) for more information.


## License

This project is licensed under GNU license. Codes here may be used for any purposed or modified.  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

