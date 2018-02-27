
# MARVEL - Metagenomic Analysis and Retrieval of Viral Extended Sequences

MARVEL is a pipeline for recovery of complete phage genomes from whole community shotgun metagenomic sequencing data.

Main script:
   * **marvel_bins.py** - Machine learning prediction of phage bins
  
Auxiliary scripts:
   * **generate_bins_from_reads.py** - Generates metagenomic bins given Illumina sequencing reads
   * **function_driven_analyses.py** - Analyze specific genes of interest within phage bins


### Dependencies

MARVEL's main scrip (marvel_bins.py) require Prokka and its dependencies to be installed:

* [Prokka](https://github.com/tseemann/prokka) - Rapid Prokaryotic genome annotation.

These Python libraries are also required:

* [Biopython](http://biopython.org/) - Handling biological sequences and records
* [Scikit Learn](http://scikit-learn.org/stable/) - Machine Learning prediction

```
pip install -U numpy,scypy,biopython
pip install scikit-learn
```
Alternatively, you may want to install the [Conda](https://anaconda.org/) package manager and just run the commands above using 'conda' instead of 'pip'.


### Installing

Getting MARVEL ready to run is as simple as clone this Github project or dowloading it to a directory inside you computer:

```
git clone https://github.com/deyvidamgarten/MARVEL
```

### Getting Started

Inside the directory where MARVEL was extracted (or cloned), you will need to download and set the models. 
This is required only once and it is simple. Just run:
```
python download_and_set_models.py
```
All set! 

Now, to run MARVEL type:
```
python marvel_bins.py -i input_directory -t num_threads
```

Change 'input_directory' for the folder where bins are stored in fasta format and 'num_threads' for the number of CPU cores to be used.

Results will be stored in the 'Results' folder inside the input directory.  

### Running the example datasets

We provide a folder with example datasets containing mocking bins of RefSeq viral and bacterial genomes.

To try these examples, run:

```
python marvel_bins.py -i example_data/bins_8k_refseq -t 12
```

### Additional scripts

MARVEL's main script receives metagenomic bins as input. However, we additionally provide a simple scrip which receives
metagenomic reads (Illumina sequencing) and generates bins.
[metaSpades](http://bioinf.spbau.ru/spades) and [Metabat2](https://bitbucket.org/berkeleylab/metabat) are used for assembling and binning, respectively.

We can't stress enough that there are several tools for assembly and binning, which should be well-chosen according to
the researcher's purposes. Here, our intention is only facilitate the use of our tool.  

```
python generate_bins_from_reads.py -1 reads-R1.fastq -2 reads-R2.fastq -t num_threads
```

### Authors




### License

This project is licensed under the Creative Commons License. Codes here may be used for any purposed or modified.


