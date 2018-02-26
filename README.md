
# MARVEL - Metagenomic Analysis and Retrieval of Viral Extended Sequences

MARVEL is a pipeline for recovery of complete phage genomes from whole community shotgun metagenomic sequencing data.
The tool was divided in three components:
   * **generate_bins_from_reads.py** - Generate bins from whole community shotgun reads
   * **marvel_bins.py** - Machine learning prediction of phage bins
   * **function_driven_analyses.py** - Analyze specific genes of interest within phage bins


### Dependencies

MARVEL's main scrip (marvel_bins.py) require Prokka and its dependencies to be installed:

* [Prokka](https://github.com/tseemann/prokka) - Rapid prokaryotic genome annotation.

Also, these Python 3 libraries need to be installed:

* [Biopython](http://biopython.org/) - Handling biological sequences and records
* [Scikit Learn](http://scikit-learn.org/stable/) - Machine Learning prediction

```
pip install numpy,scypy,biopython
pip install -U scikit-learn
```
Alternatively, you may want to install the [Conda](https://anaconda.org/) package manager and just run the commands above using 'conda' instead of 'pip'.


### Installing

Getting MARVEL ready to run is as simple as clone this Github project or dowloading it to a directory inside you computer:

```
git clone https://github.com/deyvidamgarten/MARVEL
```

### Getting Started

Inside the directory where MARVEL was extracted (or cloned), just run:

```
python3 marvel_bins.py -i input_directory -t num_threads
```

Change 'input_directory' for the folder where bins in fasta format are stored and 'num_threads' for the number of CPUs core to be used.

Results are stored in the 'Results' folder inside the input directory.  

### Running the example datasets

We provide a folder with example datasets containing mocking bins of RefSeq viral and bacterial genomes.

To try the examples, just run:

```
python3 marvel_bins.py -i example_data/ -t 12
```

### Additional scripts

MARVEL's main script receives metagenomic bins as input. However, we additionally provide a simple scrip which receives
metagenomic reads (Illumina sequencing), assemblies the data in contigs and generates bins.
[metaSpades](http://bioinf.spbau.ru/spades) and [Metabat2](https://bitbucket.org/berkeleylab/metabat) are used for assembling and binning, respectively.

We can't stress enough that there are several tools for assembly and binning, which should be well-chosen according to
the researcher's purposes. Here, our intention is only facilitate the use of our tool.  

```
python3 generate_bins_from_reads.py -1 reads-R1.fastq -2 reads-R2.fastq -t num_threads
```

### Authors




### License

This project is licensed under the Creative Commons License. Codes here may be used for any purposed or modified.


