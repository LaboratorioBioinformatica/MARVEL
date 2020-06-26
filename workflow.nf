import groovy.util.FileNameByRegexFinder

minimumExtent = 2000

hmmDb = Channel.fromPath(params.hmm_db)

fnf = new FileNameByRegexFinder()

// search the given path for fasta files
fastaFiles = fnf.getFileNames(params.bin_dir, '.*\\.(?:fna|fasta|fa)$')
        .collect{
            n = 0
            f = new File(it)
            f.toPath().eachLine{ if (!(it =~ /^>/)) { n += it.size() } }
            [fasta: f, length: n]
        }

// remove bins that are too small
fastaFiles.removeAll{it.length <= minimumExtent}

println("Accepted ${fastaFiles.size()} bins for analysis")

input_fasta = Channel.fromList(fastaFiles)

process Prokka {
        cpus 1
        publishDir params.outdir, mode: 'copy', overwrite: true, pattern: "{*.faa,*.gbk,*.gff}", saveAs: {fn -> "${prefix}/${fn}"}

        input:
        tuple nt_fasta, extent from input_fasta

        output:
        tuple prefix, nt_fasta, path("*.faa"), path("*.gbk") into prokka_files

        script:
        prefix = nt_fasta.baseName - ~/\\.\w+$/

        """
        prokka  --cpus ${task.cpus} --kingdom Viruses --centre X --compliant --gcode 11 \
                --fast --norrna --notrna --cdsrnaolap --noanno \
                --force --quiet --prefix prokka_results_$prefix --outdir . $nt_fasta
        """
}

process Hmmsearch {
        cpus 2
        publishDir params.outdir, mode: 'copy', overwrite: true, pattern: "*.tbl", saveAs: {fn -> "${prefix}/${fn}"}

        input:
        tuple prefix, nt_fasta, aa_fasta, genbank from prokka_files
        each path(hmm_db) from hmmDb

        output:
        tuple prefix, nt_fasta, aa_fasta, genbank, path("*.tbl") into hmmer_files

        """
        hmmsearch -o ${prefix}_hmmsearch.out --cpu ${task.cpus} --tblout ${prefix}_hmmsearch.tbl \
                  --noali  $hmm_db $aa_fasta
        """
}

process Predict {

    input:
    tuple prefix, nt_fasta, aa_fasta, genbank, tblout from hmmer_files

    output:
    path('predict_prob.tsv') into predict_out

    """
    marvel_predict.py $aa_fasta $genbank $tblout
    """
}

predict_out.collectFile(
        name: 'phage_predictions.tsv',
        newLine: false,
        storeDir: params.outdir,
        keepHeader: true,
        sort: true)
        .subscribe{
            println "Entries are saved to file: $it"
        }
