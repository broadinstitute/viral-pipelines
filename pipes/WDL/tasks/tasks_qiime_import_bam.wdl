version 1.0

task qiime_import_bam {
    meta {
        description: "Parsing demultiplexed fastq BAM files into qiime readable files."
    }
    input { 
        File    reads_bam
        String  sample_name
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(reads_bam, "GiB")) + 5
        String  docker     = "quay.io/qiime2/core:2022.8" 
    }
    parameter_meta {
        reads_bam: {description: "Input BAM file"}
    }

    command <<<
        set -ex -o pipefail
        #testing to see if source is correct
        # obtain the name of the qiime conda environment in the container flexible about version
        CONDA_ENV_NAME=$(conda info --envs -q | awk -F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        #https://github.com/chanzuckerberg/miniwdl/issues/603
        #manual activation of bash per suggestion: https://github.com/conda/conda/issues/7980
        conda activate ${CONDA_ENV_NAME}

        #Part 1A | BAM -> FASTQ [Simple samtools command]
        samtools fastq -1 $(pwd)/R1.fastq.gz -2 $(pwd)/R2.fastq.gz -0 /dev/null ~{reads_bam}
        #making new bash variable | regex: (_) -> (-)
        NEWSAMPLENAME=$(echo "~{sample_name}" | perl -lape 's/[_]/-/g')
        #All names added to one giant file 
        echo ${NEWSAMPLENAME} > NEWSAMPLENAME.txt
        #Make a manifest.txt that contains [1.sampleid 2.R1_fastq 3.R2.fastq]
        #> =overwrite or writes new file 
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
        #>>= appends 
        #\t= tabs next value 
        echo -e "$NEWSAMPLENAME\t$(pwd)/R1.fastq.gz\t$(pwd)/R2.fastq.gz" >> manifest.tsv
        
        #fastq -> bam (provided by qiime tools import fxn)
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path manifest.tsv \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path "~{sample_name}.qza"
    >>>

    output {
        File   reads_qza               = "~{sample_name}.qza"
        String samplename_master_sheet = read_string("NEWSAMPLENAME.txt")
    }
    runtime {
        docker: docker
        memory: "${memory_mb} MiB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}