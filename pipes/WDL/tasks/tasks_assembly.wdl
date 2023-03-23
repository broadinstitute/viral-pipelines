version 1.0

task assemble {
    input {
      File     reads_unmapped_bam
      File     trim_clip_db
      
      Int      spades_n_reads = 10000000
      Int?     spades_min_contig_len
      String?  spades_options
      
      Boolean  always_succeed = false
      
      # do this in two steps in case the input doesn't actually have "taxfilt" in the name
      String   sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt")
      
      Int?     machine_mem_gb
      String   docker = "quay.io/broadinstitute/viral-assemble:2.1.33.0"
    }
    parameter_meta{
      reads_unmapped_bam: {
        description: "Unaligned reads in BAM format."
        patterns: ["*.bam"]
        common: "required"
      }
      trim_clip_db: {
        description: ""
        patterns
      }
      

    }

    Int disk_size = 375

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
        mem_in_gb=$(/opt/viral-ngs/source/docker/calc_mem.py gb 90)

        assembly.py --version | tee VERSION

        assembly.py assemble_spades \
          ~{reads_unmapped_bam} \
          ~{trim_clip_db} \
          ~{sample_name}.assembly1-spades.fasta \
          ~{'--nReads=' + spades_n_reads} \
          ~{true="--alwaysSucceed" false="" always_succeed} \
          ~{'--minContigLen=' + spades_min_contig_len} \
          ~{'--spadesOpts="' + spades_options + '"'} \
          --memLimitGb $mem_in_gb \
          --outReads=~{sample_name}.subsamp.bam \
          --loglevel=DEBUG

        samtools view -c ~{sample_name}.subsamp.bam | tee subsample_read_count >&2
    }

    output {
        File   contigs_fasta        = "~{sample_name}.assembly1-spades.fasta"
        File   subsampBam           = "~{sample_name}.subsamp.bam"
        Int    subsample_read_count = read_int("subsample_read_count")
        String viralngs_version     = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 15]) + " GB"
        cpu: 4
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x8"
        maxRetries: 2
    }

}

task scaffold {
    input {
      File         contigs_fasta
      File         reads_bam
      Array[File]+ reference_genome_fasta

      String       aligner="muscle"
      Float?       min_length_fraction
      Float?       min_unambig
      Int          replace_length=55

      Int?         nucmer_max_gap
      Int?         nucmer_min_match
      Int?         nucmer_min_cluster
      Int?         scaffold_min_contig_len
      Float?       scaffold_min_pct_contig_aligned

      Int?         machine_mem_gb
      String       docker="quay.io/broadinstitute/viral-assemble:2.1.33.0"

      # do this in multiple steps in case the input doesn't actually have "assembly1-x" in the name
      String       sample_name = basename(basename(contigs_fasta, ".fasta"), ".assembly1-spades")
    }

    Int disk_size = 375

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_gb=$(/opt/viral-ngs/source/docker/calc_mem.py gb 90)

        assembly.py --version | tee VERSION

        assembly.py order_and_orient \
          ~{contigs_fasta} \
          ~{sep=' ' reference_genome_fasta} \
          ~{sample_name}.intermediate_scaffold.fasta \
          ~{'--min_contig_len=' + scaffold_min_contig_len} \
          ~{'--maxgap=' + nucmer_max_gap} \
          ~{'--minmatch=' + nucmer_min_match} \
          ~{'--mincluster=' + nucmer_min_cluster} \
          ~{'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned} \
          --outReference ~{sample_name}.scaffolding_chosen_ref.fasta \
          --outStats ~{sample_name}.scaffolding_stats.txt \
          --outAlternateContigs ~{sample_name}.scaffolding_alt_contigs.fasta \
          --loglevel=DEBUG

        grep '^>' ~{sample_name}.scaffolding_chosen_ref.fasta | cut -c 2- | cut -f 1 -d ' ' > ~{sample_name}.scaffolding_chosen_refs.txt

        assembly.py gapfill_gap2seq \
          ~{sample_name}.intermediate_scaffold.fasta \
          ~{reads_bam} \
          ~{sample_name}.intermediate_gapfill.fasta \
          --memLimitGb $mem_in_gb \
          --maskErrors \
          --loglevel=DEBUG

        grep -v '^>' ~{sample_name}.intermediate_gapfill.fasta | tr -d '\n' | wc -c | tee assembly_preimpute_length
        grep -v '^>' ~{sample_name}.intermediate_gapfill.fasta | tr -d '\nNn' | wc -c | tee assembly_preimpute_length_unambiguous

        assembly.py impute_from_reference \
          ~{sample_name}.intermediate_gapfill.fasta \
          ~{sample_name}.scaffolding_chosen_ref.fasta \
          ~{sample_name}.scaffolded_imputed.fasta \
          --newName ~{sample_name} \
          ~{'--replaceLength=' + replace_length} \
          ~{'--minLengthFraction=' + min_length_fraction} \
          ~{'--minUnambig=' + min_unambig} \
          ~{'--aligner=' + aligner} \
          --loglevel=DEBUG
    }

    output {
        File   scaffold_fasta                        = "~{sample_name}.scaffolded_imputed.fasta"
        File   intermediate_scaffold_fasta           = "~{sample_name}.intermediate_scaffold.fasta"
        File   intermediate_gapfill_fasta            = "~{sample_name}.intermediate_gapfill.fasta"
        Int    assembly_preimpute_length             = read_int("assembly_preimpute_length")
        Int    assembly_preimpute_length_unambiguous = read_int("assembly_preimpute_length_unambiguous")
        Array[String] scaffolding_chosen_ref_names   = read_lines("~{sample_name}.scaffolding_chosen_refs.txt")
        File   scaffolding_chosen_ref                = "~{sample_name}.scaffolding_chosen_ref.fasta"
        File   scaffolding_stats                     = "~{sample_name}.scaffolding_stats.txt"
        File   scaffolding_alt_contigs               = "~{sample_name}.scaffolding_alt_contigs.fasta"
        String viralngs_version                      = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 31]) + " GB"
        cpu: 4
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x8"
        maxRetries: 2
    }
}

task ivar_trim {
    meta {
      description: "this runs ivar trim on aligned reads, which results in soft-clipping of alignments"
    }

    input {
      File   aligned_bam
      File?  trim_coords_bed
      Int?   min_keep_length
      Int?   sliding_window
      Int?   min_quality = 1
      Int?   primer_offset
      
      Int?   machine_mem_gb
      String docker = "andersenlabapps/ivar:1.3.1"
    }

    String  bam_basename=basename(aligned_bam, ".bam")
    Int disk_size = 375

    parameter_meta {
      aligned_bam:     { description: "aligned reads in BAM format", patterns: ["*.bam"] }
      trim_coords_bed: { description: "optional primers to trim in reference coordinate space (0-based BED format)", patterns: ["*.bed"] }
      min_keep_length: { description: "Minimum length of read to retain after trimming (Default: 30)" }
      sliding_window:  { description: "Width of sliding window for quality trimming (Default: 4)" }
      min_quality:     { description: "Minimum quality threshold for sliding window to pass (Default: 20)" }
    }

    command {
        ivar version | head -1 | tee VERSION
        if [ -f "~{trim_coords_bed}" ]; then
          ivar trim -e \
            ~{'-b ' + trim_coords_bed} \
            ~{'-m ' + min_keep_length} \
            ~{'-s ' + sliding_window} \
            ~{'-q ' + min_quality} \
            ~{'-x ' + primer_offset} \
            -i ~{aligned_bam} -p trim | tee IVAR_OUT
          samtools sort -@ $(nproc) -m 1000M -o ~{bam_basename}.trimmed.bam trim.bam
        else
          echo "skipping ivar trim"
          cp "~{aligned_bam}" "~{bam_basename}.trimmed.bam"
          echo "Trimmed primers from 0% (0) of reads." > IVAR_OUT
        fi
        PCT=$(grep "Trimmed primers from" IVAR_OUT | perl -lape 's/Trimmed primers from (\S+)%.*/$1/')
        if [[ $PCT = -* ]]; then echo 0; else echo $PCT; fi > IVAR_TRIM_PCT
        grep "Trimmed primers from" IVAR_OUT | perl -lape 's/Trimmed primers from \S+% \((\d+)\).*/$1/' > IVAR_TRIM_COUNT
    }

    output {
        File   aligned_trimmed_bam         = "~{bam_basename}.trimmed.bam"
        Float  primer_trimmed_read_percent = read_float("IVAR_TRIM_PCT")
        Int    primer_trimmed_read_count   = read_int("IVAR_TRIM_COUNT")
        String ivar_version                = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 4
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x4"
        maxRetries: 2
    }
}

task ivar_trim_stats {

    input {
      File   ivar_trim_stats_tsv
      String out_basename = "ivar_trim_stats"
      String flowcell = ""

      String docker = "quay.io/broadinstitute/py3-bio:0.1.2"
    }

    Int disk_size = 50

    command <<<
      set -e
      python3<<CODE

      import json
      import pandas as pd
      import plotly.express as px

      # load and clean up data
      df = pd.read_csv("~{ivar_trim_stats_tsv}", delimiter='\t',
        names=['file', 'trim_percent', 'trim_count'])

      # make plot
      flowcell = "~{flowcell}"
      title = "ivar trim: % vs # of reads trimmed per sample"
      if flowcell:
        title += " ({})".format(flowcell)
      p = px.scatter(df,
        x='trim_count', y='trim_percent',
        title=title,
        opacity=0.7,
        hover_data=df.columns)

      # export
      out_basename = "~{out_basename}"
      df.to_csv(out_basename + ".txt", sep='\t')
      p.write_html(out_basename + ".html")
      p.write_image(out_basename + ".png")

      CODE
    >>>

    output {
      File trim_stats_html = "~{out_basename}.html"
      File trim_stats_png  = "~{out_basename}.png"
      File trim_stats_tsv  = "~{out_basename}.txt"
    }

    runtime {
        docker: docker
        memory: "1 GB"
        cpu: 1
        disks:   "local-disk " + disk_size + " HDD"
        disk:    disk_size + " GB"
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
}

task align_reads {
  meta {
    description: "Align unmapped reads to a reference genome, either using novoalign (default), minimap2, or bwa. Produces an aligned bam file (including all unmapped reads), an aligned-only bam file, both sorted and indexed, along with samtools flagstat output, fastqc stats (on mapped only reads), and some basic figures of merit."
  }

  input {
    File     reference_fasta
    File     reads_unmapped_bam

    File?    novocraft_license

    String   aligner = "minimap2"
    String?  aligner_options
    Boolean? skip_mark_dupes = false

    Int?     machine_mem_gb
    String   docker = "quay.io/broadinstitute/viral-core:2.1.33"

    String   sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")
  }

  Int disk_size = 375

  parameter_meta {
    aligner: { description: "Short read aligner to use: novoalign, minimap2, or bwa. (Default: novoalign)" }
  }
  
  command <<<
    set -ex # do not set pipefail, since grep exits 1 if it can't find the pattern

    read_utils.py --version | tee VERSION

    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    cp "~{reference_fasta}" assembly.fasta
    grep -v '^>' assembly.fasta | tr -d '\n' | wc -c | tee assembly_length

    if [ "$(cat assembly_length)" != "0" ]; then

      # only perform the following if the reference is non-empty

      if [ "~{aligner}" == "novoalign" ]; then
        read_utils.py novoindex \
          assembly.fasta \
          ~{"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG
      fi
      read_utils.py index_fasta_picard assembly.fasta --loglevel=DEBUG
      read_utils.py index_fasta_samtools assembly.fasta --loglevel=DEBUG

      read_utils.py align_and_fix \
        "~{reads_unmapped_bam}" \
        assembly.fasta \
        --outBamAll "~{sample_name}.all.bam" \
        --outBamFiltered "~{sample_name}.mapped.bam" \
        --aligner ~{aligner} \
        ~{'--aligner_options "' + aligner_options + '"'} \
        ~{true='--skipMarkDupes' false="" skip_mark_dupes} \
        --JVMmemory "$mem_in_mb"m \
        ~{"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

    else
      # handle special case of empty reference fasta -- emit empty bams (with original bam headers)
      samtools view -H -b "~{reads_unmapped_bam}" > "~{sample_name}.all.bam"
      samtools view -H -b "~{reads_unmapped_bam}" > "~{sample_name}.mapped.bam"

      samtools index "~{sample_name}.all.bam" "~{sample_name}.all.bai"
      samtools index "~{sample_name}.mapped.bam" "~{sample_name}.mapped.bai"
    fi

    cat /proc/loadavg > CPU_LOAD

    # collect figures of merit
    grep -v '^>' assembly.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    samtools view -c "~{reads_unmapped_bam}" | tee reads_provided
    samtools view -c "~{sample_name}.mapped.bam" | tee reads_aligned
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 "~{sample_name}.all.bam" | samtools flagstat - | tee ~{sample_name}.all.bam.flagstat.txt
    grep properly "~{sample_name}.all.bam.flagstat.txt" | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view "~{sample_name}.mapped.bam" | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("$(cat bases_aligned)")/"$(cat assembly_length_unambiguous)") if "$(cat assembly_length_unambiguous)">0 else print(0)" > mean_coverage

    # fastqc mapped bam
    reports.py fastqc ~{sample_name}.mapped.bam ~{sample_name}.mapped_fastqc.html --out_zip ~{sample_name}.mapped_fastqc.zip

    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    { cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes || echo 0; } > MEM_BYTES
  >>>

  output {
    File   aligned_bam                   = "~{sample_name}.all.bam"
    File   aligned_bam_idx               = "~{sample_name}.all.bai"
    File   aligned_bam_flagstat          = "~{sample_name}.all.bam.flagstat.txt"
    File   aligned_only_reads_bam        = "~{sample_name}.mapped.bam"
    File   aligned_only_reads_bam_idx    = "~{sample_name}.mapped.bai"
    File   aligned_only_reads_fastqc     = "~{sample_name}.mapped_fastqc.html"
    File   aligned_only_reads_fastqc_zip = "~{sample_name}.mapped_fastqc.zip"
    Int    reads_provided                = read_int("reads_provided")
    Int    reads_aligned                 = read_int("reads_aligned")
    Int    read_pairs_aligned            = read_int("read_pairs_aligned")
    Float  bases_aligned                 = read_float("bases_aligned")
    Float  mean_coverage                 = read_float("mean_coverage")
    Int    max_ram_gb                    = ceil(read_float("MEM_BYTES")/1000000000)
    Int    runtime_sec                   = ceil(read_float("UPTIME_SEC"))
    String cpu_load                      = read_string("CPU_LOAD")
    String viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    cpu: 8
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x8"
    preemptible: 1
    maxRetries: 2
  }
}

task refine_assembly_with_aligned_reads {
    meta {
      description: "This step refines/polishes a genome based on short read alignments, producing a new consensus genome. Uses GATK3 Unified Genotyper to produce new consensus. Produces new genome (fasta), variant calls (VCF), and figures of merit."
    }

    input {
      File     reference_fasta
      File     reads_aligned_bam
      String   sample_name

      Boolean  mark_duplicates = false
      Float    major_cutoff = 0.5
      Int      min_coverage = 3

      Int?     machine_mem_gb
      String   docker = "quay.io/broadinstitute/viral-assemble:2.1.33.0"
    }

    Int disk_size = 375

    parameter_meta {
      major_cutoff: {
        description: "If the major allele is present at a frequency higher than this cutoff, we will call an unambiguous base at that position.  If it is equal to or below this cutoff, we will call an ambiguous base representing all possible alleles at that position."
      }
      min_coverage: {
        description: "Minimum read coverage required to call a position unambiguous."
      }
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        assembly.py --version | tee VERSION

        if [ ~{true='true' false='false' mark_duplicates} == "true" ]; then
          read_utils.py mkdup_picard \
            ~{reads_aligned_bam} \
            temp_markdup.bam \
            --JVMmemory "$mem_in_mb"m \
            --loglevel=DEBUG
        else
          ln -s ~{reads_aligned_bam} temp_markdup.bam
        fi
        samtools index -@ $(nproc) temp_markdup.bam temp_markdup.bai

        ln -s ~{reference_fasta} assembly.fasta
        assembly.py refine_assembly \
          assembly.fasta \
          temp_markdup.bam \
          refined.fasta \
          --already_realigned_bam=temp_markdup.bam \
          --outVcf ~{sample_name}.sites.vcf.gz \
          --min_coverage ~{min_coverage} \
          --major_cutoff ~{major_cutoff} \
          --JVMmemory "$mem_in_mb"m \
          --loglevel=DEBUG

        file_utils.py rename_fasta_sequences \
          refined.fasta "${sample_name}.fasta" "${sample_name}"

        # collect variant counts
        if (( $(cat refined.fasta | wc -l) > 1 )); then
          bcftools filter -e "FMT/DP<${min_coverage}" -S . "${sample_name}.sites.vcf.gz" -Ou | bcftools filter -i "AC>1" -Ou > "${sample_name}.diffs.vcf"
          bcftools filter -i 'TYPE="snp"'  "${sample_name}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_snps
          bcftools filter -i 'TYPE!="snp"' "${sample_name}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_indels
        else
          # empty output
          echo "0" > num_snps
          echo "0" > num_indels
          cp "${sample_name}.sites.vcf.gz" "${sample_name}.diffs.vcf"
        fi

        # collect figures of merit
        set +o pipefail # grep will exit 1 if it fails to find the pattern
        grep -v '^>' refined.fasta | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' refined.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    }

    output {
        File   refined_assembly_fasta      = "${sample_name}.fasta"
        File   sites_vcf_gz                = "${sample_name}.sites.vcf.gz"
        Int    assembly_length             = read_int("assembly_length")
        Int    assembly_length_unambiguous = read_int("assembly_length_unambiguous")
        Int    dist_to_ref_snps            = read_int("num_snps")
        Int    dist_to_ref_indels          = read_int("num_indels")
        String viralngs_version            = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 15]) + " GB"
        cpu: 8
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x8"
        maxRetries: 2
    }
}


task refine_2x_and_plot {
    meta {
      description: "This combined task exists just to streamline the two calls to assembly.refine and one call to reports.plot_coverage that many denovo assembly workflows use. It saves on instance spin up and docker pull times, file staging time, and all steps contained here have similar hardware requirements. The more atomic WDL tasks are still available for custom workflows (see refine, refine_assembly_with_aligned_reads, align_reads, etc)."
    }

    input {
      File    assembly_fasta
      File    reads_unmapped_bam

      File?   novocraft_license

      String? refine1_novoalign_options = "-r Random -l 30 -g 40 -x 20 -t 502"
      Float?  refine1_major_cutoff = 0.5
      Int?    refine1_min_coverage = 2

      String? refine2_novoalign_options = "-r Random -l 40 -g 40 -x 20 -t 100"
      Float?  refine2_major_cutoff = 0.5
      Int?    refine2_min_coverage = 3

      String? plot_coverage_novoalign_options = "-r Random -l 40 -g 40 -x 20 -t 100 -k"

      Int?    machine_mem_gb
      String  docker = "quay.io/broadinstitute/viral-assemble:2.1.33.0"

      # do this in two steps in case the input doesn't actually have "cleaned" in the name
      String  sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")
    }

    Int disk_size = 375

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        assembly.py --version | tee VERSION

        ln -s ~{assembly_fasta} assembly.fasta
        read_utils.py novoindex \
        assembly.fasta \
        ~{"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

        # refine 1
        assembly.py refine_assembly \
          assembly.fasta \
          ~{reads_unmapped_bam} \
          ~{sample_name}.refine1.fasta \
          --outVcf ~{sample_name}.refine1.pre_fasta.vcf.gz \
          --min_coverage ~{refine1_min_coverage} \
          --major_cutoff ~{refine1_major_cutoff} \
          --novo_params="~{refine1_novoalign_options}" \
          --JVMmemory "$mem_in_mb"m \
          ~{"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG

        # refine 2
        assembly.py refine_assembly \
          ~{sample_name}.refine1.fasta \
          ~{reads_unmapped_bam} \
          ~{sample_name}.fasta \
          --outVcf ~{sample_name}.refine2.pre_fasta.vcf.gz \
          --min_coverage ~{refine2_min_coverage} \
          --major_cutoff ~{refine2_major_cutoff} \
          --novo_params="~{refine2_novoalign_options}" \
          ~{"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --JVMmemory "$mem_in_mb"m \
          --loglevel=DEBUG

        # final alignment
        read_utils.py align_and_fix \
          ~{reads_unmapped_bam} \
          ~{sample_name}.fasta \
          --outBamAll ~{sample_name}.all.bam \
          --outBamFiltered ~{sample_name}.mapped.bam \
          --aligner_options "~{plot_coverage_novoalign_options}" \
          --JVMmemory "$mem_in_mb"m \
          ~{"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG

        # collect figures of merit
        set +o pipefail # grep will exit 1 if it fails to find the pattern
        grep -v '^>' ~{sample_name}.fasta | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' ~{sample_name}.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
        samtools view -c ~{sample_name}.mapped.bam | tee reads_aligned
        # report only primary alignments 260=exclude unaligned reads and secondary mappings
        samtools view -h -F 260 ~{sample_name}.all.bam | samtools flagstat - | tee ~{sample_name}.all.bam.flagstat.txt
        grep properly ~{sample_name}.all.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
        samtools view ~{sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
        #echo $(( $(cat bases_aligned) / $(cat assembly_length) )) | tee mean_coverage
        python -c "print (float("$(cat bases_aligned)")/"$(cat assembly_length)") if "$(cat assembly_length)">0 else print(0)" > mean_coverage

        # fastqc mapped bam
        reports.py fastqc ~{sample_name}.mapped.bam ~{sample_name}.mapped_fastqc.html --out_zip ~{sample_name}.mapped_fastqc.zip

        # plot coverage
        if [ $(cat reads_aligned) != 0 ]; then
          reports.py plot_coverage \
            ~{sample_name}.mapped.bam \
            ~{sample_name}.coverage_plot.pdf \
            --outSummary "~{sample_name}.coverage_plot.txt" \
            --plotFormat pdf \
            --plotWidth 1100 \
            --plotHeight 850 \
            --plotDPI 100 \
            --plotTitle "~{sample_name} coverage plot" \
            --loglevel=DEBUG
        else
          touch ~{sample_name}.coverage_plot.pdf ~{sample_name}.coverage_plot.txt
        fi
    }

    output {
        File   refine1_sites_vcf_gz          = "~{sample_name}.refine1.pre_fasta.vcf.gz"
        File   refine1_assembly_fasta        = "~{sample_name}.refine1.fasta"
        File   refine2_sites_vcf_gz          = "~{sample_name}.refine2.pre_fasta.vcf.gz"
        File   final_assembly_fasta          = "~{sample_name}.fasta"
        File   aligned_bam                   = "~{sample_name}.all.bam"
        File   aligned_bam_idx               = "~{sample_name}.all.bai"
        File   aligned_bam_flagstat          = "~{sample_name}.all.bam.flagstat.txt"
        File   aligned_only_reads_bam        = "~{sample_name}.mapped.bam"
        File   aligned_only_reads_bam_idx    = "~{sample_name}.mapped.bai"
        File   aligned_only_reads_fastqc     = "~{sample_name}.mapped_fastqc.html"
        File   aligned_only_reads_fastqc_zip = "~{sample_name}.mapped_fastqc.zip"
        File   coverage_plot                 = "~{sample_name}.coverage_plot.pdf"
        File   coverage_tsv                  = "~{sample_name}.coverage_plot.txt"
        Int    assembly_length               = read_int("assembly_length")
        Int    assembly_length_unambiguous   = read_int("assembly_length_unambiguous")
        Int    reads_aligned                 = read_int("reads_aligned")
        Int    read_pairs_aligned            = read_int("read_pairs_aligned")
        Float  bases_aligned                 = read_float("bases_aligned")
        Float  mean_coverage                 = read_float("mean_coverage")
        String viralngs_version              = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 8
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x8"
        maxRetries: 2
    }
}

task run_discordance {
    meta {
      description: "This step evaluates discordance between sequencing runs of the same sample. The input is a merged, aligned BAM file for a single sample. If multiple runs (read groups) exist, we split the aligned reads by read group and separately evaluate consensus calls per read group using bcftools mpileup and call. A VCF is emitted that describes variation between runs."
    }

    input {
      File   reads_aligned_bam
      File   reference_fasta
      String out_basename = "run"
      Int    min_coverage = 4

      String docker = "quay.io/broadinstitute/viral-core:2.1.33"
    }

    Int disk_size = 100

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION

        # create 2-col table with read group ids in both cols
        python3 <<CODE
        import tools.samtools
        header = tools.samtools.SamtoolsTool().getHeader("~{reads_aligned_bam}")
        rgids = [[x[3:] for x in h if x.startswith('ID:')][0] for h in header if h[0]=='@RG']
        n_rgs = len(rgids)
        with open('readgroups.txt', 'wt') as outf:
          for rg in rgids:
            outf.write(rg+'\t'+rg+'\n')
        n_lbs = len(set([[x[3:] for x in h if x.startswith('LB:')][0] for h in header if h[0]=='@RG']))
        with open('num_read_groups', 'wt') as outf:
          outf.write(str(n_rgs)+'\n')
        with open('num_libraries', 'wt') as outf:
          outf.write(str(n_lbs)+'\n')
        CODE

        # bcftools call snps while treating each RG as a separate sample
        bcftools mpileup \
          -G readgroups.txt -d 10000 -a "FORMAT/DP,FORMAT/AD" \
          -q 1 -m 2 -Ou \
          -f "~{reference_fasta}" "~{reads_aligned_bam}" \
          | bcftools call \
          -P 0 -m --ploidy 1 \
          --threads $(nproc) \
          -Ov -o everything.vcf

        # mask all GT calls when less than 3 reads
        cat everything.vcf | bcftools filter -e "FMT/DP<~{min_coverage}" -S . > filtered.vcf
        cat filtered.vcf | bcftools filter -i "MAC>0" > "~{out_basename}.discordant.vcf"

        # tally outputs
        bcftools filter -i 'MAC=0' filtered.vcf | bcftools query -f '%POS\n' | wc -l | tee num_concordant
        bcftools filter -i 'TYPE="snp"'  "~{out_basename}.discordant.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_discordant_snps
        bcftools filter -i 'TYPE!="snp"' "~{out_basename}.discordant.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_discordant_indels
    }

    output {
        File   discordant_sites_vcf = "~{out_basename}.discordant.vcf"
        Int    concordant_sites     = read_int("num_concordant")
        Int    discordant_snps      = read_int("num_discordant_snps")
        Int    discordant_indels    = read_int("num_discordant_indels")
        Int    num_read_groups      = read_int("num_read_groups")
        Int    num_libraries        = read_int("num_libraries")
        String viralngs_version     = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: "3 GB"
        cpu: 2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
        maxRetries: 2
    }
}


task filter_bad_ntc_batches {
    meta {
        description: "Identify NTCs (negative control libraries) that assemble too much of a genome and fail all other genomes from the same library prep batch."
    }
    input {
        File        seqid_list
        File        demux_meta_by_sample_json
        File        assembly_meta_tsv
        Int         ntc_min_unambig
        File?       genome_status_json
    }
    Int disk_size = 50
    command <<<
        set -e
        python3<<CODE
        import csv
        import json

        # load inputs
        ntc_min_unambig = ~{ntc_min_unambig}
        with open('~{seqid_list}', 'rt') as inf:
          seqid_list = list(x.strip() for x in inf)
        num_provided = len(seqid_list)
        with open('~{demux_meta_by_sample_json}', 'rt') as inf:
          demux_meta_orig = json.load(inf)
        with open('~{assembly_meta_tsv}', 'rt') as inf:
          assembly_meta = list(csv.DictReader(inf, delimiter='\t'))
        genome_status_json = '~{default="" genome_status_json}'
        if genome_status_json:
          with open(genome_status_json, 'rt') as inf:
            fail_meta = json.load(inf)
        else:
          fail_meta = {}

        # re-index demux_meta lookup table by sample_original instead of sample_sanitized
        demux_meta = dict((v['sample_original'],v) for k,v in demux_meta_orig.items())

        # identify bad NTCs
        reject_lanes = set()
        reject_batches = set()
        for sample in assembly_meta:
          if (demux_meta[sample['sample']].get('control') == 'NTC'):
            bad_ntc = sample['assembly_length_unambiguous'] \
              and (int(sample['assembly_length_unambiguous']) >= ntc_min_unambig)
            id = sample['sample']
            lane = demux_meta[sample['sample']]['run'].split('.')[-1]
            batch = demux_meta[sample['sample']].get('batch_lib','')
            print(f"NTC:\t{id}\t{sample['assembly_length_unambiguous']}\t{bad_ntc}\t{lane}\t{batch}")
            if bad_ntc:
              if batch:
                reject_batches.add(batch)
              else:
                reject_lanes.add(lane)
        print(f"BAD BATCHES:\t{','.join(reject_batches)}")
        print(f"BAD LANES:\t{','.join(reject_lanes)}")

        # filter samples from bad batches/lanes
        fail_meta = {}
        fails = set()
        for sample in assembly_meta:
          id = sample['sample']
          if (demux_meta[id].get('batch_lib') in reject_batches) \
            or (demux_meta[id]['run'].split('.')[-1] in reject_lanes):
            fail_meta[id] = 'failed_NTC'
            fails.add(id)
        seqid_list = list(id for id in seqid_list if id not in fails)

        # write outputs
        with open('seqids.filtered.txt', 'wt') as outf:
          for id in seqid_list:
            outf.write(id+'\n')
        with open('NUM_PROVIDED', 'wt') as outf:
          outf.write(str(num_provided)+'\n')
        with open('NUM_KEPT', 'wt') as outf:
          outf.write(str(len(seqid_list))+'\n')
        with open('REJECT_BATCHES', 'wt') as outf:
          for x in sorted(reject_batches):
            outf.write(x+'\n')
        with open('REJECT_LANES', 'wt') as outf:
          for x in sorted(reject_lanes):
            outf.write(x+'\n')
        with open('genome_status.json', 'wt') as outf:
          json.dump(fail_meta, outf, indent=2)

        CODE
    >>>
    runtime {
        docker: "python:slim"
        memory: "2 GB"
        cpu:    1
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File           seqids_kept      = "seqids.filtered.txt"
        Int            num_provided     = read_int("NUM_PROVIDED")
        Int            num_kept         = read_int("NUM_KEPT")
        Array[String]  reject_batches   = read_lines("REJECT_BATCHES")
        Array[String]  reject_lanes     = read_lines("REJECT_LANES")
        File           fail_meta_json   = "genome_status.json"
    }
}
