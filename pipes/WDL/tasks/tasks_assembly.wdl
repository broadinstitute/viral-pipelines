version 1.0

task assemble {
    input {
      File     reads_unmapped_bam
      File     trim_clip_db

      Int?     trinity_n_reads=250000
      Int?     spades_n_reads=10000000
      Int?     spades_min_contig_len=0

      String?  assembler="trinity"  # trinity, spades, or trinity-spades
      Boolean? always_succeed=false

      # do this in two steps in case the input doesn't actually have "taxfilt" in the name
      String   sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt")

      Int?     machine_mem_gb
      String   docker="quay.io/broadinstitute/viral-assemble:2.1.4.0"
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
        mem_in_gb=$(/opt/viral-ngs/source/docker/calc_mem.py gb 90)

        assembly.py --version | tee VERSION

        if [[ "${assembler}" == "trinity" ]]; then
          assembly.py assemble_trinity \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            ${'--n_reads=' + trinity_n_reads} \
            ${true='--alwaysSucceed' false="" always_succeed} \
            --JVMmemory "$mem_in_mb"m \
            --outReads=${sample_name}.subsamp.bam \
            --loglevel=DEBUG

        elif [[ "${assembler}" == "spades" ]]; then
          assembly.py assemble_spades \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            ${'--nReads=' + spades_n_reads} \
            ${true="--alwaysSucceed" false="" always_succeed} \
            ${'--minContigLen=' + spades_min_contig_len} \
            --memLimitGb $mem_in_gb \
            --outReads=${sample_name}.subsamp.bam \
            --loglevel=DEBUG

        elif [[ "${assembler}" == "trinity-spades" ]]; then
          assembly.py assemble_trinity \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-trinity.fasta \
            ${'--n_reads=' + trinity_n_reads} \
            --JVMmemory "$mem_in_mb"m \
            --outReads=${sample_name}.subsamp.bam \
            ${true='--always_succeed' false='' always_succeed} \
            --loglevel=DEBUG
          assembly.py assemble_spades \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            --contigsUntrusted=${sample_name}.assembly1-trinity.fasta \
            ${'--nReads=' + spades_n_reads} \
            ${true='--alwaysSucceed' false='' always_succeed} \
            ${'--minContigLen=' + spades_min_contig_len} \
            --memLimitGb $mem_in_gb \
            --loglevel=DEBUG

        else
          echo "unrecognized assembler ${assembler}" >&2
          exit 1
        fi

        samtools view -c ${sample_name}.subsamp.bam | tee subsample_read_count >&2
    }

    output {
        File   contigs_fasta        = "${sample_name}.assembly1-${assembler}.fasta"
        File   subsampBam           = "${sample_name}.subsamp.bam"
        Int    subsample_read_count = read_int("subsample_read_count")
        String viralngs_version     = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 15]) + " GB"
        cpu: 4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

}

task scaffold {
    input {
      File         contigs_fasta
      File         reads_bam
      Array[File]+ reference_genome_fasta

      String?      aligner
      Float?       min_length_fraction
      Float?       min_unambig
      Int?         replace_length=55

      Int?         nucmer_max_gap
      Int?         nucmer_min_match
      Int?         nucmer_min_cluster
      Float?       scaffold_min_pct_contig_aligned

      Int?         machine_mem_gb
      String       docker="quay.io/broadinstitute/viral-assemble:2.1.4.0"

      # do this in multiple steps in case the input doesn't actually have "assembly1-x" in the name
      String       sample_name = basename(basename(basename(contigs_fasta, ".fasta"), ".assembly1-trinity"), ".assembly1-spades")
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_gb=$(/opt/viral-ngs/source/docker/calc_mem.py gb 90)

        assembly.py --version | tee VERSION

        assembly.py order_and_orient \
          ${contigs_fasta} \
          ${sep=' ' reference_genome_fasta} \
          ${sample_name}.intermediate_scaffold.fasta \
          ${'--maxgap=' + nucmer_max_gap} \
          ${'--minmatch=' + nucmer_min_match} \
          ${'--mincluster=' + nucmer_min_cluster} \
          ${'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned} \
          --outReference ${sample_name}.scaffolding_chosen_ref.fasta \
          --outStats ${sample_name}.scaffolding_stats.txt \
          --outAlternateContigs ${sample_name}.scaffolding_alt_contigs.fasta \
          --loglevel=DEBUG

        grep '^>' ${sample_name}.scaffolding_chosen_ref.fasta | cut -c 2- | tr '\n' '\t' > ${sample_name}.scaffolding_chosen_ref.txt

        assembly.py gapfill_gap2seq \
          ${sample_name}.intermediate_scaffold.fasta \
          ${reads_bam} \
          ${sample_name}.intermediate_gapfill.fasta \
          --memLimitGb $mem_in_gb \
          --maskErrors \
          --loglevel=DEBUG

        grep -v '^>' ${sample_name}.intermediate_gapfill.fasta | tr -d '\n' | wc -c | tee assembly_preimpute_length
        grep -v '^>' ${sample_name}.intermediate_gapfill.fasta | tr -d '\nNn' | wc -c | tee assembly_preimpute_length_unambiguous

        assembly.py impute_from_reference \
          ${sample_name}.intermediate_gapfill.fasta \
          ${sample_name}.scaffolding_chosen_ref.fasta \
          ${sample_name}.scaffolded_imputed.fasta \
          --newName ${sample_name} \
          ${'--replaceLength=' + replace_length} \
          ${'--minLengthFraction=' + min_length_fraction} \
          ${'--minUnambig=' + min_unambig} \
          ${'--aligner=' + aligner} \
          --loglevel=DEBUG
    }

    output {
        File   scaffold_fasta                        = "${sample_name}.scaffolded_imputed.fasta"
        File   intermediate_scaffold_fasta           = "${sample_name}.intermediate_scaffold.fasta"
        File   intermediate_gapfill_fasta            = "${sample_name}.intermediate_gapfill.fasta"
        Int    assembly_preimpute_length             = read_int("assembly_preimpute_length")
        Int    assembly_preimpute_length_unambiguous = read_int("assembly_preimpute_length_unambiguous")
        String scaffolding_chosen_ref_name           = read_string("${sample_name}.scaffolding_chosen_ref.txt")
        File   scaffolding_chosen_ref                = "${sample_name}.scaffolding_chosen_ref.fasta"
        File   scaffolding_stats                     = "${sample_name}.scaffolding_stats.txt"
        File   scaffolding_alt_contigs               = "${sample_name}.scaffolding_alt_contigs.fasta"
        String viralngs_version                      = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 15]) + " GB"
        cpu: 4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task ivar_trim {
    meta {
      description: "this runs ivar trim on aligned reads, which results in soft-clipping of alignments"
    }

    input {
      File    aligned_bam
      File?   trim_coords_bed
      Int?    min_keep_length
      Int?    sliding_window
      Int?    min_quality=1

      Int?    machine_mem_gb
      String  docker="andersenlabapps/ivar:1.2.2"
    }

    String  bam_basename=basename(aligned_bam, ".bam")

    parameter_meta {
      aligned_bam:     { description: "aligned reads in BAM format", patterns: ["*.bam"] }
      trim_coords_bed: { description: "optional primers to trim in reference coordinate space (0-based BED format)", patterns: ["*.bed"] }
      min_keep_length: { description: "Minimum length of read to retain after trimming (Default: 30)" }
      sliding_window:  { description: "Width of sliding window for quality trimming (Default: 4)" }
      min_quality:     { description: "Minimum quality threshold for sliding window to pass (Default: 20)" }
    }

    command {
        ivar version | head -1 | tee VERSION
        if [ -f "${trim_coords_bed}" ]; then
          ivar trim -e \
            ${'-b ' + trim_coords_bed} \
            ${'-m ' + min_keep_length} \
            ${'-s ' + sliding_window} \
            ${'-q ' + min_quality} \
            -i ${aligned_bam} -p trim
          samtools sort -@ $(nproc) -m 1000M -o ${bam_basename}.trimmed.bam trim.bam
        else
          echo "skipping ivar trim"
          cp "${aligned_bam}" "${bam_basename}.trimmed.bam"
        fi
    }

    output {
        File   aligned_trimmed_bam = "${bam_basename}.trimmed.bam"
        String ivar_version        = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x4"
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

    String   aligner="minimap2"
    String?  aligner_options
    Boolean? skip_mark_dupes=false

    String   docker="quay.io/broadinstitute/viral-core:2.1.8"

    String   sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")
  }

  parameter_meta {
    aligner: { description: "Short read aligner to use: novoalign, minimap2, or bwa. (Default: novoalign)" }
  }
  
  command {
    set -ex # do not set pipefail, since grep exits 1 if it can't find the pattern

    read_utils.py --version | tee VERSION

    cp ${reference_fasta} assembly.fasta
    grep -v '^>' assembly.fasta | tr -d '\n' | wc -c | tee assembly_length

    if [ "$(cat assembly_length)" != "0" ]; then

      # only perform the following if the reference is non-empty

      if [ "${aligner}" == "novoalign" ]; then
        read_utils.py novoindex \
          assembly.fasta \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG
      fi
      read_utils.py index_fasta_picard assembly.fasta --loglevel=DEBUG
      read_utils.py index_fasta_samtools assembly.fasta --loglevel=DEBUG

      read_utils.py align_and_fix \
        ${reads_unmapped_bam} \
        assembly.fasta \
        --outBamAll "${sample_name}.all.bam" \
        --outBamFiltered "${sample_name}.mapped.bam" \
        --aligner ${aligner} \
        ${'--aligner_options "' + aligner_options + '"'} \
        ${true='--skipMarkDupes' false="" skip_mark_dupes} \
        --JVMmemory=3g \
        ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

    else
      # handle special case of empty reference fasta -- emit empty bams (with original bam headers)
      samtools view -H -b "${reads_unmapped_bam}" > "${sample_name}.all.bam"
      samtools view -H -b "${reads_unmapped_bam}" > "${sample_name}.mapped.bam"

      samtools index "${sample_name}.all.bam" "${sample_name}.all.bai"
      samtools index "${sample_name}.mapped.bam" "${sample_name}.mapped.bai"
    fi

    cat /proc/loadavg > CPU_LOAD

    # collect figures of merit
    grep -v '^>' assembly.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    samtools view -c ${reads_unmapped_bam} | tee reads_provided
    samtools view -c ${sample_name}.mapped.bam | tee reads_aligned
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${sample_name}.all.bam | samtools flagstat - | tee ${sample_name}.all.bam.flagstat.txt
    grep properly ${sample_name}.all.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("$(cat bases_aligned)")/"$(cat assembly_length_unambiguous)") if "$(cat assembly_length_unambiguous)">0 else print(0)" > mean_coverage

    # fastqc mapped bam
    reports.py fastqc ${sample_name}.mapped.bam ${sample_name}.mapped_fastqc.html --out_zip ${sample_name}.mapped_fastqc.zip

    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
  }

  output {
    File   aligned_bam                   = "${sample_name}.all.bam"
    File   aligned_bam_idx               = "${sample_name}.all.bai"
    File   aligned_bam_flagstat          = "${sample_name}.all.bam.flagstat.txt"
    File   aligned_only_reads_bam        = "${sample_name}.mapped.bam"
    File   aligned_only_reads_bam_idx    = "${sample_name}.mapped.bai"
    File   aligned_only_reads_fastqc     = "${sample_name}.mapped_fastqc.html"
    File   aligned_only_reads_fastqc_zip = "${sample_name}.mapped_fastqc.zip"
    Int    reads_provided                = read_int("reads_provided")
    Int    reads_aligned                 = read_int("reads_aligned")
    Int    read_pairs_aligned            = read_int("read_pairs_aligned")
    Float  bases_aligned                 = read_float("bases_aligned")
    Float  mean_coverage                 = read_float("mean_coverage")
    Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
    Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
    String cpu_load = read_string("CPU_LOAD")
    String viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu: 8
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x8"
    preemptible: 1
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

      Boolean? mark_duplicates=false
      Float?   major_cutoff=0.5
      Int?     min_coverage=3

      Int?     machine_mem_gb
      String   docker="quay.io/broadinstitute/viral-assemble:2.1.4.0"
    }

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

        if [ ${true='true' false='false' mark_duplicates} == "true" ]; then
          read_utils.py mkdup_picard \
            ${reads_aligned_bam} \
            temp_markdup.bam \
            --JVMmemory "$mem_in_mb"m \
            --loglevel=DEBUG
        else
          ln -s ${reads_aligned_bam} temp_markdup.bam
        fi
        samtools index -@ $(nproc) temp_markdup.bam temp_markdup.bai

        ln -s ${reference_fasta} assembly.fasta
        assembly.py refine_assembly \
          assembly.fasta \
          temp_markdup.bam \
          refined.fasta \
          --already_realigned_bam=temp_markdup.bam \
          --outVcf ${sample_name}.sites.vcf.gz \
          --min_coverage ${min_coverage} \
          --major_cutoff ${major_cutoff} \
          --JVMmemory "$mem_in_mb"m \
          --loglevel=DEBUG

        file_utils.py rename_fasta_sequences \
          refined.fasta "${sample_name}.fasta" "${sample_name}"

        # collect variant counts
        bcftools filter -e "FMT/DP<${min_coverage}" -S . "${sample_name}.sites.vcf.gz" -Ou | bcftools filter -i "AC>1" -Ou > "${sample_name}.diffs.vcf"
        bcftools filter -i 'TYPE="snp"'  "${sample_name}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_snps
        bcftools filter -i 'TYPE!="snp"' "${sample_name}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_indels

        # collect figures of merit
        set +o pipefail # grep will exit 1 if it fails to find the pattern
        grep -v '^>' refined.fasta | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' refined.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    }

    output {
        File   refined_assembly_fasta       = "${sample_name}.fasta"
        File   sites_vcf_gz                 = "${sample_name}.sites.vcf.gz"
        Int    assembly_length              = read_int("assembly_length")
        Int    assembly_length_unambiguous  = read_int("assembly_length_unambiguous")
        Int    dist_to_ref_snps             = read_int("num_snps")
        Int    dist_to_ref_indels           = read_int("num_indels")
        String viralngs_version             = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 8
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task refine {
    meta {
      description: "This step refines/polishes a genome based on unmapped short reads, aligning them (with either novoalign or bwa), and producing a new consensus genome. Uses GATK3 Unified Genotyper to produce new consensus. Produces new genome (fasta), variant calls (VCF), and figures of merit."
    }

    input {
      File    assembly_fasta
      File    reads_unmapped_bam

      File?   novocraft_license

      String? novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100"
      Float?  major_cutoff=0.5
      Int?    min_coverage=1

      Int?    machine_mem_gb
      String  docker="quay.io/broadinstitute/viral-assemble:2.1.4.0"

      String  assembly_basename=basename(basename(assembly_fasta, ".fasta"), ".scaffold")
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        assembly.py --version | tee VERSION

        ln -s ${assembly_fasta} assembly.fasta
        read_utils.py novoindex \
        assembly.fasta \
        ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

        assembly.py refine_assembly \
          assembly.fasta \
          ${reads_unmapped_bam} \
          ${assembly_basename}.refined.fasta \
          --outVcf ${assembly_basename}.sites.vcf.gz \
          --min_coverage ${min_coverage} \
          --major_cutoff ${major_cutoff} \
          --novo_params="${novoalign_options}" \
          --JVMmemory "$mem_in_mb"m \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG
    }

    output {
        File   refined_assembly_fasta  = "${assembly_basename}.refined.fasta"
        File   sites_vcf_gz            = "${assembly_basename}.sites.vcf.gz"
        String viralngs_version        = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 8
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
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

      String? refine1_novoalign_options="-r Random -l 30 -g 40 -x 20 -t 502"
      Float?  refine1_major_cutoff=0.5
      Int?    refine1_min_coverage=2

      String? refine2_novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100"
      Float?  refine2_major_cutoff=0.5
      Int?    refine2_min_coverage=3

      String? plot_coverage_novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100 -k"

      Int?    machine_mem_gb
      String  docker="quay.io/broadinstitute/viral-assemble:2.1.4.0"

      # do this in two steps in case the input doesn't actually have "cleaned" in the name
      String  sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        assembly.py --version | tee VERSION

        ln -s ${assembly_fasta} assembly.fasta
        read_utils.py novoindex \
        assembly.fasta \
        ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

        # refine 1
        assembly.py refine_assembly \
          assembly.fasta \
          ${reads_unmapped_bam} \
          ${sample_name}.refine1.fasta \
          --outVcf ${sample_name}.refine1.pre_fasta.vcf.gz \
          --min_coverage ${refine1_min_coverage} \
          --major_cutoff ${refine1_major_cutoff} \
          --novo_params="${refine1_novoalign_options}" \
          --JVMmemory "$mem_in_mb"m \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG

        # refine 2
        assembly.py refine_assembly \
          ${sample_name}.refine1.fasta \
          ${reads_unmapped_bam} \
          ${sample_name}.fasta \
          --outVcf ${sample_name}.refine2.pre_fasta.vcf.gz \
          --min_coverage ${refine2_min_coverage} \
          --major_cutoff ${refine2_major_cutoff} \
          --novo_params="${refine2_novoalign_options}" \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --JVMmemory "$mem_in_mb"m \
          --loglevel=DEBUG

        # final alignment
        read_utils.py align_and_fix \
          ${reads_unmapped_bam} \
          ${sample_name}.fasta \
          --outBamAll ${sample_name}.all.bam \
          --outBamFiltered ${sample_name}.mapped.bam \
          --aligner_options "${plot_coverage_novoalign_options}" \
          --JVMmemory "$mem_in_mb"m \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG

        # collect figures of merit
        set +o pipefail # grep will exit 1 if it fails to find the pattern
        grep -v '^>' ${sample_name}.fasta | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' ${sample_name}.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
        samtools view -c ${sample_name}.mapped.bam | tee reads_aligned
        # report only primary alignments 260=exclude unaligned reads and secondary mappings
        samtools view -h -F 260 ${sample_name}.all.bam | samtools flagstat - | tee ${sample_name}.all.bam.flagstat.txt
        grep properly ${sample_name}.all.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
        samtools view ${sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
        #echo $(( $(cat bases_aligned) / $(cat assembly_length) )) | tee mean_coverage
        python -c "print (float("$(cat bases_aligned)")/"$(cat assembly_length)") if "$(cat assembly_length)">0 else print(0)" > mean_coverage

        # fastqc mapped bam
        reports.py fastqc ${sample_name}.mapped.bam ${sample_name}.mapped_fastqc.html --out_zip ${sample_name}.mapped_fastqc.zip

        # plot coverage
        if [ $(cat reads_aligned) != 0 ]; then
          reports.py plot_coverage \
            ${sample_name}.mapped.bam \
            ${sample_name}.coverage_plot.pdf \
            --outSummary "${sample_name}.coverage_plot.txt" \
            --plotFormat pdf \
            --plotWidth 1100 \
            --plotHeight 850 \
            --plotDPI 100 \
            --plotTitle "${sample_name} coverage plot" \
            --loglevel=DEBUG
        else
          touch ${sample_name}.coverage_plot.pdf ${sample_name}.coverage_plot.txt
        fi
    }

    output {
        File refine1_sites_vcf_gz          = "${sample_name}.refine1.pre_fasta.vcf.gz"
        File refine1_assembly_fasta        = "${sample_name}.refine1.fasta"
        File refine2_sites_vcf_gz          = "${sample_name}.refine2.pre_fasta.vcf.gz"
        File final_assembly_fasta          = "${sample_name}.fasta"
        File aligned_bam                   = "${sample_name}.all.bam"
        File aligned_bam_idx               = "${sample_name}.all.bai"
        File aligned_bam_flagstat          = "${sample_name}.all.bam.flagstat.txt"
        File aligned_only_reads_bam        = "${sample_name}.mapped.bam"
        File aligned_only_reads_bam_idx    = "${sample_name}.mapped.bai"
        File aligned_only_reads_fastqc     = "${sample_name}.mapped_fastqc.html"
        File aligned_only_reads_fastqc_zip = "${sample_name}.mapped_fastqc.zip"
        File coverage_plot                 = "${sample_name}.coverage_plot.pdf"
        File coverage_tsv                  = "${sample_name}.coverage_plot.txt"
        Int  assembly_length               = read_int("assembly_length")
        Int  assembly_length_unambiguous   = read_int("assembly_length_unambiguous")
        Int  reads_aligned                 = read_int("reads_aligned")
        Int  read_pairs_aligned            = read_int("read_pairs_aligned")
        Float bases_aligned                 = read_float("bases_aligned")
        Float mean_coverage                = read_float("mean_coverage")
        String viralngs_version            = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 8
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task run_discordance {
    meta {
      description: "This step evaluates discordance between sequencing runs of the same sample. The input is a merged, aligned BAM file for a single sample. If multiple runs (read groups) exist, we split the aligned reads by read group and separately evaluate consensus calls per read group using bcftools mpileup and call. A VCF is emitted that describes variation between runs."
    }

    input {
      File     reads_aligned_bam
      File     reference_fasta
      String   out_basename = "run"
      Int      min_coverage=4

      String   docker="quay.io/broadinstitute/viral-core:2.1.8"
    }

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION

        # create 2-col table with read group ids in both cols
        python3 <<CODE
        import tools.samtools
        header = tools.samtools.SamtoolsTool().getHeader("${reads_aligned_bam}")
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
          -f "${reference_fasta}" "${reads_aligned_bam}" \
          | bcftools call \
          -P 0 -m --ploidy 1 \
          --threads $(nproc) \
          -Ov -o everything.vcf

        # mask all GT calls when less than 3 reads
        cat everything.vcf | bcftools filter -e "FMT/DP<${min_coverage}" -S . > filtered.vcf
        cat filtered.vcf | bcftools filter -i "MAC>0" > "${out_basename}.discordant.vcf"

        # tally outputs
        bcftools filter -i 'MAC=0' filtered.vcf | bcftools query -f '%POS\n' | wc -l | tee num_concordant
        bcftools filter -i 'TYPE="snp"'  "${out_basename}.discordant.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_discordant_snps
        bcftools filter -i 'TYPE!="snp"' "${out_basename}.discordant.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_discordant_indels
    }

    output {
        File   discordant_sites_vcf = "${out_basename}.discordant.vcf"
        Int    concordant_sites  = read_int("num_concordant")
        Int    discordant_snps   = read_int("num_discordant_snps")
        Int    discordant_indels = read_int("num_discordant_indels")
        Int    num_read_groups   = read_int("num_read_groups")
        Int    num_libraries     = read_int("num_libraries")
        String viralngs_version  = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: "3 GB"
        cpu: 2
        disks: "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
}
