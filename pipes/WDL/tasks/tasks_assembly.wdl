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
      String   docker = "quay.io/broadinstitute/viral-assemble:2.4.3.1"
    }
    parameter_meta{
      reads_unmapped_bam: {
        description: "Unaligned reads in BAM format.",
        patterns: ["*.bam"],
        category: "required"
      }
      trim_clip_db: {
        description: "Trimmomatic clip database.",
        category: "required"
      }
      spades_n_reads: {
        description: "Subsample reads threshold prior to assembly. Default set to 10000000",
        category: "required"
      }
      spades_min_contig_len: {
        description: "Minimum length of output contig.",
        category: "other"
      }
      spades_options: {
        description: "Display additional options to pass the SPAdes assembler.",
        category: "other"
      }
      always_succeed: {
        description: "In the event that assembly fails for any reason, output an empty contigs file instead of an error code.",
        cateogory: "other"
      }
      contigs_fasta: {
        description:"De novo RNA-seq assembly of contigs with the SPAdes assembler in FASTA format.",
        patterns: ["*.fasta"],
        category: "other"
      }
      subsampBam:{
        description: "Subsample your reads to speed up the process of the run.",
        category: "other"
      }
      subsample_read_count:{
        description: "Number of reads that your subsample contains.",
        category: "other"
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
        memory: select_first([machine_mem_gb, 63]) + " GB"
        cpu: 4
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x8"
        maxRetries: 2
    }

}

task select_references {
  meta {
    description: "Evaluate reference genomes based on ANI similarity to provided contigs. This will emit an ANI-rank-ordered table of references and will also cluster reference genomes based on ANI similarity to each other, picking only the top hit (based on ANI similairty to contigs) of each cluster. Default parameters for skani are tuned for viral genomes (-m 50 -s 85 --no-learned-ani --slow --robust --no-marker-index). This tool can tolerate a large number of reference genomes as input."
  }
  input {
    Array[File]   reference_genomes_fastas
    File          contigs_fasta

    Int?          skani_m
    Int?          skani_s
    Int?          skani_c
    Int?          skani_n

    String        docker = "quay.io/broadinstitute/viral-assemble:2.4.3.1"
    Int           machine_mem_gb = 4
    Int           cpu = 2
    Int           disk_size = 100
  }
  String contigs_basename = basename(basename(contigs_fasta, '.fasta'), '.assembly1-spades')

  command <<<
    set -e

    # run skani, find top hits, cluster references
    assembly.py skani_contigs_to_refs \
      "~{contigs_fasta}" \
      "~{sep='" "' reference_genomes_fastas}" \
      "~{contigs_basename}.refs_skani_dist.full.tsv" \
      "~{contigs_basename}.refs_skani_dist.top.tsv" \
      "~{contigs_basename}.ref_clusters.tsv" \
      ~{'-m ' + skani_m} \
      ~{'-s ' + skani_s} \
      ~{'-c ' + skani_c} \
      ~{'-n ' + skani_n} \
      --loglevel=DEBUG

    # create basename-only version of ref_clusters output file
    # create tar-bundles of ref_clusters fastas, since Cromwell doesn't delocalize files in a Array[Array[File]] = read_tsv
    python3 <<CODE
    import os, os.path, shutil, tarfile
    os.mkdir("clusters")
    with open("~{contigs_basename}.ref_clusters.tsv", 'r') as inf:
      with open("~{contigs_basename}.ref_clusters.basenames.tsv", 'w') as outf:
        for line in inf:
          fnames = line.strip().split('\t')
          assert fnames
          basefnames = list([f[:-6] if f.endswith('.fasta') else f for f in map(os.path.basename, fnames)])
          outf.write('\t'.join(basefnames) + '\n')

          with tarfile.open(os.path.join("clusters", basefnames[0] + "-" + str(len(basefnames))  + ".tar.gz"), "w:gz") as tarball:
            for f in fnames:
              shutil.copy(f, ".")
              tarball.add(os.path.basename(f))
              os.unlink(os.path.basename(f))
    CODE

    # create top-hits output files
    cut -f 1 "~{contigs_basename}.refs_skani_dist.top.tsv" | tail +2 > TOP_FASTAS
    for f in $(cat TOP_FASTAS); do basename "$f" .fasta; done > TOP_FASTAS_BASENAMES
  >>>

  output {
    Array[File]          matched_reference_clusters_fastas_tars = glob("clusters/*.tar.gz")
    Array[Array[String]] matched_reference_clusters_basenames = read_tsv("~{contigs_basename}.ref_clusters.basenames.tsv")
    Array[String]        top_matches_per_cluster_basenames = read_lines("TOP_FASTAS_BASENAMES")
    Array[File]          top_matches_per_cluster_fastas = read_lines("TOP_FASTAS")
    File                 skani_dist_full_tsv = "~{contigs_basename}.refs_skani_dist.full.tsv"
    File                 skani_dist_top_tsv  = "~{contigs_basename}.refs_skani_dist.top.tsv"
  }

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu:    cpu
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TESs
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 2
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
      Boolean      allow_incomplete_output = false

      Int?          skani_m
      Int?          skani_s
      Int?          skani_c

      Int?         nucmer_max_gap
      Int?         nucmer_min_match
      Int?         nucmer_min_cluster
      Int?         scaffold_min_contig_len
      Float?       scaffold_min_pct_contig_aligned

      Int?         machine_mem_gb
      String       docker="quay.io/broadinstitute/viral-assemble:2.4.3.1"

      # do this in multiple steps in case the input doesn't actually have "assembly1-x" in the name
      String       sample_name = basename(basename(contigs_fasta, ".fasta"), ".assembly1-spades")
    }
    parameter_meta {
      reads_bam: {
        description: "Reads in BAM format.",
        patterns: ["*.bam"],
        category: "required"
      }

      contigs_fasta: {
        description: "De novo contigs in fasta format",
        patterns: ["*.fasta"],
        category: "required"
      }

      reference_genome_fasta: {
        description: "Reference genomes to scaffold against. Multiple reference genomes may be provided, one per fasta file, and this task will attempt each of them and select the reference that produces the most complete output. Each reference genome should be in a single fasta file with all segments/chromosomes contained in the file as separate sequences in the correct order. Output genomes from this task will contain the same number of sequences as the selected input reference genome, and will be emitted in the same order.",
        patterns: ["*.fasta"],
        category: "required"
      }
      aligner: {
        description: "Alignment tools used to align the reference sequence to aligned contigs. Possible options: muscle, mafft, mummer (= nucmer), set to muscle for default.",
        cateogory: "advanced"
      }
      min_length_fraction: {
        description: "This step will fail with a PoorAssemblyError if the total end-to-end genome length in the output genome (inclusive of interior Ns) is less than this fraction of the length of the reference genome selected. Valid values are fractions from 0 to 1, default value is 0.5.",
        category: "common"
      }
      min_unambig: {
        description: "This step will fail with a PoorAssemblyError if the total number of unambiguous bases in the output genome (exclusive of interior Ns) is less than this fraction of its end-to-end length (inclusive of interior Ns). Valid values are fractions from 0 to 1, default value is 0.5.",
        category: "common"
      }
      replace_length: {
        description: "The first and last replace_length base pairs of each segment in the output genome will be replaced with the equivalent sequences in the reference genome as a mechanism to handle common assembly errors in repetitive or inverted regions that are common to chromosome/segment ends. Valid values are any non-negative integer. Default is 55 bp.",
        category: "advanced"
      }
      nucmer_max_gap: {
        description: "When scaffolding contigs to the reference via nucmer, this specifies the -g parameter to nucmer (the maximum allowed gap between adjacent matches in a cluster). Our default is 200 (up from nucmer default of 90), mummer documentation suggests it is valid to increase up to 1000 to allow for more diversity.",
        category: "advanced"
      }
      nucmer_min_match: {
        description: "When scaffolding contigs to the reference via nucmer, this specifies the -l parameter to nucmer (the minimal size of a maximal exact match). Our default is 10 (down from nucmer default of 20) to allow for more divergence.",
        category: "advanced"
      }
      nucmer_min_cluster:{
        description: "When scaffolding contigs to the reference via nucmer, this specifies the -c parameter to nucmer (minimum cluster length). Our default is the nucmer default of 65 bp.",
        category: "advanced"
      }
      scaffold_min_contig_len: {
        description: "Any sequences in contigs_fasta that are shorter than this length will be ignored for scaffolding.",
        category: "advanced"
      }
      scaffold_min_pct_contig_aligned: {
        description: "Any contig alignments to the reference scaffold that account for less than this fraction of the contig's length will be rejected for scaffolding. Valid values are fractions from 0 to 1; the default value is 0.3.",
        category: "advanced"
      }
      scaffold_fasta: {
        description: "This is the output draft genome after scaffolding contigs to references and imputing missing sequence from those references. This resulting genome is a hybrid of sequences from the de novo assembly and imputed reference sequence, and *requires* polishing with reads to be considered a valid consensus sequence. This is the final output of this task that should be used for polishing.",
        patterns: ["*.fasta"],
        category: "other"
      }
      intermediate_scaffold_fasta: {
        description: "This is the output draft genome after scaffolding contigs to reference genomes but prior to imputation with reference sequence or gapfilling with reads. The only unambiguous bases are from the contigs_fasta file.",
        patterns: ["*.fasta"],
        category: "other"
      }
    }

    Int disk_size = 375

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_gb=$(/opt/viral-ngs/source/docker/calc_mem.py gb 90)

        assembly.py --version | tee VERSION

        # use skani to choose top hit
        assembly.py skani_contigs_to_refs \
          "~{contigs_fasta}" \
          "~{sep='" "' reference_genome_fasta}" \
          "~{sample_name}.refs_skani_dist.full.tsv" \
          "~{sample_name}.refs_skani_dist.top.tsv" \
          "~{sample_name}.ref_clusters.tsv" \
          ~{'-m ' + skani_m} \
          ~{'-s ' + skani_s} \
          ~{'-c ' + skani_c} \
          --loglevel=DEBUG

        # sometimes skani fails; if so, just fall-back to sending all refs downstream
        if [[ $(wc -l <"~{sample_name}.refs_skani_dist.full.tsv") -ge 2 ]]; then
          # skani reference selection worked: just try one reference
          CHOSEN_REF_FASTA=$(cut -f 1 "~{sample_name}.refs_skani_dist.full.tsv" | tail +2 | head -1)
          basename "$CHOSEN_REF_FASTA" .fasta > CHOSEN_REF_BASENAME

          assembly.py order_and_orient \
            "~{contigs_fasta}" \
            "$CHOSEN_REF_FASTA" \
            "~{sample_name}".intermediate_scaffold.fasta \
            ~{'--min_contig_len=' + scaffold_min_contig_len} \
            ~{'--maxgap=' + nucmer_max_gap} \
            ~{'--minmatch=' + nucmer_min_match} \
            ~{'--mincluster=' + nucmer_min_cluster} \
            ~{'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned} \
            --outReference "~{sample_name}".scaffolding_chosen_ref.fasta \
            --outStats "~{sample_name}".scaffolding_stats.txt \
            --outAlternateContigs ~{sample_name}.scaffolding_alt_contigs.fasta \
            ~{true='--allow_incomplete_output' false="" allow_incomplete_output} \
            --loglevel=DEBUG

          cut -f 3 "~{sample_name}.refs_skani_dist.full.tsv" | tail +2 | head -1 > SKANI_ANI
          cut -f 4 "~{sample_name}.refs_skani_dist.full.tsv" | tail +2 | head -1 > SKANI_REF_AF
          cut -f 5 "~{sample_name}.refs_skani_dist.full.tsv" | tail +2 | head -1 > SKANI_CONTIGS_AF

        else
          # skani reference selection failed: try all references
          echo "0" > SKANI_ANI
          echo "0" > SKANI_REF_AF
          echo "0" > SKANI_CONTIGS_AF
          echo "" > CHOSEN_REF_BASENAME

          assembly.py order_and_orient \
            "~{contigs_fasta}" \
            "~{sep='" "' reference_genome_fasta}" \
            "~{sample_name}".intermediate_scaffold.fasta \
            ~{'--min_contig_len=' + scaffold_min_contig_len} \
            ~{'--maxgap=' + nucmer_max_gap} \
            ~{'--minmatch=' + nucmer_min_match} \
            ~{'--mincluster=' + nucmer_min_cluster} \
            ~{'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned} \
            --outReference "~{sample_name}".scaffolding_chosen_ref.fasta \
            --outStats "~{sample_name}".scaffolding_stats.txt \
            --outAlternateContigs ~{sample_name}.scaffolding_alt_contigs.fasta \
            ~{true='--allow_incomplete_output' false="" allow_incomplete_output} \
            --loglevel=DEBUG

        fi
        grep '^>' "~{sample_name}".scaffolding_chosen_ref.fasta | cut -c 2- | cut -f 1 -d ' ' > "~{sample_name}".scaffolding_chosen_refs.txt

        assembly.py gapfill_gap2seq \
          "~{sample_name}".intermediate_scaffold.fasta \
          "~{reads_bam}" \
          "~{sample_name}".intermediate_gapfill.fasta \
          --memLimitGb $mem_in_gb \
          --maskErrors \
          --loglevel=DEBUG

        set +e +o pipefail
        grep -v '^>' "~{sample_name}".intermediate_gapfill.fasta | tr -d '\n' | wc -c | tee assembly_preimpute_length
        grep -v '^>' "~{sample_name}".intermediate_gapfill.fasta | tr -d '\nNn' | wc -c | tee assembly_preimpute_length_unambiguous
        grep '^>' "~{sample_name}".intermediate_gapfill.fasta | wc -l | tee assembly_num_segments_recovered
        grep '^>' "~{sample_name}".scaffolding_chosen_ref.fasta | wc -l | tee reference_num_segments_required
        grep -v '^>' "~{sample_name}".scaffolding_chosen_ref.fasta | tr -d '\n' | wc -c | tee reference_length
        set -e -o pipefail

        if ~{true='true' false='false' allow_incomplete_output} && ! cmp -s assembly_num_segments_recovered reference_num_segments_required
        then
          # draft assembly does not have enough segments--and that's okay (allow_incomplete_output=true)
          file_utils.py rename_fasta_sequences \
            "~{sample_name}".intermediate_gapfill.fasta \
            "~{sample_name}".scaffolded_imputed.fasta \
            "~{sample_name}" --suffix_always --loglevel=DEBUG
        else
          # draft assembly must have the right number of segments (fail if not)
          assembly.py impute_from_reference \
            "~{sample_name}".intermediate_gapfill.fasta \
            "~{sample_name}".scaffolding_chosen_ref.fasta \
            "~{sample_name}".scaffolded_imputed.fasta \
            --newName "~{sample_name}" \
            ~{'--replaceLength=' + replace_length} \
            ~{'--minLengthFraction=' + min_length_fraction} \
            ~{'--minUnambig=' + min_unambig} \
            ~{'--aligner=' + aligner} \
            --loglevel=DEBUG
        fi
    }

    output {
        File   scaffold_fasta                        = "~{sample_name}.scaffolded_imputed.fasta"
        File   intermediate_scaffold_fasta           = "~{sample_name}.intermediate_scaffold.fasta"
        File   intermediate_gapfill_fasta            = "~{sample_name}.intermediate_gapfill.fasta"
        Int    assembly_preimpute_length             = read_int("assembly_preimpute_length")
        Int    assembly_preimpute_length_unambiguous = read_int("assembly_preimpute_length_unambiguous")
        Int    assembly_num_segments_recovered       = read_int("assembly_num_segments_recovered")
        Int    reference_num_segments_required       = read_int("reference_num_segments_required")
        Int    reference_length                      = read_int("reference_length")
        Array[String] scaffolding_chosen_ref_names   = read_lines("~{sample_name}.scaffolding_chosen_refs.txt")
        String scaffolding_chosen_ref_basename       = read_string("CHOSEN_REF_BASENAME")
        File   scaffolding_chosen_ref                = "~{sample_name}.scaffolding_chosen_ref.fasta"
        File   scaffolding_stats                     = "~{sample_name}.refs_skani_dist.full.tsv"
        File   scaffolding_alt_contigs               = "~{sample_name}.scaffolding_alt_contigs.fasta"
        Float  skani_ani                             = read_float("SKANI_ANI")
        Float  skani_ref_aligned_frac                = read_float("SKANI_REF_AF")
        Float  skani_contigs_aligned_frac            = read_float("SKANI_CONTIGS_AF")
        String viralngs_version                      = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 63]) + " GB"
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

      String  bam_basename=basename(aligned_bam, ".bam")
      Int disk_size = 375
    }
    parameter_meta {
      aligned_bam:{
        description: "aligned reads in BAM format",
        patterns: ["*.bam"],
        category: "required"
      }
      trim_coords_bed:{
        description: "optional primers to trim in reference coordinate space (0-based BED format)",
        patterns: ["*.bed"],
        category: "advanced"
      }
      min_keep_length:{
        description: "Minimum length of read to retain after trimming (Default: 30)",
        category: "advanced"
      }
      sliding_window:  {
        description: "Width of sliding window for quality trimming (Default: 4)",
        category: "advanced"
      }
      min_quality: {
        description: "Minimum quality threshold for sliding window to pass (Default: 20)",
        category: "advanced"
      }
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
    parameter_meta {
      ivar_trim_stats_tsv: {
        description: "Number of trimmed sequences based on a quality threshold set above.",
        category: "required"
      }
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
    Boolean  skip_mark_dupes = false

    Int?     machine_mem_gb
    String   docker = "quay.io/broadinstitute/viral-core:2.4.1"

    String   sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")
  }

  Int disk_size = 375

  parameter_meta {
    reference_fasta: {
      description: "Reference genome, in FASTA format, pre-indexed by Novoindex",
      category: "required"
    }
    reads_unmapped_bam: {
      description: "Unaligned reads in BAM format.",
      category: "required"
    }
    aligner: { 
      description: "Short read aligner to use novoalign, minimap2, or bwa. (Default novoalign)",
      category: "advanced"
      }
    skip_mark_dupes: {
      description: "If specific, duplicate reads will not be marked in the resulting output file.",
      category: "advanced"
    }
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
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    File   aligned_bam                   = "~{sample_name}.all.bam"
    File   aligned_bam_idx               = "~{sample_name}.all.bai"
    File   aligned_bam_flagstat          = "~{sample_name}.all.bam.flagstat.txt"
    File   aligned_only_reads_bam        = "~{sample_name}.mapped.bam"
    File   aligned_only_reads_bam_idx    = "~{sample_name}.mapped.bai"
    File   aligned_only_reads_fastqc     = "~{sample_name}.mapped_fastqc.html"
    File   aligned_only_reads_fastqc_zip = "~{sample_name}.mapped_fastqc.zip"
    Int    reference_length              = read_int("assembly_length")
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
      String   out_basename = basename(reads_aligned_bam, '.bam')
      String   sample_name = out_basename

      Boolean  mark_duplicates = false
      Float    major_cutoff = 0.5
      Int      min_coverage = 3

      Int      machine_mem_gb = 15
      String   docker = "quay.io/broadinstitute/viral-assemble:2.4.3.1"
    }

    Int disk_size = 375

    parameter_meta {
      reference_fasta:{
        description: "Reference genome, in FASTA format, pre-indexed by Novoindex",
        category: "required"
      }
      reads_aligned_bam: {
        description: "Aligned reads in BAM format.",
        patterns: ["*.bam"],
        category: "required"
      }  
      mark_duplicates:{
        description: "Instead of removing duplicates, simply marks them.",
        category: "advanced"
      }
      major_cutoff: {
        description: "If the major allele is present at a frequency higher than this cutoff, we will call an unambiguous base at that position.  If it is equal to or below this cutoff, we will call an ambiguous base representing all possible alleles at that position.",
        category: "advanced"
      }
      min_coverage: {
        description: "Minimum read coverage required to call a position unambiguous.",
        category: "advanaced"
      }
    }

    command <<<
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
          --outVcf "~{out_basename}.sites.vcf.gz" \
          --min_coverage ~{min_coverage} \
          --major_cutoff ~{major_cutoff} \
          --JVMmemory "$mem_in_mb"m \
          --loglevel=DEBUG

        # trim edges if they're too ambiguous or else table2asn will reject if first or last line of seq > 40% ambig
        # (we don't really know how long the last line will be so use both the --ten and --fifty rules
        # and basically disable the --maxfrac rule as we will handle that elsewhere)
        if [ -s refined.fasta ]; then
          fasta-trim-terminal-ambigs.pl --3rules \
            --ten 4 \
            --fifty 15 \
            --maxfrac 0.9 \
            refined.fasta > trimmed.fasta
        else
          # VADR perl script fails on empty input
          touch trimmed.fasta
        fi

        # rename for external consumption
        file_utils.py rename_fasta_sequences \
          trimmed.fasta "~{out_basename}.fasta" "~{sample_name}"

        # collect variant counts
        if (( $(cat trimmed.fasta | wc -l) > 1 )); then
          bcftools filter -e "FMT/DP<~{min_coverage}" -S . "~{out_basename}.sites.vcf.gz" -Ou | bcftools filter -i "AC>1" -Ou > "~{out_basename}.diffs.vcf"
          bcftools filter -i 'TYPE="snp"'  "~{out_basename}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_snps
          bcftools filter -i 'TYPE!="snp"' "~{out_basename}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_indels
        else
          # empty output
          echo "0" > num_snps
          echo "0" > num_indels
          cp "~{out_basename}.sites.vcf.gz" "~{out_basename}.diffs.vcf"
        fi

        # collect figures of merit
        set +o pipefail # grep will exit 1 if it fails to find the pattern
        grep -v '^>' trimmed.fasta | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' trimmed.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    >>>

    output {
        File   refined_assembly_fasta      = "~{out_basename}.fasta"
        File   sites_vcf_gz                = "~{out_basename}.sites.vcf.gz"
        Int    assembly_length             = read_int("assembly_length")
        Int    assembly_length_unambiguous = read_int("assembly_length_unambiguous")
        Int    dist_to_ref_snps            = read_int("num_snps")
        Int    dist_to_ref_indels          = read_int("num_indels")
        String viralngs_version            = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
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

      String docker = "quay.io/broadinstitute/viral-core:2.4.1"
    }
    parameter_meta {
      reads_aligned_bam: {
        description: "Unaligned reads in BAM Format",
        category: "required"
      }
      reference_fasta: {
        description: "Reference assembled genome in FASTA format ",
        category: "required"
      }
      discordant_sites_vcf:{
        description:"The SNPs between runs of the same sample. ",
        category: "other"
      }
    }
  

    Int disk_size = 100

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION

        python3 <<CODE
        # create 2-col table with read group ids in both cols
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

        # detect empty fasta situation and manually create empty VCF
        import os.path
        if (os.path.getsize('~{reference_fasta}') == 0):
          sample_name = [[x[3:] for x in h if x.startswith('SM:')][0] for h in header if h[0]=='@RG'][0]
          with open('everything.vcf', 'wt') as outf:
              outf.write('##fileformat=VCFv4.3')
              outf.write('##ALT=<ID=*,Description="Represents allele(s) other than observed.">')
              outf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">')
              outf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
              outf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">')
              outf.write('\t'.join(('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample_name))+'\n')
        CODE

        if [ ! -f everything.vcf ]; then
          # bcftools call snps while treating each RG as a separate sample

          # first index the ref fasta, iff it contains non-N sequence
          if [ $(grep -v '^>' "~{reference_fasta}" | tr -d '\nNn' | wc -c) != "0" ]; then
            samtools faidx "~{reference_fasta}"
          fi

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

        else
          # handle empty case
          cp everything.vcf "~{out_basename}.discordant.vcf"
          echo 0 > num_concordant
          echo 0 > num_discordant_snps
          echo 0 > num_discordant_indels
        fi
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
    parameter_meta {
      seqid_list:{
        description: "List of sequence ID tags for orginal samples.",
        category: "required"
      }
      demux_meta_by_sample_json: {
        description: "JSON file that specifies demux meta",
        cateogry: "required"      
      }
      assembly_meta_tsv: {
        description: "File containing list of assembled reads in table format",
        category: "required"
      }
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

task wgsim {
    meta {
        description: "Generate synthetic Illumina paired-end reads from reference sequences using wgsim via assembly.py simulate_illumina_reads."
    }

    input {
        String coverage_string
        File   reference_fasta
        File?  coverage_bed
        String out_basename = "simulated_reads"

        Int?   read_length
        Int?   outer_distance
        Float? mutation_rate
        Float? indel_fraction
        Float? indel_extend_probability
        Int?   random_seed

        Int    machine_mem_gb = 7
        String docker = "quay.io/broadinstitute/viral-assemble:2.4.3.1"
    }

    parameter_meta {
        coverage_string: {
            description: "Either a simple coverage value (e.g., '10X') or a space-separated string of colon-separated pairs where each pair consists of a GenBank accession and a coverage value (e.g., 'KJ660346.2:12.5x NC_004296.1:0.9X'). Ignored if coverage_bed is provided.",
            category: "required"
        }
        reference_fasta: {
            description: "Reference genome sequences in FASTA format from which to generate reads.",
            patterns: ["*.fasta"],
            category: "required"
        }
        coverage_bed: {
            description: "Optional BED file specifying coverage per region. If provided, this overrides coverage_string.",
            patterns: ["*.bed"],
            category: "advanced"
        }
        read_length: {
            description: "Length of simulated reads (default: 150bp).",
            category: "advanced"
        }
        outer_distance: {
            description: "Outer distance between paired reads (default: 500bp).",
            category: "advanced"
        }
        mutation_rate: {
            description: "Base mutation rate (default: 0.001).",
            category: "advanced"
        }
        indel_fraction: {
            description: "Fraction of mutations that are indels (default: 0.15).",
            category: "advanced"
        }
        indel_extend_probability: {
            description: "Probability an indel is extended (default: 0.3).",
            category: "advanced"
        }
        random_seed: {
            description: "Random seed for reproducibility.",
            category: "advanced"
        }
    }

    Int disk_size = 100

    command <<<
        set -e
        assembly.py --version | tee VERSION

        # Run simulate_illumina_reads
        assembly.py simulate_illumina_reads \
            "~{reference_fasta}" \
            "~{out_basename}.bam" \
            ~{if defined(coverage_bed) then '--coverage_bed "' + coverage_bed + '"' else '"' + coverage_string + '"'} \
            ~{'--read_length ' + read_length} \
            ~{'--outer_distance ' + outer_distance} \
            ~{'--mutation_rate ' + mutation_rate} \
            ~{'--indel_fraction ' + indel_fraction} \
            ~{'--indel_extend_probability ' + indel_extend_probability} \
            ~{'--random_seed ' + random_seed} \
            --loglevel=DEBUG

        # Count reads in output BAM
        samtools view -c "~{out_basename}.bam" | tee READ_COUNT
    >>>

    output {
        File   simulated_reads_bam = "~{out_basename}.bam"
        Int    read_count          = read_int("READ_COUNT")
        String viralngs_version    = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu: 2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
}
