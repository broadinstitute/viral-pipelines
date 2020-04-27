version 1.0

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports

workflow multi_sample_assemble_kraken {

    input {
        Array[File]+ raw_uBAMs
        Array[File]+ reference_genome_fasta
        File spikein_db
        File trim_clip_db
        Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? bwaDbs
    }
    
    scatter(raw_reads in raw_uBAMs) {
        call reports.spikein_report as spikein {
            input:
                reads_bam = raw_reads,
                spikein_db = spikein_db
        }
        call taxon_filter.deplete_taxa as deplete {
            input:
                raw_reads_unmapped_bam = raw_reads,
                bmtaggerDbs = bmtaggerDbs,
                blastDbs = blastDbs,
                bwaDbs = bwaDbs
        }
        call assembly.assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = deplete.cleaned_bam,
                trim_clip_db = trim_clip_db,
                always_succeed = true
        }
        call assembly.scaffold {
            input:
                contigs_fasta = spades.contigs_fasta,
                reads_bam = deplete.cleaned_bam,
                reference_genome_fasta = reference_genome_fasta
        }
        call assembly.refine_2x_and_plot {
            input:
                assembly_fasta = scaffold.scaffold_fasta,
                reads_unmapped_bam = deplete.cleaned_bam
        }        
    }

    call metagenomics.krakenuniq as krakenuniq {
        input:
            reads_unmapped_bam = raw_uBAMs
    }

    call reports.spikein_summary as spike_summary {
        input:
            spikein_count_txt = spikein.report
    }

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }
}
