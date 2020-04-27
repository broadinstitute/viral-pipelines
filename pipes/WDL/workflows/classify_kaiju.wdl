version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kaiju {
    call metagenomics.kaiju
}
