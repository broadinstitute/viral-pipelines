version 1.0

import "../tasks/tasks_interhost.wdl" as interhost

workflow trimal {
    call interhost.trimal_clean_msa
}
