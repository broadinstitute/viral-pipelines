import "../tasks/tasks_assembly.wdl" as assembly

workflow assemble_refbased_old {
  call assembly.refine_2x_and_plot
}
