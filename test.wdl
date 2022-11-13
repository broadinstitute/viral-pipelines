#####OVERALL STRUCTURE-----------------------------------------------------------------------------------------
#basic structure for WDL
#call component 'calls' task
workflow test_worfklow { 
    File ref1
    File input1
    String name 

call task_A {
    input: ref= ref1, in= input_1, id=name 
} 
call task_B{
    input:ref= ref1, in=task_A.out
}
}
##WORKFLOW -----------------------------------------------------------------------------------------
##Workflow: required component
## contains "call" action which in turn contains "task": workflow>>call>>>task
workflow myWorkflowName {
    call my_task 
}
##parameters that you can add to workflow: meta, runtime, parameter_meta-- 
#runtime- How long the Docker images runs for 
#parameter_meta- Descriptions of Inputs and Outputs 
#Meta- Author && email 
##EX: 
workflow augur_export_only {
    meta {
        description: "Convert a newick formatted phylogenetic tree with other config settings and node values into a json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

####CALL -----------------------------------------------------------------------------------------
#look at the ways we can type out the call-command
# in its simplest form 
call my_task

# with input variables
call my_task{
    input: task_var1= workflow_var1, task_var2= workflow_var2, ...
}
##with an alias and input variables
call my_task as task_alias {
    input: task_var1= workflow_var1, task_var2= workflow_var2, ...
    }
##TASK -----------------------------------------------------------------------------------------
    ## Task = "To do" fxns 
    ##can be done using parameters

    task my_task {
        [ input definitions ]
        command { ... }
        output { ... }
    }

## the command block = literal command line to run (basically any command that you could otherwise run in a terminal shell) with placeholders (e.g. ${inputfile})
##all variable placeholders MUST be defined in the task input 
#more information on task: 4 parts= command + runtime + output 
workflow myWorkflowName {
call my_task {
        [ input definitions ]
        command { ... }
        output { ... }
}
}
##EX:command {
    java -jar myExecutable.jar \
        INPUT=${input_file} \
        OUTPUT=${output_basename}.txt
}
