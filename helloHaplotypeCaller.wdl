#test WDL file to create workflow called helloHaplotypeCaller that consist of a single task that calls a pre-exisitng GATK's HaplotypeCaller. 
#purpose to perform variant discovery on high-throughput sequencing data. 

#declare your workflow
#"task" get declared by using the 'call' action
workflow helloHaplotypeCaller {
    call HaplotypeCaller
}

