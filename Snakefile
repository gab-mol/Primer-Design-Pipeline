configfile: "config/config.yml"

GENES = config["genes"]
GENERA = config["genera"]

rule all:
    input:
        expand("data/{gen}_{gene}_primers.txt", gen=GENERA, gene=GENES)

rule exploration:
    output: 
        "data/ids_{gen}_{gene}.json"
    script:
        "scripts/exploration.py"

rule download:
    output:
        "data/{gen}_{gene}_seqs.fasta"
    input: 
        "data/ids_{gen}_{gene}.json"
    script:
        "scripts/download.py"

rule alignment:
    output:
        "data/{gen}_{gene}_seqs_align.fasta"
    input:
        "{gen}_{gene}_seqs.fasta"
    shell:
        "clustalo -i {input} -o {output} --auto --verbose"

rule consensus:
    output:
        "data/{gen}_{gene}_cons.fasta"
    input:
        "data/{gen}_{gene}_seqs_align.fasta"
    shell:
        "em_cons -sequence {input} -outseq {output} -name {wildcards.gen}_{wildcards.gene}_consensus"

rule primerinput:
    output:
        "data/{gen}_{gene}_primer_input.txt"
    input:
        "data/{gen}_{gene}_cons.fasta"
    script:
        "scripts/primerinput.py"

rule primer:
    output:
        "data/{gen}_{gene}_primers.txt"
    input:
        "data/{gen}_{gene}_primer_input.txt"
    shell:
        "primer3_core < {input} > {output}"