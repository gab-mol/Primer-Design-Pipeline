configfile: "config/.cfg_validated.yml"

def ensure_list(x):
    return x if isinstance(x, list) else [x]

terms = config["entrez"]
GENES = ensure_list(terms["genes"])
ORGANISM = ensure_list(terms["organism"])

rule all:
    input:
        expand("data/{org}_{gene}_primers.txt", org=ORGANISM, gene=GENES)

rule download:
    output:
        "data/{org}_{gene}_seqs.fasta"
    input: 
        "data/ids_{org}_{gene}.json"
    script:
        "scripts/download.py"

rule alignment:
    output:
        "data/{org}_{gene}_seqs_align.fasta"
    input:
        "data/{org}_{gene}_seqs.fasta"
    shell:
        "clustalo -i {input} -o {output} --auto --verbose"

rule consensus:
    output:
        "data/{org}_{gene}_cons.fasta"
    input:
        "data/{org}_{gene}_seqs_align.fasta"
    shell:
        "em_cons -sequence {input} -outseq {output} -name {wildcards.org}_{wildcards.gene}_consensus"

rule primerinput:
    output:
        "data/{org}_{gene}_primer_input.txt"
    input:
        "data/{org}_{gene}_cons.fasta"
    script:
        "scripts/primerinput.py"

rule primer:
    output:
        "data/{org}_{gene}_primers.txt"
    input:
        "data/{org}_{gene}_primer_input.txt"
    shell:
        "primer3_core < {input} > {output}"