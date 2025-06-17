'''
Prepare input file for `primer3_core`.
Disambiguate (n) in sequence.
'''
from Bio import SeqIO

snk = snakemake # type: ignore
CONS_PATH = snk.input[0]
PRIM_INPUT_PATH = snk.output[0]

cons = SeqIO.read(CONS_PATH, "fasta")
sequence = str(cons.seq).upper().replace("N", "")
sequence = sequence.replace("\n", "").replace("\r", "")

with open(PRIM_INPUT_PATH, "w") as f:
    f.write(f"""SEQUENCE_ID={cons.id}
SEQUENCE_TEMPLATE={sequence}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=24
PRIMER_PRODUCT_SIZE_RANGE=100-300
=
""")
