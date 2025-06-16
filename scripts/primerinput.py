'''
Paso vital para alimentar `primer3_core`.

Eliminar ambiguedades (marcadascomo "n" en consenso).
Generar el archivo input con formato correcto para el 
comando de Primer3.
'''
from Bio import SeqIO
from os.path import dirname, join, exists, abspath

DIR = join(dirname(__file__), "..", "data")
CONS_PATH = join(DIR, "viola_cons.fasta")
PRIM_INPUT_PATH = join(DIR, "viola_primer_input.txt")

registro = SeqIO.read(CONS_PATH, "fasta")
secuencia = str(registro.seq).upper().replace("N", "")  # remover Ns
secuencia = secuencia.replace("\n", "").replace("\r", "")

with open(PRIM_INPUT_PATH, "w") as f:
    f.write(f"""SEQUENCE_ID={registro.id}
SEQUENCE_TEMPLATE={secuencia}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=24
PRIMER_PRODUCT_SIZE_RANGE=100-300
=
""")
