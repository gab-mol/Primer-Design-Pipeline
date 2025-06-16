from os.path import dirname, join, exists, abspath
import subprocess

DIR = join(dirname(__file__), "..", "data")
SEQS_ALIGN_PATH = join(DIR, "viola_seqs_align.fasta")
PATH_CONS = join(DIR, "viola_cons.fasta")

if not exists(PATH_CONS):
    subprocess.run([
        "em_cons", 
        "-sequence", SEQS_ALIGN_PATH, 
        "-outseq", PATH_CONS, 
        "-name", "viola_matK_consensus"
    ])
else:
    print("Las secuencia consenso ya fue creada:\n", abspath(PATH_CONS))