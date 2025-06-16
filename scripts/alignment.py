from os.path import dirname, join, exists, abspath
import subprocess
from Bio import AlignIO

DIR = join(dirname(__file__), "..", "data")
SEQS_PATH = join(DIR, "viola_seqs.fasta")
SEQS_ALIGN_PATH = join(DIR, "viola_seqs_align.fasta")

# Instalado con apt: clustalo (Clustal Omega)
if not exists(SEQS_ALIGN_PATH):
    command = [
        "clustalo",
        "-i", SEQS_PATH,
        "-o", SEQS_ALIGN_PATH,
        "--auto", "--verbose"
    ]
    subprocess.run(command, check=True)
else:
    print("secuencias ya alineadas en:", abspath(SEQS_ALIGN_PATH))

alignment = AlignIO.read(SEQS_ALIGN_PATH, "fasta")
print(f"{len(alignment)} secuencias alineadas")
print(f"Longitud del alineamiento: {alignment.get_alignment_length()}")