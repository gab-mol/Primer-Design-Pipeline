'''
El mejor par de primers es el primer0 (ver "primers.txt)
valores de PENALTY más bajos.
PRIMER_PAIR_0_PENALTY=0.428421
PRIMER_LEFT_0_PENALTY=0.253347
PRIMER_RIGHT_0_PENALTY=0.175073

Secuencias:
PRIMER_LEFT_0_SEQUENCE=TCCAAGCATTCCCTCTCCCT
PRIMER_RIGHT_0_SEQUENCE=ATCAGCCCGAGTCGGTTTAC

Posiciones, largos:
PRIMER_LEFT_0=613,20
PRIMER_RIGHT_0=828,20
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os.path import dirname, join, exists, abspath

DIR = join(dirname(__file__), "..", "data")
PATH_CONS = join(DIR, "viola_cons.fasta")
PATH_AMP = join(DIR, "viola_ampl.fasta")

# Cargar la secuencia (puede venir del archivo fasta que usaste)
record = SeqIO.read(PATH_CONS, "fasta")
sequence = str(record.seq)

# Datos de primers
left_start = 613
left_len = 20
right_start = 828
right_len = 20

left_primer = sequence[left_start:left_start + left_len]
right_primer = sequence[right_start:right_start + right_len]

amplicon = sequence[left_start:right_start + right_len]

# Mostrar amplicón con primers resaltados
print("Amplicón:")
print(amplicon)

print("\nPrimer forward (5'→3'):", left_primer)
print("Primer reverse (5'→3'):", right_primer[::-1])  # invertir para visualización

print("\nPosiciones:")
print(f"Forward primer: {left_start}–{left_start+left_len-1}")
print(f"Reverse primer: {right_start}–{right_start+right_len-1}")

if not exists(PATH_AMP):
    record = SeqRecord(
        Seq(amplicon),
        id="viola_amplicon",
        description="PCR product using best primer pair"
    )

    with open(PATH_AMP, "w") as handle:
        SeqIO.write(record, handle, "fasta")
else:
    print("\nAmplicón ya se halla guardado en:\n", abspath(PATH_AMP))