from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

from os.path import dirname, join, exists, abspath

DIR = join(dirname(__file__), "..", "data")
PATH_AMP = join(DIR, "viola_ampl.fasta")
PATH_XML_BLAST = join(DIR, "blast_result.xml")

amplicon = SeqIO.read(PATH_AMP, "fasta")

# Enviamos el amplicón para BLASTn contra nt (puede tardar varios minutos)

# Guardamos los resultados
if not exists(PATH_XML_BLAST):
    result = NCBIWWW.qblast("blastn", "nt", amplicon.seq)
    with open(PATH_XML_BLAST, "w") as handle:
        handle.write(result.read())
        result.close()
else:
    print("\nBúsqueda ya realizada, resultados en:\n", abspath(PATH_XML_BLAST))

# Parseamos los resultados
with open(PATH_XML_BLAST) as result:
    blast_record = NCBIXML.read(result)

for alignment in blast_record.alignments:
    print("Organismo:", alignment.hit_def)
    print("Score:", alignment.hsps[0].score)
    print("E-value:", alignment.hsps[0].expect)
    print("---")
