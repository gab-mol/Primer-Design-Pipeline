'''
Search for conserved sequences for Viola spp. in NCBI-nucleotide.

Selected sequence: "matK" (Maturase K)
Plastidial gene for intron maturase, a protein essential for the in vivo splicing of Group II introns.
'''
from Bio import Entrez
import os
import json
import dotenv

dotenv.load_dotenv()
Entrez.email = os.getenv("email")

snk = snakemake # type: ignore
ID_PATH = snk.output[0]

# NCBI search term
term = 'Viola[Organism] AND matk[Gene] AND ("200"[SLEN] : "2000"[SLEN])'
stream = Entrez.esearch(db="nucleotide", term=term, retmax=20)
record = Entrez.read(stream)
stream.close()
ids = record["IdList"]

with open(ID_PATH, "w", encoding="UTF-8") as f:
    f.write(json.dumps({"id_List":list(ids)}))
