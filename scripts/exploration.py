'''
Search for conserved sequences for Viola spp. in NCBI-nucleotide.
'''
from Bio import Entrez
import os
import json
import dotenv

dotenv.load_dotenv()
Entrez.email = os.getenv("email")

snk = snakemake # type: ignore
ID_PATH = snk.output[0]
search_args = snk.config["entrez"]

genes = search_args["genes"]
organisms = search_args["organisms"]
min_len = search_args["min_len"]
max_len = search_args["max_len"]

# NCBI search term(s)
terms = [
    f"{organism}[Organism] AND {gene}[Gene] AND (\"{min_len}\"[SLEN] : \"{max_len}\"[SLEN])"
    for organism in organisms
    for gene in genes
]

stream = Entrez.esearch(db="nucleotide", term=terms, retmax=20)
record = Entrez.read(stream)
stream.close()
ids = record["IdList"]

with open(ID_PATH, "w", encoding="UTF-8") as f:
    f.write(json.dumps({"id_List":list(ids)}))
