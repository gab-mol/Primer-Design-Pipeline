'''
Search for conserved sequences for Viola spp. in NCBI-nucleotide.
'''
from Bio import Entrez
from os import getenv
from os.path import join, abspath
from datetime import datetime
import json
import dotenv

dotenv.load_dotenv()
Entrez.email = getenv("email")

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

# Logging
log_name = f"{snk.wildcards.org}_{snk.wildcards.gene}_search.log"
PATH_INFO = abspath(join("..", "log", log_name))

timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
summary = Entrez.read(handle)

with open(PATH_INFO, "a", encoding="UTF-8") as f:
    f.write(f"Search performed at {timestamp}\n")
    f.write(f"Sequences found, {len(ids)}:\n"+"#"*40)
    for entry in summary:
        entry_format = f"""
ID: {entry['Id']}
Title: {entry['Title']}
Length: {entry['Length']}
Update date: {entry['UpdateDate']}
{"-" * 40}"""
        f.write(entry_format)