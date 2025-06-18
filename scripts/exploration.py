'''
Search in NCBI-nucleotide. Store IDs in JSON format.
'''
from Bio import Entrez
from os import getenv, mkdir
from os.path import join, abspath, exists
from datetime import datetime
import json
import dotenv

dotenv.load_dotenv()
Entrez.email = getenv("email")

snk = snakemake # type: ignore
ID_PATH = snk.output[0]
search_args = snk.config["entrez"] 

organism = snk.wildcards.org
gene = snk.wildcards.gene
min_len = int(search_args["min_len"])
max_len = int(search_args["max_len"])

# NCBI search term
terms = f"{organism}[Organism] AND {gene}[Gene] AND (\"{min_len}\"[SLEN] : \"{max_len}\"[SLEN])"

stream = Entrez.esearch(db="nucleotide", term=terms, retmax=20)
record = Entrez.read(stream)
stream.close()
ids = record["IdList"]

with open(ID_PATH, "w", encoding="UTF-8") as f:
    f.write(json.dumps({"id_List":list(ids)}))

# Logging
log_name = f"{snk.wildcards.org}_{snk.wildcards.gene}_search.log"
LOG_DIR = "logs"
PATH_INFO = abspath(join(LOG_DIR, log_name))
if not exists(LOG_DIR): mkdir(LOG_DIR)

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