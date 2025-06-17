'''
Download sequences as .fasta from ID list
'''
from Bio import Entrez
from os.path import exists
from os import getenv, remove
import json
from time import sleep
import dotenv

dotenv.load_dotenv()

snk = snakemake # type: ignore
ID_PATH = snk.input[0]
VIOLA_SEQS_PATH = snk.output[0]

with open(ID_PATH, "r") as f:
    id_file = json.loads(f.read())

ids = list(id_file["id_List"])

if exists(VIOLA_SEQS_PATH): remove(VIOLA_SEQS_PATH)

Entrez.email = getenv("email")

for id in ids:
    try:
        stream = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", 
                            retmode="text")
        fasta = stream.read()
        stream.close()
        sleep(3)
        with open(VIOLA_SEQS_PATH, "a") as f:
            f.write(fasta)
        print(f"ID Downloaded: {id}")
    except Exception as e:
        print(f"Download failed: {id}\n")
        raise e

print(f"Consolidated fasta: {VIOLA_SEQS_PATH}")