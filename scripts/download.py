'''
Download sequences as .fasta from ID list
'''
from Bio import Entrez
from os.path import dirname, join, exists
from os import getenv, remove
import json
from time import sleep
import dotenv

dotenv.load_dotenv()

DIR = join(dirname(__file__), "..", "data")
ID_PATH = join(DIR, "ids_matk_viola_spp.json")

with open(ID_PATH, "r") as f:
    id_file = json.loads(f.read())

ids = list(id_file["id_List"])
filename = "viola_seqs.fasta"
VIOLA_SEQS_PATH = join(DIR, filename)

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

print(f"Consolidated fasta: {filename}")