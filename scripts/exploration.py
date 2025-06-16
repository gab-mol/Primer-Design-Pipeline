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

# NCBI search term for Viola matk gene (200 to 2000 bp)
term = 'Viola[Organism] AND matk[Gene] AND ("200"[SLEN] : "2000"[SLEN])'

DIR = os.path.dirname(__file__)
ID_PATH = os.path.join(DIR, "..", "data", "ids_matk_viola_spp.json")

if not os.path.exists(ID_PATH):
    stream = Entrez.esearch(db="nucleotide", term=term, retmax=20)
    record = Entrez.read(stream)
    ids = record["IdList"]

    with open(ID_PATH, "w", encoding="UTF-8") as f:
        f.write(json.dumps({"id_List":list(ids)}))
else:
    print("Already exists:", ID_PATH)

# Save Summary
PATH_INFO = os.path.join(DIR, "..", "data", "viola_search.txt")
if not os.path.exists(PATH_INFO):
    handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
    summary = Entrez.read(handle)

    with open(PATH_INFO, "a", encoding="UTF-8") as f:
        f.write(f"Se encontraron {len(ids)} secuencias:\n"+"#"*40)
        for entry in summary:
            entry_format = f"""
ID: {entry['Id']}
Title: {entry['Title']}
Length: {entry['Length']}
Update date: {entry['UpdateDate']}
{"-" * 40}"""
            f.write(entry_format)
else:
    print("Already exists:", PATH_INFO)

