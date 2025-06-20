'''
Search in NCBI-nucleotide. Store IDs in JSON format.
'''
from Bio import Entrez
from os import getenv, mkdir, remove
from os.path import join, abspath, exists
import sys
from datetime import datetime
import json
import dotenv
import yaml
from time import sleep
import subprocess

dotenv.load_dotenv()
Entrez.email = getenv("email")

CONFIG_PATH = join("config", "config.yml")
_CFG_V_PATH = join("config", ".cfg_validated.yml")

founds = 0
validated_gen: list[str] = list()

def search(organism, gene, min_len, max_len, LOG_DIR = "logs"):

    global founds, validated_gen

    ID_PATH = join("data", f"ids_{organism}_{gene}.json")

    # NCBI search term
    terms = f"{organism}[Organism] AND {gene}[Gene] AND (\"{min_len}\"[SLEN] : \"{max_len}\"[SLEN])"
    
    stream = Entrez.esearch(db="nucleotide", term=terms, retmax=20)
    sleep(3)
    record = Entrez.read(stream)
    stream.close()
    ids = record["IdList"]

    results = True if len(ids) > 0 else False

    if results:
        with open(ID_PATH, "w", encoding="UTF-8") as f:
            f.write(json.dumps({"id_List":list(ids)}))

        handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
        summary = Entrez.read(handle)
        founds += 1
        validated_gen.append(gene)
        print("IDs FOUND    --  ", terms)
    else:
        print("NO RESULTS !    --  ", terms)

    # Logging
    timestamp = datetime.now()
    log_name = f"{organism}_{gene}_search_{timestamp.strftime('%Y-%m-%d_%H.%M.%S')}.log"
    PATH_INFO = abspath(join(LOG_DIR, log_name))
    if not exists(LOG_DIR): mkdir(LOG_DIR)

    with open(PATH_INFO, "a", encoding="UTF-8") as f:

        if results:

            f.write(f"Search performed at {timestamp.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Sequences found, {len(ids)}:\n"+"#"*40)
            for entry in summary:
                entry_format = f"""
        ID: {entry['Id']}
        Title: {entry['Title']}
        Length: {entry['Length']}
        Update date: {entry['UpdateDate']}
        {"-" * 40}"""
                f.write(entry_format)

        else:
            f.write(f"Search performed at {timestamp}\n")
            f.write("\n\nNo sequences were found for the search terms !!!\n")
            f.write(f"TERM:\n\t{terms}\n\n")
            f.write("\nNothing to download.\nSkipping...")

if __name__ == "__main__":
    print("""
          
          ##################################
          ##### Primer-Design-Pipeline #####
          ##################################\n"""
    )

    print("Reading configuration (config/config.yml) ...\n")
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    ncor = int(config["cores"])
    entrez_cfg = config["entrez"]
    organism = entrez_cfg["organism"]
    genes = entrez_cfg["genes"]
    min_len = entrez_cfg["min_len"]
    max_len = entrez_cfg["max_len"]

    print("Searching...\n")
    for gen in genes:
        search(organism, gen, min_len, max_len)

    if founds < 1:
        print('''
            
        No sequences were found for any combination of the search terms. !!!
        
            >>>  Aborting Pipeline execution <<<
    ''')
        sys.exit(1)
    else:
        
        config["entrez"]["genes"] = validated_gen
        with open(_CFG_V_PATH, "w", encoding="UTF-8") as f:
            yaml.dump(config, f, sort_keys=False)

        print("\nLaunching Pipeline...")
        stat_code = subprocess.run(["snakemake", "--cores", f"{ncor}"])

        if stat_code.returncode == 0:
            print("\nPipeline completed successfully.")
            if exists(_CFG_V_PATH):
                remove(_CFG_V_PATH)
        else:
            print("Pipeline failed !!!")