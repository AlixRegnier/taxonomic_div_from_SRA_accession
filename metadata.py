#!/usr/bin/env python3
"""
Retrieve SRA metadata using Entrez to plot taxonomic division
Currently, only corresponding ENA taxonomic division and organism are stored

accession -> taxid -> { div, organism }

Requires: pip install ratelimit,biopython,matplotlib
"""

import subprocess
import json
import sys
import pickle
import re
import concurrent.futures
from typing import List
from collections import OrderedDict
import urllib.error
import matplotlib.pyplot as plt
from ratelimit import limits, sleep_and_retry

import traceback

from Bio import Entrez

Entrez.api_key = "4459a438546c82f6e448b428a73e72eb0908"
Entrez.email = "alix.regnier@inria.fr"

MAX_RATE = 9

class Div:
    BCT = "BCT"
    ENV = "ENV"
    HUMAN = "HUMAN"
    INV = "INV"
    MAM = "MAM"
    MICE = "MICE"
    PHG = "PHG"
    PLN = "PLN"
    PRI = "PRI"
    ROD = "ROD"
    SYN = "SYN"
    UNKNOWN = "UNKNOWN"
    VRL = "VRL"
    VRT = "VRT"
    UNRESOLVED = "UNRESOLVED"

DIV_LOGANDIV = {
    "BCT" :  Div.BCT,
    "ENV" :  Div.ENV,
    "HUM" :  Div.HUMAN,
    "INV" :  Div.INV,
    "MAM" :  Div.MAM,
    "MICE" : Div.MICE,
    "PHG" :  Div.PHG,
    "PLN" :  Div.PLN,
    "PRI" :  Div.PRI,
    "ROD" :  Div.ROD,
    "SYN":   Div.SYN,
    "VRL" :  Div.VRL,
    "VRT" :  Div.VRT,
    "MUS" :  Div.MICE,
    #UNKNOWN
    "PRO":   Div.UNKNOWN, #No corresponding group in Logan (Archea/Bact)
    "UNC":   Div.UNKNOWN, #Unclassified
    "UNA":   Div.UNKNOWN, #Unannotated
    "FUN":   Div.UNKNOWN, #No corresponding group in Logan (Fungi)
    "ARC":   Div.UNKNOWN  #No corresponding group in Logan (Archea)
}

UNRESOLVED_TAXID = -2 #Error couldn't be handled
UNKNOWN_TAXID = -1    #Was probably deleted from SRA

EXCEPTION_DIV = {
    UNRESOLVED_TAXID : Div.UNRESOLVED,
    UNKNOWN_TAXID    : Div.UNKNOWN,
    9604             : Div.HUMAN, #Hominidae
    10066            : Div.MICE   #Muridae
}

accession_logandiv = dict()


#Load dictionaries
print(":: Loading dicts...")
try:
    with open("pkl/taxid_div.pkl", "rb") as f:
        known_taxid_div = pickle.load(f)
except:
    known_taxid_div = { 9604: Div.HUMAN, 10066: Div.MICE }

try:
    with open("pkl/taxid_organism.pkl", "rb") as f:
        known_taxid_organism = pickle.load(f)
except:
    known_taxid_organism = { 9604: "Hominidae", 10066: "Muridae"}

try:
    with open("pkl/accession_taxid.pkl", "rb") as f:
        known_accession_taxid = pickle.load(f)
except:
    known_accession_taxid = dict()

def ex(d, depth=0):
    if type(d) is dict:
        for k in d:
            print(depth*"\t", k, sep="")
            ex(d[k], depth+1)
    elif type(d) is list:
        for i in range(len(d)):
            print(depth*"\t", f"[{i}]", sep="")
            ex(d[i], depth+1)
    else:
        print(depth*"\t", f"{type(d)}: {d}")

#import tabulate
#d = [(taxid, known_taxid_organism[taxid], known_taxid_div[taxid]) for taxid in known_taxid_div]
#print(tabulate.tabulate(d))

RE_EFETCH = re.compile(r"<TAXON_ID>([0-9]+)</TAXON_ID>")

@sleep_and_retry
@limits(calls=MAX_RATE, period=1.0)
def get_taxids_from_accessions(accessions: List[str], db="sra", batch_size=1000) -> List[int]:
    taxids = OrderedDict()
    unresolved_accessions = []

    for accession in accessions:
        if accession in known_accession_taxid:
            taxids[accession] = known_accession_taxid[accession]
        else: #By default, accession taxid is set to UNRESOLVED
            taxids[accession] = UNRESOLVED_TAXID
            known_accession_taxid[accession] = UNRESOLVED_TAXID
            unresolved_accessions.append(accession)

    #Resolve only accessions that were never seen
    for i in range(0, len(unresolved_accessions), batch_size):
        batch = unresolved_accessions[i:i+batch_size]
        handle = None

        try:
            if len(batch) == 1:
                handle = Entrez.efetch(db=db, id=batch[0])
            else:
                handle = Entrez.efetch(db=db, id=",".join(batch), retmax=batch_size)

        except urllib.error.HTTPError as e:
            if e.code == 400:
                for accession in batch:
                    taxids[accession] = UNKNOWN_TAXID
                    known_accession_taxid[accession] = UNKNOWN_TAXID
            else:
                raise

        response:str = handle.read().decode("utf-8")
        handle.close()

        matches = RE_EFETCH.findall(response)
        
        #If NCBI does shit (when it doesn't return all queried accessions metadata)
        if len(matches) == 0:
            pass
        elif len(matches) != len(batch):
            pos=0 #Search starting index
            j = 0 #Index on batch
            k = 0 #Index on matches

            while j < len(batch):
                str_index = response.find(batch[j], pos)
                if str_index != -1:
                    taxids[batch[j]] = int(matches[k])
                    known_accession_taxid[batch[j]] = int(matches[k])
                    
                    pos = str_index+1
                    k += 1
                j += 1
        else:
            for j in range(len(batch)):
                taxids[batch[j]] = int(matches[j])
                known_accession_taxid[batch[j]] = int(matches[j])

    return [taxids[accession] for accession in accessions]

def get_division_from_taxid(taxid: int) -> dict:

    #Check if taxid has an overriden division
    if taxid in EXCEPTION_DIV:
        return EXCEPTION_DIV[taxid]
    
    #Check if corresponding division is known
    if taxid in known_taxid_div:
        return known_taxid_div[taxid]
    
    try:
        result = subprocess.run(
            ["curl", "https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/" + str(taxid)],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse JSON output
        j = json.loads(result.stdout)

        #Retrieve scientific name, division
        known_taxid_div[taxid] = DIV_LOGANDIV[j["division"]]
        known_taxid_organism[taxid] = j["scientificName"]

    except:
        known_taxid_div[taxid] = Div.UNRESOLVED
        known_taxid_organism[taxid] = "UNKNOWN"
    finally:
        return known_taxid_div[taxid]


def get_divisions_from_accessions(accessions: List[str], batch_size=1000) -> str:
    divisions = [None]*len(accessions)

    for i, taxid in enumerate(get_taxids_from_accessions(accessions, batch_size=batch_size)):
       divisions[i] = get_division_from_taxid(taxid)

    return divisions



def camembert(div: str, accessions: List[str], filename: str = None):
    divs = { get_division_from_taxid(known_accession_taxid[accession]) for accession in accessions}
    divs.add(div)

    div_counts = dict(zip(divs, [0]*len(divs)))

    for accession in accessions:
        div = get_division_from_taxid(known_accession_taxid[accession])
        div_counts[div] += 1

    category = list(divs)
    value = [div_counts[div] for div in divs]

    _, ax = plt.subplots()
    ax.pie(value, labels=category, autopct="%1.1f%%")

    if filename is None:
        plt.show()
    else:
        plt.savefig(f"camembert/{filename}", dpi=300)


def get_index_div(index_name: str) -> str:
    a = index_name.find("_")
    b = index_name.find("_", a+1)

    if a == -1 or b == -1:
        raise Exception("Invalid index name, should be only like 'GENOMIC_BCT_24*', (*) anything.")
    
    return index_name[a+1:b]

def main():

    if len(sys.argv) != 2:
        print("Usage: python metadata.py <index_name>")
        exit(2)
    
    index_name = sys.argv[1]
    index_div = get_index_div(index_name)
    fof = f"index_data/{index_name}/kmtricks.fof"

    accessions = []

    print(f":: Retrieving accessions from {index_name}...")
    with open(fof, "r") as f:
        for line in f:
            accession = line[:line.find(':')]
            accessions.append(accession)
            accession_logandiv[accession] = index_div

    #Get unknown accessions
    unknown_accessions = list(set(accessions) - known_accession_taxid.keys())

    processed = len(accessions) - len(unknown_accessions)
    batch_size = 100

    if len(unknown_accessions) == 0:
        print(f":: All accessions were found in dicts.", end="")
    else:
        print(f":: Retrieving accessions metadata from NCBI (batch size: {batch_size})...")
        print(f"\r\t{processed} / {len(accessions)} ({int(processed/len(accessions)*100)}%)", end="")

    try:
        futures_batches = OrderedDict()
        with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_RATE) as executor:
            for i in range(0, len(unknown_accessions), batch_size):
                batch = unknown_accessions[i:i+batch_size]
                future = executor.submit(get_divisions_from_accessions, batch, batch_size)
                futures_batches[future] = batch[:]

            for future in concurrent.futures.as_completed(futures_batches):
                try:
                    future.result() #Discard
                except urllib.error.HTTPError as e:
                    if e.code == 429:
                        print(f"\nToo many requests were done to NCBI.")
                    else:
                        print(traceback.format_exc(), end="\n\n")
                except Exception as e:
                    print(f"\nSomething went wrong: {e}")
                    print(traceback.format_exc(), end="\n\n")
                finally:
                    processed += len(futures_batches[future])
                    print(f"\r\t{processed} / {len(accessions)} ({int(processed/len(accessions)*100)}%)", end="")

    except KeyboardInterrupt:
        pass
    except Exception as e:
        with open("err.txt", "a") as f:
            f.write(f"#Index: {fof}, error: {e}\n")
    finally:
        if len(unknown_accessions) == 0:
            print("\n:: Dicts don't need to be updated.")
        else:
            print("\n:: Saving dicts...")
            with open("pkl/taxid_div.pkl", "wb") as f:
                pickle.dump(known_taxid_div, f)

            with open("pkl/taxid_organism.pkl", "wb") as f:
                pickle.dump(known_taxid_organism, f)

            with open("pkl/accession_taxid.pkl", "wb") as f:
                pickle.dump(known_accession_taxid, f)
        
    camembert(index_div, accessions, f"{index_name}.png")

if __name__ == "__main__":
    main()