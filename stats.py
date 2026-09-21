#!/usr/bin/env python3
"""
Retrieve SRA metadata using Entrez to plot taxonomic division
Currently, only corresponding ENA taxonomic division and organism are stored

accession -> taxid -> { div, organism }

Requires: pip3 install ratelimit,biopython,matplotlib
Having installed: curl
"""

import subprocess
import json
import sys
import pickle
import os

from random import sample
from typing import List
from collections import OrderedDict
import urllib.error
import matplotlib.pyplot as plt
from typing import Dict, List
import pandas as pd
import seaborn

import traceback

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
    TIMEOUT = "TIMEOUT"

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
    "PRO":   Div.UNKNOWN, #See in get_division_from_taxid(), will be overriden to check if BCT
    
    #UNKNOWN
    "UNC":   Div.UNKNOWN, #Unclassified
    "UNA":   Div.UNKNOWN, #Unannotated
    "FUN":   Div.UNKNOWN, #No corresponding group in Logan (Fungi)
    "ARC":   Div.UNKNOWN  #No corresponding group in Logan (Archea)
}

UNKNOWN_TAXID = -1    #Was probably deleted from SRA
UNRESOLVED_TAXID = -2 #Error couldn't be handled
TIMEOUT_TAXID = -3

BAD_DIV_CODES = (Div.UNKNOWN, Div.UNRESOLVED, Div.TIMEOUT)

BAD_TAXID_CODES = dict(zip(
    (UNKNOWN_TAXID, UNRESOLVED_TAXID, TIMEOUT_TAXID), 
    BAD_DIV_CODES
))

EXCEPTION_DIV = {
    9604  : Div.HUMAN, #Hominidae
    10066 : Div.MICE   #Muridae
}

#Load dictionaries
print(":: Loading dicts...")
try:
    with open("pkl/taxid_div.pkl", "rb") as f:
        known_taxid_div = pickle.load(f)
except:
    known_taxid_div = EXCEPTION_DIV.copy()

try:
    with open("pkl/taxid_organism.pkl", "rb") as f:
        known_taxid_organism = pickle.load(f)
except:
    known_taxid_organism = EXCEPTION_DIV.copy()

try:
    with open("pkl/accession_taxid.pkl", "rb") as f:
        known_accession_taxid = pickle.load(f)
except:
    known_accession_taxid = dict()

def get_index_tech(index_name : str) -> str:
    a = index_name.find("_")

    if a == -1:
        raise Exception("Invalid index name, should be only like 'GENOMIC_BCT_24*', (*) anything.")

    return index_name[:a]

def get_index_div(index_name: str) -> str:
    a = index_name.find("_")
    b = index_name.find("_", a+1)

    if a == -1 or b == -1:
        raise Exception("Invalid index name, should be only like 'GENOMIC_BCT_24*', (*) anything.")
    
    return index_name[a+1:b]

def get_index_span(index_name: str) -> str:
    a = index_name.find("_")
    b = index_name.find("_", a+1)
    c = index_name.find("_", b+1)

    if b == -1:
        raise Exception("Invalid index name, should be only like 'GENOMIC_BCT_24*', (*) anything.")

    if c == -1:
        return index_name[b+1:]
    return index_name[b+1:c]


def main():
    MAX_SPAN = 40

    divs = frozenset([  
        Div.BCT,
        Div.ENV,
        Div.HUMAN,
        Div.INV,
        Div.MAM,
        Div.MICE,
        Div.PHG,
        Div.PLN,
        Div.PRI,
        Div.ROD,
        Div.SYN,
        Div.VRL,
        Div.VRT,
        Div.UNKNOWN,
        Div.UNRESOLVED,
        Div.TIMEOUT,
    ])

    #Chain: DIV -> DIV -> COUNT
    div_predicted_div = { div : dict(zip(divs, [0]*len(divs))) for div in (divs) }

    #Chain: SPAN -> DIV -> DIV -> COUNT
    span_div_predicted_div = { str(i) : { div : dict(zip(divs, [0]*len(divs))) for div in (divs) } for i in range(MAX_SPAN) }

    #Chain: SPAN -> DIV -> COUNT
    span_predicted_div_count = { str(i) : dict(zip(divs, [0]*len(divs))) for i in range(MAX_SPAN) }

    index_div_count = dict(zip(divs, [0]*(len(divs))))
    divs = set()

    for index_name in os.listdir("index_data/"):

        if get_index_tech(index_name) not in {'TRANSCRIPTOMIC', 'SYNTHETIC', 'GENOMICSINGLECELL', 'TRANSCRIPTOMICSINGLECELL', 'VIRALRNA', 'OTHER', 'GENOMIC'}:
            continue

        if get_index_div(index_name) == "UNKNOWN":
            continue

        accessions = []
        fof = f"index_data/{index_name}/kmtricks.fof"
        with open(fof, "r") as f:
            for line in f:
                accession = line[:line.find(':')]

                if accession in known_accession_taxid:
                    accessions.append(accession)

        div = get_index_div(index_name)
        span = get_index_span(index_name)

        index_div_count[div] += 1

        #Count
        for accession in accessions:
            taxid = known_accession_taxid[accession]

            predicted_div = Div.UNKNOWN
            if taxid in known_taxid_div:
                predicted_div = known_taxid_div[taxid]

            if predicted_div not in BAD_DIV_CODES:
                div_predicted_div[div][predicted_div] += 1

            span_predicted_div_count[span][predicted_div] += 1

            span_div_predicted_div[span][div][predicted_div] += 1

            divs.add(predicted_div)

    for div in div_predicted_div:
        s = 0
        for predicted_div in div_predicted_div:
            s += div_predicted_div[div][predicted_div]

        for predicted_div in div_predicted_div:
            div_predicted_div[div][predicted_div] = round(div_predicted_div[div][predicted_div] / max(s, 1), 3)
    
    df = pd.DataFrame.from_dict(div_predicted_div, orient="index")

    # Ensure rows and columns use the same ordering
    fields = sorted((set(div_predicted_div) | {k for row in div_predicted_div.values() for k in row}) - set(BAD_DIV_CODES))
    df = df.reindex(index=fields, columns=fields, fill_value=0)

    seaborn.heatmap(df, annot=True, fmt="g", cmap="Blues")
    plt.xlabel("Predicted taxonomic division")
    plt.ylabel("Taxonomic division")
    plt.show()

    for span in range(MAX_SPAN):
        span = str(span)
        for div in divs:
            total = 0
            total_without_errors = 0

            for predicted_div in divs - set(BAD_DIV_CODES):
                total += span_div_predicted_div[span][div][predicted_div]
                total_without_errors += span_div_predicted_div[span][div][predicted_div]

            for predicted_div in set(BAD_DIV_CODES):
                total += span_div_predicted_div[span][div][predicted_div]

            total_without_errors = max(1, total_without_errors)
            total = max(1, total)

            span_div_predicted_div[span][div]["accuracy%"] = min((span_div_predicted_div[span][div][div]) / total_without_errors * 100.0, 100.0)
            span_div_predicted_div[span][div]["accuracy_with_errors%"] = min((span_div_predicted_div[span][div][div]) / total * 100.0, 100.0)
            span_div_predicted_div[span][div]["unknown%"] = min((span_div_predicted_div[span][div][Div.UNKNOWN]) / total * 100.0, 100.0)
            span_div_predicted_div[span][div]["unresolved%"] = min((span_div_predicted_div[span][div][Div.UNRESOLVED]) / total * 100.0, 100.0)
            span_div_predicted_div[span][div]["timeout%"] = min((span_div_predicted_div[span][div][Div.TIMEOUT]) / total * 100.0, 100.0)

    lspan = []
    laccuracy = []
    laccuracy_with_errors = []
    ldiv = []
    lunresolved = []
    lunknown = []
    ltimeout = []

    for span in range(MAX_SPAN):
        span = str(span)
        for div in { Div.PLN, Div.BCT, Div.HUMAN, Div.MICE}: #divs:
            lspan.append(span)
            ldiv.append(div)
            laccuracy.append(span_div_predicted_div[span][div]["accuracy%"])
            laccuracy_with_errors.append(span_div_predicted_div[span][div]["accuracy_with_errors%"])
            lunknown.append(span_div_predicted_div[span][div]["unknown%"])
            lunresolved.append(span_div_predicted_div[span][div]["unresolved%"])
            ltimeout.append(span_div_predicted_div[span][div]["timeout%"])

    df = pd.DataFrame({
        "span" : lspan,
        "accuracy" : laccuracy,
        "accuracy_with_errors" : laccuracy_with_errors,
        "div" : ldiv,
        "timeout" : ltimeout,
        "unresolved" : lunresolved,
        "unknown" : lunknown
    })

    df["span"] = df["span"].astype(int)
    df = df.sort_values(["div", "span"])

    seaborn.lineplot(
        data=df,
        x=df["span"].astype(str),
        y="accuracy",
        hue="div",
        marker="o"
    )
    plt.xlabel("Span")
    plt.ylabel("Accuracy (%)")
    plt.title("Accuracy of taxonomic division classification according to span and div (errors are skipped)")
    plt.show()

    seaborn.lineplot(
            data=df,
            x=df["span"].astype(str),
            y="accuracy_with_errors",
            hue="div",
            marker="o"
        )
    plt.xlabel("Span")
    plt.ylabel("Accuracy (%)")
    plt.title("Accuracy of taxonomic division classification according to span and div (errors are kept)")
    plt.show()

    seaborn.lineplot(
            data=df,
            x=df["span"].astype(str),
            y="unknown",
            hue="div",
            marker="o"
        )
    plt.xlabel("Span")
    plt.ylabel("Unknown (%)")
    plt.title("Proportion of unknown taxonomic division according to span and index div")
    plt.show()

    seaborn.lineplot(
                data=df,
                x=df["span"].astype(str),
                y="timeout",
                hue="div",
                marker="o"
            )
    plt.xlabel("Span")
    plt.ylabel("Timeout (%)")
    plt.title("Proportion of timeout taxonomic division according to span and index div")
    plt.show()

    seaborn.lineplot(
                data=df,
                x=df["span"].astype(str),
                y="unresolved",
                hue="div",
                marker="o"
            )
    plt.xlabel("Span")
    plt.ylabel("Unresolved (%)")
    plt.title("Proportion of unresolved taxonomic division according to span and index div")
    plt.show()


if __name__ == "__main__":
    main()
