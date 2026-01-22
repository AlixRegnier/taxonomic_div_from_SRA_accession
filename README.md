# Get taxonomy from SRA accession

## Requirements

* In ``metadata.py``, you need to put your API key and your NCBI account mail address.
* Make sure your Python environment has all needed packages with:
```bash
#Isolated environment (recommended)
python3 -m venv env_metadata
source env_metadata/bin/activate

#Installing needed packages
pip3 install -r requirements.txt
```

* Currently the FOF hierarchy is the following
```
index_data/
├── GENOMIC_BCT_10
│   └── kmtricks.fof
├── GENOMIC_BCT_11
│   └── kmtricks.fof
├── GENOMIC_BCT_12
│   └── kmtricks.fof
└── GENOMIC_BCT_13
    └── kmtricks.fof
```
## Getting accessions taxonomy from a single FOF in ``index_data``

```bash
mkdir -p pkl camembert
python3 metadata.py GENOMIC_BCT_10
```

## Getting accessions taxonomy from all FOF in ``index_data``
```bash
./run.sh
``` 

## Getting camembert (pie chart)
For each FOF, a simple pie chart is generated from accessions division in ``./camembert`` 