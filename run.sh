#!/bin/bash

#Make sure to have activate an python environment (pip install -r requirements.txt)
source env_metadata/bin/activate

#Create needed subdirectories
mkdir -p pkl
mkdir -p logs
mkdir -p camembert

i=0
echo "START ($i / 2869)"
for f in index_data/*; do
	index_name=${f#index_data/}

	python3 metadata.py $index_name > logs/$index_name.log 2>&1

	i=$(($i+1))
	echo "$index_name ($i / 2869)"
done
