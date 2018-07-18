#!/bin/bash

python create_sample_weights.py
python create_individual_label_files.py

cd $DL/scripts/
./makehdf5.sh
./train_functional.sh
