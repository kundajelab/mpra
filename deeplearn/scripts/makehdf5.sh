#!/usr/bin/env bash
export YAMLDIR=../yamls/regression_jul3_strands_augmented_makehdf5/*
export HDF5DIR=../hdf5files/regression_jul3_strands_augmented_rep2only/

if [ ! -d $HDFDIR ]; then
    mkdir $HDF5DIR
fi

make_hdf5 --yaml_configs $YAMLDIR --output_dir $HDF5DIR
