#!/usr/bin/env bash
export YAMLDIR=../yamls/regression_jul12_normalized_makehdf5/*
export HDF5DIR=../hdf5files/regression_jul12_normalized/

if [ ! -d $HDFDIR ]; then
    mkdir $HDF5DIR
fi

make_hdf5 --yaml_configs $YAMLDIR --output_dir $HDF5DIR
