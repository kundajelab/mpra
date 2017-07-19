#!/usr/bin/env bash
export YAMLDIR=../yamls/sharpr_regression_znormed_jul18_makehdf5/*
export HDF5DIR=../hdf5files/sharpr_regression_znormed_jul18/

if [ ! -d "$HDFDIR" ]; then
    mkdir $HDF5DIR
fi

make_hdf5 --yaml_configs $YAMLDIR --output_dir $HDF5DIR
