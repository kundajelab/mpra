#!/usr/bin/env bash
export YAMLDIR=../yamls/upweightends_aug10_makehdf5/*
export HDF5DIR=../hdf5files/upweightends_aug10/

if [ ! -d "$HDFDIR" ]; then
    mkdir $HDF5DIR
fi

make_hdf5 --yaml_configs $YAMLDIR --output_dir $HDF5DIR
