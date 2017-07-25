#!/usr/bin/env bash
export YAMLDIR=../yamls/atac_xferlearn_jul24_makehdf5/*
export HDF5DIR=../hdf5files/atac_xferlearn_jul24/

if [ ! -d "$HDFDIR" ]; then
    mkdir $HDF5DIR
fi

make_hdf5 --yaml_configs $YAMLDIR --output_dir $HDF5DIR
