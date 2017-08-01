#!/usr/bin/env bash
export YAMLDIR=../yamls/promoter_tasks_jul31_makehdf5/*
export HDF5DIR=../hdf5files/minP_jul31/

if [ ! -d "$HDFDIR" ]; then
    mkdir $HDF5DIR
fi

make_hdf5 --yaml_configs $YAMLDIR --output_dir $HDF5DIR
