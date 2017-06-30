#!/usr/bin/env bash
export YAMLDIR="classificationJun26_makehdf5/0.1subsample/*"
export HDF5DIR="classificationJun26/0.1subsample/"
mkdir ../hdf5files/$HDF5DIR
make_hdf5 --yaml_configs ../yamls/$YAMLDIR --output_dir ../hdf5files/$HDF5DIR
