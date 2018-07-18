echo "Normalizing data and combining into a single data matrix file"
python normalize_and_combine_data.py

echo "Z-score normalizing labels column-wise"
python zscore_normalize_labels.py

echo "Creating train/val/test splits, regression labels, and sequence features"
python create_splits_labels_features.py

export RUN_NAME=counts_poisson_sep7

echo "Gzipping files"
gzip -f $DL/labels/labels_${RUN_NAME}.txt
gzip -f $DL/splits/${RUN_NAME}/*.txt
gzip -f $DL/features/sequences_${RUN_NAME}.fa

echo "Making HDF5 files from splits, labels, features"
cd $DL/scripts/

export YAMLDIR=../yamls/${RUN_NAME}_makehdf5/*
export HDF5DIR=../hdf5files/${RUN_NAME}/

mkdir $HDF5DIR

make_hdf5 --yaml_configs $YAMLDIR --output_dir $HDF5DIR
