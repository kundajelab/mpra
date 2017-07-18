echo "Normalizing data and combining into a single data matrix file"
python normalize_and_combine_data.py

echo "Creating train/val/test splits, regression labels, and sequence features"
python create_splits_labels_features.py

export RUN_NAME=sharpr_regression_jul17

echo "Gzipping files"
gzip -f $DL/labels/labels_$RUN_NAME.txt
gzip -f $DL/splits/$RUN_NAME/*.txt
gzip -f $DL/features/sequences_$RUN_NAME.fa

echo "Making HDF5 files from splits, labels, features"
cd $DL/scripts/
./makehdf5.sh
