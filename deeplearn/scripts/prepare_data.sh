echo "Extracting fragments"
# extract relevant fragments and print them out to a file in $DL/data/
python ./extract_fragments.py

echo "Creating splits and labels"
# create splits and labels in $DL/splits/<runName>/ and $DL/labels/ respectively
python ./create_splits_labels.py

echo "Creating sequence features"
bedtools getfasta -fi ../../../data/hg19.fa -bed ../data/classificationSignals0.05NegativeSubsample.bed -s -fo ../features/sequencesJun26Classification0.05NegativeSubsample.fa
bedtools getfasta -fi ../../../data/hg19.fa -bed ../data/classificationSignals0.1NegativeSubsample.bed -s -fo ../features/sequencesJun26Classification0.1NegativeSubsample.fa
bedtools getfasta -fi ../../../data/hg19.fa -bed ../data/classificationSignals0.2NegativeSubsample.bed -s -fo ../features/sequencesJun26Classification0.2NegativeSubsample.fa
bedtools getfasta -fi ../../../data/hg19.fa -bed ../data/classificationSignals0.4NegativeSubsample.bed -s -fo ../features/sequencesJun26Classification0.4NegativeSubsample.fa

echo "Padding sequences"
python ./pad_sequences.py

echo "Gzipping splits, labels, and features"
gzip ../splits/*Jun26*/*.txt
gzip ../labels/*.txt
gzip ../features/padded*.fa
