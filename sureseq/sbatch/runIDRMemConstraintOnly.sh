#!/bin/bash
#SBATCH -t 48:00:00                    # Runtime in HH:MM:SS
#SBATCH -p owners                      # Partition to submit to
#SBATCH --mem=512G                     # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o IDRsbatchJune20_%j.out      # File to which STDOUT will be written
#SBATCH -e IDRsbatchJune20_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=therajiv@gmail.com # Email to which notifications will be sent

idr --samples ../fixedIdrSureFragmentSignalMappedReplicate1.bed ../fixedIdrSureFragmentSignalMappedReplicate2.bed --input-file-type bed --rank 5 --plot -o idrSingleBpSbatchRun.txt
