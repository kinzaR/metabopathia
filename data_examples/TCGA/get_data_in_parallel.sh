#!/bin/sh
module load parallel
if [ ! -d logs ]; then
	mkdir logs 
fi
if [ ! -d processed_data ]; then
	mkdir processed_data 
fi
parallel='parallel -j 30 --joblog logs/.runtask_2.log --resume'
$parallel ./get_data_per_cancer.sh '--output_folder' 'processed_data' '--cancer' '{1}' ::: "$@"

### sbatch  --job-name getData --mem=100000  --error '.err_parallel.job' --output '.out_parallel.job' get_data_in_parallel.sh {BLCA,BRCA,COAD,HNSC,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,THCA,UCEC}
# Then I realized that for batch effect they used more data than the mentioned data , so in the next command I will get the missed data only 
### sbatch  --job-name getMissedData --mem=100000  --error '.err_parallel_missed.job' --output '.out_parallel_missed.job' get_data_in_parallel.sh {ESCA,OV,STAD}

