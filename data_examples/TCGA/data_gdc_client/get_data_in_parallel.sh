#!/bin/sh
module load parallel
parallel='parallel -j 30 --joblog logs/.runtask.log --resume'
$parallel 'get_data_per_cancer.sh' '--output_folder' 'processed_data' '--cancer' '{1}' ::: "$@"