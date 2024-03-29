## Bash script detects non-exchangeability for each population in 1KG dataset
## Created by Alan Aw on 7/26/2021
DIR_TO_POP_ID_LIST='pop_id_list.txt'
(while read p; do
    echo "Checking exchangeability of individuals from $p"
    Rscript test_exchangeability.R /Users/alanaw/Documents/research/pgs/030121/pop_files/$p/1000G_$p /Users/alanaw/Documents/research/pgs/101220/plink2 1000
done < $DIR_TO_POP_ID_LIST) > output_log.txt
