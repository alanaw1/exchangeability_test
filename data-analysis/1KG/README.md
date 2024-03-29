# Detecting Exchangeability from the Terminal

## For R Users

This directory contains scripts and data files used to run our test of exchangeability from the terminal. As described in our pedagogical vignette, this demonstration involves the use of the [1000 Genomes Phase 3 dataset](https://www.cog-genomics.org/plink/2.0/resources). 

We encourage interested users to repurpose our scripts for their own analyses. In particular, the R script, written to be executable from the terminal, can be modified easily.

**What is Provided**

| File Name      | Function |
| ----------- | ----------- |
| `1kg_exchange_test.bash`     | Automates execution of test of exchangeability for all 26 populations of 1KG       |
| `test_exchangeability.R`   | Executes test of exchangeability from terminal, given three arguments: directory to 1KG population files (BED/BIM/FAM), directory to PLINK2, resampling number for permutation test      |
| `pop_id_list.txt` | List of 1KG populations for automating execution of exchangeability test across all 26 populations | 

Additionally, example output log files are provided in the subdirectory `examples`.  

**What isn't Provided**

- PLINK2. If you don't have PLINK2 already installed, you have to [download](https://www.cog-genomics.org/plink/2.0/) PLINK2 and place it in an accessible directory.
- The 1000 Genomes population-specific datasets (26 of them), which are large and inconvenient for provision by this repo. We are currently exploring other options to make this publicly available. 

**How to Run**  

1. Download this directory. 
2. Download the 1000G datasets and run PLINK commands as described [here](https://alanaw1.github.io/flintyR/articles/extras.html#running-our-test-from-terminal-1) (*requires knowledge of PLINK*).
3. Place the 1000G datasets in your favourite directory. Modify line 6 of `1kg_exchange_test.bash` to point to the directory.
4. Run the following command in the terminal (*recommended: use a job scheduler, or a terminal multiplexer like [screen](https://blog.thibaut-rousseau.com/2015/12/04/screen-terminal-multiplexer.html)*)

    ```
    bash 1kg_exchange_test.bash
    ```

    Just a heads up, this script takes about 73 minutes to run to completion on a standard MacBook Pro. You can see this by printing the terminal output after the second line below is entered.

    ```
    grep "Time difference of" output_log_manhattan.txt | cut -d " " -f4 > time.txt

    ( echo 0 3k ; sed 's/$/ +/' time.txt ; echo p) | dc # 4393.699899 seconds

    rm time.txt
    ``` 

5. There will be an output log file (named `output_log.txt`) that should look like `examples/R_output_log_manhattan.txt`. 

## For Python Users

This directory contains scripts and data files used to run our test of exchangeability from the terminal. As described in our Examples page, this demonstration involves the use of the [1000 Genomes Phase 3 dataset](https://www.cog-genomics.org/plink/2.0/resources). 

We encourage interested users to repurpose our scripts for their own analyses. In particular, the Python script, written to be executable from the terminal, can be modified easily.

**What is Provided**

| File Name      | Function |
| ----------- | ----------- |
| `1kg_exchange_test_python.bash`     | Automates execution of test of exchangeability for all 26 populations of 1KG       |
| `test_exchangeability.py`   | Executes test of exchangeability from terminal, given three arguments: directory to 1KG population files (BED/BIM/FAM), directory to PLINK2, resampling number for permutation test      |
| `pop_id_list.txt` | List of 1KG populations for automating execution of exchangeability test across all 26 populations | 

Additionally, an example output log file is provided in the subdirectory `examples`.  

**How to Run**  

1. Download this directory. 
2. Download the 1000G datasets and run PLINK commands as described [here](https://alanaw1.github.io/flintyR/articles/extras.html#running-our-test-from-terminal-1) (*requires knowledge of PLINK*).
3. Place the 1000G datasets in your favourite directory. Modify line 6 of `1kg_exchange_test_python.bash` to point to the directory.
4. Activate the virtual environment `flinty`.  
5. Run the following command in the terminal (*recommended: use a job scheduler, or a terminal multiplexer like [screen](https://blog.thibaut-rousseau.com/2015/12/04/screen-terminal-multiplexer.html)*)

    ```
    bash 1kg_exchange_test_python.bash
    ```
6. There will be an output log file (named `python_output_log.txt`) that should look like `examples/python_output_log_manhattan.txt`. 