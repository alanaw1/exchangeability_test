# Wrangling World Values Survey Data

This directory contains the scripts and data files used to generate the WVS data we analysed in our work. The subdirectory `MM_2020` contains allele dimensions for WVS questions and WVS variables included in their study ([Muthukrishna et al., 2020](https://journals.sagepub.com/doi/full/10.1177/0956797620916782)). 

**How to Run**  

1. Download this directory. 
2. Go to [WVS Wave 6](https://www.worldvaluessurvey.org/WVSDocumentationWV6.jsp) and fill in the registration form to download the data file. (We used `WV6 Data R v20180912` for our analysis.) 
3. Place the downloaded WVS data file (`F00007762-WV6_Data_R_v20180912`) in the same directory as the scripts.
4. Run the following command in the terminal:

    ```
    Rscript gen_WVS_df.R F00007762-WV6_Data_R_v20180912.Rds MM_2020/included-variables.csv
    ```

5. There will be two output files, one named `all_MM_num_filtered_recode.RData`. This is the dataframe we load and analyse in our vignette, [Analysis of World Values Survey Data](https://alanaw1.github.io/flintyR/articles/wvs.html).