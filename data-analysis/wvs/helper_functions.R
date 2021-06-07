#' Get Raw Responses
#'
#' Extracts labels attribute from a numerically encoded feature of the WVS tibble. 
#' Used to convert the original raw WVS tibble, with encoded responses, into raw 
#' responses to actual survey questions. 
#'
#' Dependencies: None
#' @param x A vector of class labelled 
#' @return A vector of same length as x, whose entries are the raw responses to WVS
getRawResponses <- function(x) {
  # Extract mapping from raw labels to encoded
  codes <- attr(x, "labels")
  
  if (is.null(codes)) {
    print(paste("No labels with", attr(x, "label")))
    return(x)
  }
  
  df <- data.frame(y = as.numeric(x), 
                   stringsAsFactors = FALSE)
  
  replacement <- data.frame(
    y = as.numeric(codes),
    newform = names(codes),
    stringsAsFactors = FALSE
  )
  
  # Join encoded to mapped raw response for each individual 
  df2 <- df %>% 
    left_join(replacement, by="y") %>%
    mutate(newform = ifelse(is.na(newform), 
                            as.character(x), newform))
  
  # Check that extracted transform has same length as original 
  if(length(x) == length(df2$newform) &
     sum(is.na(x)) == sum(is.na(df2$newform))) {
    return(df2$newform)
  } else {
    print(paste("Problem with", attr(x, "label")))
    return(x)
  }
}

#' Keep Features Used in Cultural Distances App (with Recoding)
#'
#' Removes features in WVS tibble that do not appear in MM's codebook. Additionally,
#' automatically recodes the individual responses according to the codebook.
#'
#' Dependencies: getRawResponses
#' @param wvs_obj The WVS data loaded from the RDS file
#' @param culdist_metadata The codebook loaded from the csv file of Muthukrishna et al.
#' @return WVS data including only features appearing in codebook, with recoding
keepCulDistFeatures <- function(wvs_obj, culdist_metadata) {
  # Obtain question labels
  cat(">>> Obtaining Question Labels...\n")
  q_labels <- wvs_obj %>% map_chr(~attributes(.)$label) 
  my_df <- data.frame(question = q_labels,
                      codebookID = names(q_labels),
                      row.names = NULL) 
  
  # Find questions with HEX-coded symbols
  cat(">>> Finding Questions with HEX-coded Symbols...\n")
  blacklist <- grep("\x92|\xb4|\xf3|\x93|\x94", culdist_metadata$WVS.label, perl = TRUE)
  
  # Convert cultural distance codebook questions to lower casing
  WVS.label <- c()
  for (i in 1:length(culdist_metadata$WVS.label)) {
    if (i %in% blacklist) {
      #cat("String ", i, " in blacklist...\n", sep = "")
      WVS.label <- c(WVS.label, 
                     culdist_metadata$WVS.label[i])
    } else {
      #cat("String ", i, " not in blacklist...\n", sep = "")
      WVS.label <- c(WVS.label, 
                     tolower(culdist_metadata$WVS.label[i]))
    }
  }
  culdist_metadata$WVS.label <- WVS.label
  
  # Convert WVS questions to lower casing
  target_q <- tolower(my_df$question)
  
  # Check each question for inclusion in MM codebook, and
  # include the WVS variable 
  cat(">>> Finding Matching WVS Questions in MM Codebook...\n")
  target_q_local <- target_q # local copy to be modified and used
  match_counter <- 0
  for (i in 1:length(target_q)) {
    for (j in 1:length(culdist_metadata$WVS.label)) {
      if (grepl(culdist_metadata$WVS.label[j], target_q[i], perl = TRUE)) {
        cat("Match found --- ", culdist_metadata$WVS.label[j], 
            " for i = ", i, " and j = ", j, "\n",
            sep = "")
        match_counter <- match_counter + 1 
        target_q_local[i] <- culdist_metadata$WVS.label[j]
        break
      }
    }
  }
  
  wvs_ID_and_inclusion <- sapply(target_q_local, function(x) {
    if (nrow(culdist_metadata %>% filter(WVS.label == as.character(x))) > 0) {
      output <- culdist_metadata %>% filter(WVS.label == as.character(x))
      return(c(as.character(output$WVS.variable), 
               as.character(output$Included.in.cultural.distance)))
    } else {
      return(c(NA,NA))
    }  
  })
  wvs_ID_and_inclusion <- t(wvs_ID_and_inclusion)
  colnames(wvs_ID_and_inclusion) <- c("WVS.variable", "Included") 
  
  # Combine dataframes to obtain (Question, Codebook ID, WVS variable, Inclusion) df
  QCWI_df <- cbind(my_df, wvs_ID_and_inclusion)
  rownames(QCWI_df) <- NULL
  
  # Use getRawResponses to convert tibble of 
  # encoded responses to raw responses
  wv <- wvs_obj %>% 
    mutate_all(getRawResponses) 
  
  # Restrict WVS tibble to only questions included in MM codebook
  # Resulting tibble should have features like country, wave, regions
  # appearing as themselves, but the question responses are encoded numerics
  q_included <- QCWI_df %>% filter(Included == "Yes")
  wv_num_filtered <- cbind(wv[,1:4],(wvs_obj[,-(1:4)])[,q_included$codebookID]) 
  wv_num_filtered_ <- wv_num_filtered %>% 
    `colnames<-`(c(q_labels[c("V1","V2","V2A","V3")], 
                   as.character(q_included$WVS.variable)))
  wv_num_filtered_[,4:dim(wv_num_filtered_)[2]] <- 
    apply(wv_num_filtered_[,4:dim(wv_num_filtered_)[2], 
                           drop = FALSE], 2, as.numeric)
  
  # Convert using the recoding scheme provided by MM codebook
  cat(">>> Recoding Individual Responses Using MM Codebook...\n")
  wv_num_filtered_ <- wv_num_filtered_ %>% 
    mutate(A001 = dplyr::case_when(A001 %in% 1:2 ~ 1, 
                                   A001 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(A002 = dplyr::case_when(A002 %in% 1:2 ~ 1, 
                                   A002 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(A003 = dplyr::case_when(A003 %in% 1:2 ~ 1, 
                                   A003 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(A004 = dplyr::case_when(A004 %in% 1:2 ~ 1, 
                                   A004 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(A005 = dplyr::case_when(A005 %in% 1:2 ~ 1, 
                                   A005 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(A006 = dplyr::case_when(A006 %in% 1:2 ~ 1, 
                                   A006 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(A029 = sapply(A029, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A030 = sapply(A030, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A032 = sapply(A032, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A034 = sapply(A034, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A035 = sapply(A035, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A038 = sapply(A038, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A040 = sapply(A040, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A041 = sapply(A041, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A042 = sapply(A042, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(A165 = sapply(A165, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(D059 = dplyr::case_when(D059 %in% 1:2 ~ 1, 
                                   D059 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E001 = sapply(E001, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(E003 = sapply(E003, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(E005 = sapply(E005, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(E015 = sapply(E015, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(E018 = sapply(E018, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(B008 = sapply(B008, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(E023 = dplyr::case_when(E023 %in% 1:2 ~ 1, 
                                   E023 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E025 = sapply(E025, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(E026 = sapply(E026, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(E033 = dplyr::case_when(E033 %in% 1:4 ~ 1, 
                                   E033 %in% 5:6 ~ 2, 
                                   E033 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E035 = dplyr::case_when(E035 %in% 1:4 ~ 1, 
                                   E035 %in% 5:6 ~ 2, 
                                   E035 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E036 = dplyr::case_when(E036 %in% 1:4 ~ 1, 
                                   E036 %in% 5:6 ~ 2, 
                                   E036 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E037 = dplyr::case_when(E037 %in% 1:4 ~ 1, 
                                   E037 %in% 5:6 ~ 2, 
                                   E037 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E039 = dplyr::case_when(E039 %in% 1:4 ~ 1, 
                                   E039 %in% 5:6 ~ 2, 
                                   E039 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E114 = dplyr::case_when(E114 %in% 1:2 ~ 1, 
                                   E114 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E116 = dplyr::case_when(E116 %in% 1:2 ~ 1, 
                                   E116 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(E117 = dplyr::case_when(E117 %in% 1:2 ~ 1, 
                                   E117 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F001 = sapply(F001, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(F028 = dplyr::case_when(F028 %in% 1:2 ~ 1, 
                                   F028 %in% 3:5 ~ 2, 
                                   F028 %in% 6:7 ~ 3,
                                   TRUE ~ NA_real_)) %>%
    mutate(F034 = sapply(F034, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(F050 = sapply(F050, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(F053 = sapply(F053, function(x) ifelse(x < 0, NA, x))) %>%
    mutate(F063 = dplyr::case_when(F063 %in% 1:4 ~ 1, 
                                   F063 %in% 5:6 ~ 2, 
                                   F063 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F114 = dplyr::case_when(F114 %in% 1:4 ~ 1, 
                                   F114 %in% 5:6 ~ 2, 
                                   F114 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F115 = dplyr::case_when(F115 %in% 1:4 ~ 1, 
                                   F115 %in% 5:6 ~ 2, 
                                   F115 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F116 = dplyr::case_when(F116 %in% 1:4 ~ 1, 
                                   F116 %in% 5:6 ~ 2, 
                                   F116 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F117 = dplyr::case_when(F117 %in% 1:4 ~ 1, 
                                   F117 %in% 5:6 ~ 2, 
                                   F117 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F118 = dplyr::case_when(F118 %in% 1:4 ~ 1, 
                                   F118 %in% 5:6 ~ 2, 
                                   F118 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F119 = dplyr::case_when(F119 %in% 1:4 ~ 1, 
                                   F119 %in% 5:6 ~ 2, 
                                   F119 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F120 = dplyr::case_when(F120 %in% 1:4 ~ 1, 
                                   F120 %in% 5:6 ~ 2, 
                                   F120 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F121 = dplyr::case_when(F121 %in% 1:4 ~ 1, 
                                   F121 %in% 5:6 ~ 2, 
                                   F121 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F123 = dplyr::case_when(F123 %in% 1:4 ~ 1, 
                                   F123 %in% 5:6 ~ 2, 
                                   F123 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(F122 = dplyr::case_when(F122 %in% 1:4 ~ 1, 
                                   F122 %in% 5:6 ~ 2, 
                                   F122 %in% 7:10 ~ 3, 
                                   TRUE ~ NA_real_)) %>%
    mutate(G006 = dplyr::case_when(G006 %in% 1:2 ~ 1, 
                                   G006 %in% 3:4 ~ 2, 
                                   TRUE ~ NA_real_))
  
  colnames(wv_num_filtered_) <- q_labels[c("V1","V2","V2A","V3",
                                           as.character(q_included$codebookID))]
  
  # Return WVS df with only MM codebook questions included,
  # and with recoding performed according to MM codebook
  return(wv_num_filtered_)
}
