## ---------------------------
##
## Script name: DMA
##
## Purpose of script: Differential Metabolomics Analysis
##
## Author: Dimitrios Prymidis and Christina Schmidt
##
## Date Created: 2022-10-28
##
## Copyright (c) Dimitrios Prymidis and Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

#' This script allows you to perform differential metabolite analysis to obtain a Log2FC, pval, padj and tval comparing two conditions.
#'
#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. includes experimental Experimental_design and outlier column.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates".
#' @param Conditon1 Input needs to contain a column named "Condition" including the Condition1 that will be compared to Condition2, e.g. "KO".
#' @param Conditon2 Input needs to contain a column named "Condition" including Condition2 that is compared to Condition1, e.g. "WT".
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value (t.test or wilcox.test) \strong{"t-test"}
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{"fdr"}
#' @param OutputName String which is added to the output files of the DMA.
#' @param Input_Pathways \emph{Optional: } DF which contains the pathway information for each metabolite. \strong{NULL}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used \strong{FALSE}
#' @param plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{TRUE}
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf, jpeg, tiff, bmp. \strong{Default = svg}
#'
#' @keywords Differential Metabolite Analysis, Multiple Hypothesis testing, Normality testing
#' @export


########################################################
### ### ### Differential Metabolite Analysis ### ### ###
########################################################

DMA <-function(Input_data,
               Experimental_design,
               Condition1,
               Condition2,
               STAT_pval ="t.test",
               STAT_padj="fdr",
               Input_Pathways = NULL,
               OutputName='',
               CoRe=FALSE,
               Plot = TRUE,
               Save_as = "svg"){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "gtools", "EnhancedVolcano")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
    }
  suppressMessages(library(tidyverse))

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions

  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
    } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
      stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
    } else{
      Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
      if((any(Test_num) ==  FALSE) ==  TRUE){
        stop("Input_data needs to be of class numeric")
      } else{
        Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
        if(nrow(Test_match) ==  0){
          stop("row.names Input_data need to match row.names Experimental_design.")
        } else{
           Input_data <- Input_data
        }
      }
    }

      C1 <- Input_data %>%
        filter(Experimental_design$Conditions %in% Condition1) %>%
        select_if(is.numeric)#only keep numeric columns with metabolite values
      C2 <- Input_data %>%
        filter(Experimental_design$Conditions %in% Condition2) %>%
        select_if(is.numeric)

      if(nrow(C1)==1){
        stop("There is only one sample available for ", Condition1, ", so no statistical test can be performed.")
      } else if(nrow(C2)==1){
        stop("There is only one sample available for ", Condition2, ", so no statistical test can be performed.")
      }else if(nrow(C1)==0){
        stop("There is no sample available for ", Condition1, ".")
      }else if(nrow(C2)==0){
        stop("There is no sample available for ", Condition2, ".")
      }

  #2. General parameters
  STAT_pval_options <- c("t.test", "wilcox.test","chisq.test", "cor.test")
  if(STAT_pval %in% STAT_pval_options == FALSE){
    stop("Check input. The selected STAT_pval option for Hypothesis testing is not valid. Please select one of the folowwing: ",paste(STAT_pval_options,collapse = ", "),"." )
  }
  STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(STAT_padj %in% STAT_padj_options == FALSE){
    stop("Check input. The selected STAT_padj option for multiple Hypothesis testing correction is not valid. Please select one of the folowwing: ",paste(STAT_padj_options,collapse = ", "),"." )
  }
  if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for analysis of Consuption/Release experiment or =FALSE if not.")
  }
  if(is.logical(Plot) == FALSE){
    stop("Check input. The plot value should be either =TRUE if a Volcano plot presenting the DMA results is to be exported or =FALSE if not.")
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }

  #3. Input_Pathways
  if(is.null(Input_Pathways) == FALSE){
    if('Metabolite' %in% colnames(Input_Pathways) == FALSE){
      warning("The provided file Input_Pathways must have 2 columns named: `Metabolite` and `Pathway`.")
    }else if('Pathway' %in% colnames(Input_Pathways) == FALSE){
     warning("The provided file Input_Pathways must have 2 columns named: `Metabolite` and `Pathway`.")
    }
  }

  ## ------------ Check Missingness ------------- ##
  # If missing value imputation has not been performed the input data will most likely contain NA or 0 values for some metabolites, which will lead to Log2FC = NA.
  # Here we will check how many metabolites this affects in C1 and C2, and weather all replicates of a metabolite are affected.
  C1_Miss <- replace(C1, C1==0, NA)
  C1_Miss <- C1_Miss[, (colSums(is.na(C1_Miss)) > 0), drop = FALSE]

  C2_Miss <- replace(C2, C2==0, NA)
  C2_Miss <- C2_Miss[, (colSums(is.na(C2_Miss)) > 0), drop = FALSE]

  if((ncol(C1_Miss)>0 & ncol(C2_Miss)==0)){
    message("In `Condition1` ",paste0(toString(Condition1)), ", NA/0 values exist in ", ncol(C1_Miss), " Metabolite(s): ", paste0(colnames(C1_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(C1_Miss)
    } else if(ncol(C1_Miss)==0 & ncol(C2_Miss)>0){
    message("In `Condition2` ",paste0(toString(Condition2)), ", NA/0 values exist in ", ncol(C2_Miss), " Metabolite(s): ", paste0(colnames(C2_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(C2_Miss)
    } else if(ncol(C1_Miss)>0 & ncol(C2_Miss)>0){
    message("In `Condition1` ",paste0(toString(Condition1)), ", NA/0 values exist in ", ncol(C1_Miss), " Metabolite(s): ", paste0(colnames(C1_Miss), collapse = ", "), " and in `Condition2`",paste0(toString(Condition2)), " ",ncol(C2_Miss), " Metabolite(s): ", paste0(colnames(C2_Miss), collapse = ", "),". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- c(colnames(C1_Miss), colnames(C2_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
    } else{
    message("There are no NA/0 values")
    Metabolites_Miss <- c(colnames(C1_Miss), colnames(C2_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
    }

  ## ------------ Check data normality and statistical test chosen ----------- ##
  # Before Hypothesis testing, we have to decide whether to use a parametric or a non parametric test. we can test the data normality using the Shapiro test.

  # 1. Load the data and perform the shapiro.test on each metabolite:
  Input_data_NA <- replace(Input_data, Input_data==0, NA)#Shapiro test ignores NAs!

  shaptest <-   Input_data_NA %>% # Select data
    filter(Experimental_design$Conditions %in% Condition1 | Experimental_design$Conditions %in% Condition2)%>%
    select_if(is.numeric)

  temp<- as.vector(sapply(shaptest, function(x) var(x)) == 0)#  we have to remove features with zero variance if there are any.
  shaptest <- shaptest[,!temp]

  shaptestres <- as.data.frame(sapply(shaptest, function(x) shapiro.test(x))) # do the test for each metabolite
  shaptestres <- as.data.frame(t(shaptestres))

  # 2. Give feedback to the user if the chosen test fits the data distribution. The data are normal if the p-value of the shapiro.test > 0.05.
  Norm <- format((round(sum(shaptestres$p.value > 0.05)/dim(shaptest)[2],4))*100, nsmall = 2) # Percentage of normally distributed metabolites across samples
  NotNorm <- format((round(sum(shaptestres$p.value < 0.05)/dim(shaptest)[2],4))*100, nsmall = 2) # Percentage of not-normally distributed metabolites across samples
  message(Norm, " % of the metabolites follow a normal distribution and ", NotNorm, " % of the metabolites are not-normally distributed according to the shapiro test. `shapiro.test` ignores missing values in the calculation.")

  if (Norm > 50 & STAT_pval =="wilcox.test"){
      warning(Norm, " % of the metabolites follow a normal distribution but 'wilcox.test' for non parametric Hypothesis testing was chosen. Please consider selecting a parametric test (t.test) instead.")
    }else if(NotNorm > 50 & STAT_pval =="t.test"){
    message(NotNorm, " % of the metabolites follow a not-normal distribution but 't.test' for parametric Hypothesis testing was chosen. Please consider selecting a non-parametric test (wilcox.test) instead.")
    }else if((STAT_pval =="wilcox-test" & nrow(C1)<5)|(STAT_pval =="wilcox-test" & nrow(C2)<5)){# check number of samples for wilcoxons test
      warning("Number of samples measured per condition is <5 in at least one of the two conditions, which is small for using wilcox.test. Consider using another test.")
    }

  ## ------------ Create Results output folder ----------- ##
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name) # Make Results folder
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  Results_folder_DMA_folder <- file.path(Results_folder,"DMA") # Make DMA results folder
  if (!dir.exists(Results_folder_DMA_folder)) {dir.create(Results_folder_DMA_folder)}
  Results_folder_Conditions <- file.path(Results_folder_DMA_folder,paste0(toString(Condition1),"_vs_",toString(Condition2))) # Make comparison folder
  if (!dir.exists(Results_folder_Conditions)) {dir.create(Results_folder_Conditions)}

  ################################################################################################################################################################################################
  ############### Calculate Log2FC, pval, padj, tval ###############

  ## ------------  Calculate Log2FC ----------- ##
  # For C1_Mean and C2_Mean use 0 to obtain values, leading to Log2FC=NA if mean = 0 (If one value is NA, the mean will be NA even though all other values are available.)
  C1_Zero <- C1
  C1_Zero[is.na(C1_Zero)] <- 0
  Mean_C1 <- C1_Zero %>%
    summarise_all("mean")

  C2_Zero <- C2
  C2_Zero[is.na(C2_Zero)] <- 0
  Mean_C2 <- C2_Zero %>%
    summarise_all("mean")

  if(CoRe==TRUE){#Calculate Log2FC by taking into account the distance between the means:
    #CoRe values can be negative and positive, which can lead to problems when calculating the Log2FC and hence the Log2FC(A versus B)=(log2(A+x)-log2(B+x)) for A or B being negative, with x being a constant that is adapted to the size range of the respective metabolite.
    Mean_C1_t <- as.data.frame(t(Mean_C1))%>%
      rownames_to_column("Metabolite")
    Mean_C2_t <- as.data.frame(t(Mean_C2))%>%
      rownames_to_column("Metabolite")
    Mean_Merge <-merge(Mean_C1_t, Mean_C2_t, by="Metabolite", all=TRUE)%>%
      rename("C1"=2,
             "C2"=3)
    Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss#Column to enable the check if mean values of 0 are due to missing values (NA/0) and not by coincidence

    Mean_Merge <- Mean_Merge%>%#Now we can adapt the values to take into account the distance
      mutate(C1_Adapted = case_when(C1 < 0 & C2 > 0 ~ paste(C1*-1),#Here we have a negative value and transform it into a positive value
                                    C1 > 0 & C2 < 0 ~ paste(C1+(C2*-2)),#Here we have a positive value, but the other condition has a negative value that is transformed to a positive value, hence we add the distance
                                    C1 > 0 & C2 > 0 ~ paste(C1),#Both values are positive, so no action needed
                                    C1 < 0 & C2 < 0 ~ paste(C1),#Both values are negative, so no action needed
                                    C2 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                    C1 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                    C2 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                    C1 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                    TRUE ~ 'NA'))%>%
      mutate(C2_Adapted = case_when(C2 < 0 & C1 > 0 ~ paste(C2*-1),#Here we have a negative value and transform it into a positive value
                                    C2 > 0 & C1 < 0 ~ paste(C2+(C1*-2)),#Here we have a positive value, but the other condition has a negative value that is transformed to a positive value, hence we add the distance
                                    C2 > 0 & C1 > 0 ~ paste(C2),#Both values are positive, so no action needed
                                    C2 < 0 & C1 < 0 ~ paste(C2),#Both values are negative, so no action needed
                                    C1 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                    C2 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                    C1 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                    C2 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                    TRUE ~ 'NA'))%>%
      mutate(C1_Adapted = as.numeric(C1_Adapted))%>%
      mutate(C2_Adapted = as.numeric(C2_Adapted))

    if(any((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))==TRUE){
      X <- Mean_Merge%>%
        subset((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))
      message("We added +1 to the mean value of metabolite(s) ", paste0(X$Metabolite, collapse = ", "), ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
    }

    #Calculate the Log2FC
    Mean_Merge$FC_C1vC2 <- Mean_Merge$C1_Adapted/Mean_Merge$C2_Adapted #FoldChange
    Mean_Merge$Log2FC <- gtools::foldchange2logratio(Mean_Merge$FC_C1vC2, base=2)
    Log2FC_C1vC2 <-Mean_Merge[,c(1,8)]

    }else if(CoRe==FALSE){
      #Mean values could be 0, which can not be used to calculate a Log2FC and hence the Log2FC(A versus B)=(log2(A+x)-log2(B+x)) for A and/or B being 0, with x being set to 1
      Mean_C1_t <- as.data.frame(t(Mean_C1))%>%
        rownames_to_column("Metabolite")
      Mean_C2_t <- as.data.frame(t(Mean_C2))%>%
        rownames_to_column("Metabolite")
      Mean_Merge <-merge(Mean_C1_t, Mean_C2_t, by="Metabolite", all=TRUE)%>%
        rename("C1"=2,
               "C2"=3)
      Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss#Column to enable the check if mean values of 0 are due to missing values (NA/0) and not by coincidence

      Mean_Merge <- Mean_Merge%>%#Now we can adapt the values to take into account the distance
        mutate(C1_Adapted = case_when(C2 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C1 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C2 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      C1 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      TRUE ~ paste(C1)))%>%
        mutate(C2_Adapted = case_when(C1 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C2 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                      C1 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      C2 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                      TRUE ~ paste(C2)))%>%
        mutate(C1_Adapted = as.numeric(C1_Adapted))%>%
        mutate(C2_Adapted = as.numeric(C2_Adapted))

      if(any((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))==TRUE){
        X <- Mean_Merge%>%
          subset((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))
        message("We added +1 to the mean value of metabolite(s) ", paste0(X$Metabolite, collapse = ", "), ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
      }

      #Calculate the Log2FC
      Mean_Merge$FC_C1vC2 <- Mean_Merge$C1_Adapted/Mean_Merge$C2_Adapted #FoldChange
      Mean_Merge$Log2FC <- gtools::foldchange2logratio(Mean_Merge$FC_C1vC2, base=2)
      Log2FC_C1vC2 <-Mean_Merge[,c(1,8)]

      }else{
        stop("Please choose CoRe= TRUE or CoRe=FALSE.")
        }

  ## ------------ Perform Hypothesis testing ----------- ##
  # For C1 and C2 we use 0, since otherwise we can not perform the statistical testing.
  C1[is.na(C1)] <- 0
  C2[is.na(C2)] <- 0

  #### 1. p.value
  T_C1vC2 <-mapply(STAT_pval, x= as.data.frame(C2), y = as.data.frame(C1), SIMPLIFY = F)

  VecPVAL_C1vC2 <- c()
  for(i in 1:length(T_C1vC2)){
    p_value <- unlist(T_C1vC2[[i]][3])
    VecPVAL_C1vC2[i] <- p_value
  }
  Metabolite <- colnames(C2)
  PVal_C1vC2 <- data.frame(Metabolite, VecPVAL_C1vC2)

  #we set p.val= NA, for metabolites that had 1 or more replicates with NA/0 values and remove them prior to p-value adjustment
  PVal_C1vC2$`NA/0` <- PVal_C1vC2$Metabolite %in% Metabolites_Miss
  PVal_C1vC2 <-PVal_C1vC2%>%
    mutate(p.val = case_when(`NA/0`== TRUE ~ NA,
                             TRUE ~ paste(VecPVAL_C1vC2)))
  PVal_C1vC2$p.val = as.numeric(as.character(PVal_C1vC2$p.val))

  #### 2. p. adjusted
  #Split data for p.value adjustment to exclude NA
  PVal_NA <- PVal_C1vC2[is.na(PVal_C1vC2$p.val), c(1,4)]
  PVal_C1vC2 <-PVal_C1vC2[!is.na(PVal_C1vC2$p.val), c(1,4)]

  #perform adjustment
  VecPADJ_C1vC2 <- p.adjust((PVal_C1vC2[,2]),method = STAT_padj, n = length((PVal_C1vC2[,2]))) #p-adjusted
  Metabolite <- PVal_C1vC2[,1]
  PADJ_C1vC2 <- data.frame(Metabolite, VecPADJ_C1vC2)
  STAT_C1vC2 <- merge(PVal_C1vC2,PADJ_C1vC2, by="Metabolite")
  STAT_C1vC2 <- merge(Log2FC_C1vC2,STAT_C1vC2, by="Metabolite")
  names(STAT_C1vC2)[names(STAT_C1vC2) == "VecPADJ_C1vC2"] <- "p.adj"

  #### 3. t.value
  STAT_C1vC2$t.val <- qnorm((1 - STAT_C1vC2$p.val / 2)) * sign(STAT_C1vC2$Log2FC) # calculate and add t-value
  STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value

  #Add Metabolites that have p.val=NA back into the DF for completeness.
  if(nrow(PVal_NA)>0){
  PVal_NA <- merge(Log2FC_C1vC2,PVal_NA, by="Metabolite", all.y=TRUE)
  PVal_NA$p.adj <- NA
  PVal_NA$t.val <- NA
  STAT_C1vC2 <- rbind(STAT_C1vC2, PVal_NA)
  }


  ################################################################################################################################################################################################

  ## ------------ Add information on groups to DMA results----------- ##
  if(CoRe==TRUE){#add consumption release info to the result
    temp1 <- Mean_C1
    temp2 <- Mean_C2
    #Add Info:
    CoRe_info <- rbind(temp1, temp2,rep(0,length(temp1)))
    for (i in 1:length(temp1)){
      if (temp1[i]>0 & temp2[i]>0){
        CoRe_info[3,i] <- "Released"
      }else if (temp1[i]<0 & temp2[i]<0){
        CoRe_info[3,i] <- "Consumed"
      }else if(temp1[i]>0 & temp2[i]<0){
        CoRe_info[3,i] <- paste("Released in" ,Condition1 , "and Consumed",Condition2 , sep=" ")
      } else if(temp1[i]<0 & temp2[i]>0){
        CoRe_info[3,i] <- paste("Consumed in" ,Condition1 , " and Released",Condition2 , sep=" ")
      }else{
        CoRe_info[3,i] <- "No Change"
      }
    }

    CoRe_info <- t(CoRe_info) %>% as.data.frame()
    CoRe_info <- rownames_to_column(CoRe_info, "Metabolite")
    names(CoRe_info)[2] <- paste("Mean", "CoRe", Condition1, sep="_")
    names(CoRe_info)[3] <- paste("Mean", "CoRe", Condition2, sep="_")
    names(CoRe_info)[4] <- "CoRe"
    STAT_C1vC2 <- merge(STAT_C1vC2,CoRe_info,by= "Metabolite")
    STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value
  }

  ## ------------ Add pathway information to DMA results ----------- ##
  if(is.null(Input_Pathways)!=TRUE & 'Metabolite' %in% colnames(Input_Pathways) & 'Pathway' %in% colnames(Input_Pathways)){
      STAT_C1vC2$Metabolite %in% Input_Pathways$Metabolite
      STAT_C1vC2<- merge(STAT_C1vC2,Input_Pathways,by="Metabolite", all.x=T)
      STAT_C1vC2$Pathway[  is.na(STAT_C1vC2$Pathway)] <- "unknown"    # Merge the pathways to DMA result. All non matching metabolites get NA that are changed into "unknown".
    }else{
      warning("No column Pathway was added to the output. The pathway data must have 2 columns named: Metabolite and Pathway")
    }
  DMA_Output <- STAT_C1vC2


  xlsDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(Condition1),"_vs_",toString(Condition2),"_", OutputName, ".xlsx"))   # Save the DMA results table
  writexl::write_xlsx(DMA_Output,xlsDMA, col_names = TRUE) # save the DMA result DF

  if(Plot == TRUE){ # Make a simple Volcano plot
    VolcanoPlot<- EnhancedVolcano::EnhancedVolcano(DMA_Output,
                                   lab = DMA_Output$Metabolite,#Metabolite name
                                   x = "Log2FC",#Log2FC
                                   y = "p.adj",#p-value or q-value
                                   xlab = bquote(~Log[2]~ FC),
                                   ylab = bquote(~-Log[10]~p.adj),#(~-Log[10]~adjusted~italic(P))
                                   pCutoff = 0.05,
                                   FCcutoff = 0.5,#Cut off Log2FC, automatically 2
                                   pointSize = 3,
                                   labSize = 2,
                                   titleLabSize = 16,
                                   # colCustom = c("black", "grey", "grey", "red"),
                                   colAlpha = 0.7,
                                   title= paste0(toString(Condition1)," versus ",toString(Condition2)),
                                   subtitle = bquote(italic("Differential Metabolite Analysis")),
                                   caption = paste0("total = ", nrow(DMA_Output), " Metabolites"),
                                   xlim =  c(min(DMA_Output$Log2FC[is.finite(DMA_Output$Log2FC )])-0.2,max(DMA_Output$Log2FC[is.finite(DMA_Output$Log2FC )])+0.2),
                                   ylim = c(0,(ceiling(-log10(Reduce(min,DMA_Output$p.adj))))),
                                   cutoffLineType = "dashed",
                                   cutoffLineCol = "black",
                                   cutoffLineWidth = 1,
                                   legendLabels=c('No changes',paste(0.5,"< |Log2FC|"),paste("p.adj <",0.05) , paste('p.adj<',0.05,' &',0.5,"< |Log2FC|")),
                                   legendPosition = 'right',
                                   legendLabSize = 8,
                                   legendIconSize =4,
                                   gridlines.major = FALSE,
                                   gridlines.minor = FALSE,
                                   drawConnectors = FALSE)
    OutputPlotName = paste0(OutputName,"_padj_",0.05,"Log2FC_",0.5)

    volcanoDMA <- file.path(Results_folder_Conditions,paste0( "Volcano_Plot_",toString(Condition1),"-versus-",toString(Condition2),"_", OutputPlotName,".",Save_as))
    ggsave(volcanoDMA,plot=VolcanoPlot, width=10, height=8) # save the voplcano plot

    plot(VolcanoPlot)
  }
  assign(paste0("DMA_",toString(Condition1),"_vs_",toString(Condition2)), DMA_Output, envir=.GlobalEnv)
}




