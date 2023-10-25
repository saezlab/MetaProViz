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

#' This script allows you to perform differential metabolite analysis to obtain a Log2FC, pval, padj and tval comparing two or multiple conditions.
#'
#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Input_SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param Input_SettingsInfo \emph{Optional: } Named vector including the information about the conditions column c(conditions="ColumnName_Plot_SettingsFile"). Can additionally pass information on numerator or denominator c(numerator = "ColumnName_Plot_SettingsFile", denominator  = "ColumnName_Plot_SettingsFile") for specifying which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all). Using =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator and are presented in the right side on the Volcano plot (For CoRe the Log2Distance). \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test", "cor.test" or lmFit (=limma), for one-vs-all or all-vs-all comparison choose aov (=annova), kruskal.test or lmFit (=limma) \strong{Default = "lmFit"}
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param OutputName String which is added to the output files of the DMA.
#' @param Input_SettingsFile_Metab \emph{Optional: } DF which contains the metadata information , i.e. pathway information, retention time,..., for each metabolite. \strong{Default = NULL}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used \strong{Default = FALSE}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{Default = "csv"}
#' @param plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @keywords Differential Metabolite Analysis, Multiple Hypothesis testing, Normality testing
#' @export


########################################################
### ### ### Differential Metabolite Analysis ### ### ###
########################################################

DMA <-function(Input_data,
               Input_SettingsFile_Sample,
               Input_SettingsInfo = c(conditions="Conditions", numerator = NULL, denominator  = NULL),
               STAT_pval ="lmFit",
               STAT_padj="fdr",
               Input_SettingsFile_Metab = NULL,
               OutputName='',
               CoRe=FALSE,
               Save_as_Plot = "svg",
               Save_as_Results = "csv",
               Plot = TRUE,
               Folder_Name = NULL
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "gtools", "EnhancedVolcano")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install(new.packages)
  }
  suppressMessages(library(tidyverse))

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(class(Input_data) != "data.frame"){
    stop("Input_data should be a data.frame. It's currently a ", paste(class(Input_data), ".",sep = ""))
  }
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Input_SettingsFile_Sample, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE)
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Input_SettingsFile_Sample.")
      } else{
        Input_data <- Input_data
      }
    }
  }

  #2.  Input_SettingsFile_Metab
  if(is.null(Input_SettingsFile_Metab) == FALSE){
    if('Metabolite' %in% colnames(Input_SettingsFile_Metab) == FALSE){
      warning("The provided file Input_SettingsFile_Metab must have a columns named: `Metabolite`.")
    }
  }

  ## ------------ Check Input SettingsInfo ----------- ##
  #3. Input_SettingsInfo
  if(Input_SettingsInfo[["conditions"]] %in% Input_SettingsInfo==TRUE){
    if(Input_SettingsInfo[["conditions"]] %in% colnames(Input_SettingsFile_Sample)== FALSE){
      stop("The ",Input_SettingsInfo[["conditions"]], " column selected as Conditions in Input_SettingsInfo was not found in Input_SettingsFile_Sample. Please check your input.")
    }else{# if true rename to Conditions
      Input_SettingsFile_Sample<- Input_SettingsFile_Sample%>%
        dplyr::rename("Conditions"= paste(Input_SettingsInfo[["conditions"]]) )
    }
  }else{
    stop("You have to provide a Input_SettingsInfo for conditions.")
  }

  ##########################
  if("denominator" %in% names(Input_SettingsInfo)==TRUE){
    if(Input_SettingsInfo[["denominator"]] %in% Input_SettingsFile_Sample$Conditions==FALSE){
      stop("The ",Input_SettingsInfo[["denominator"]], " column selected as denominator in Input_SettingsInfo was not found in Input_SettingsFile_Sample. Please check your input.")
    }else{
      denominator <- Input_SettingsInfo[["denominator"]]
    }
  }
  if("numerator" %in% names(Input_SettingsInfo)==TRUE){
    if(Input_SettingsInfo[["numerator"]] %in% Input_SettingsFile_Sample$Conditions  == FALSE){
      stop("The ",Input_SettingsInfo[["numerator"]], " column selected as numerator in Input_SettingsInfo was not found in Input_SettingsFile_Sample. Please check your input.")
    }else{
      numerator <- Input_SettingsInfo[["numerator"]]
    }
  }
  if("denominator" %in% names(Input_SettingsInfo)==FALSE  & "numerator" %in% names(Input_SettingsInfo) ==TRUE){
    stop("Check input. The selected denominator option is empty while ",paste(Input_SettingsInfo[["numerator"]])," has been selected as a numerator. Please add a denminator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison." )
  }

  ## ------------ Check Denominator/numerator ----------- ##
  #4.  Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("denominator" %in% names(Input_SettingsInfo)==FALSE  & "numerator" %in% names(Input_SettingsInfo) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = Input_SettingsFile_Sample$Conditions
    denominator <-unique(Input_SettingsFile_Sample$Conditions)
    numerator <-unique(Input_SettingsFile_Sample$Conditions)
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("denominator" %in% names(Input_SettingsInfo)==TRUE  & "numerator" %in% names(Input_SettingsInfo)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = Input_SettingsFile_Sample$Conditions
    denominator <- Input_SettingsInfo[["denominator"]]
    numerator <-unique(Input_SettingsFile_Sample$Conditions)
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }else if("denominator" %in% names(Input_SettingsInfo)==TRUE  & "numerator" %in% names(Input_SettingsInfo)==TRUE){
    # one-vs-one: Generate the comparisons
    comparisons <- matrix(c(denominator, numerator))
    #Settings:
    MultipleComparison = FALSE
    all_vs_all = FALSE
  }

  #5. Check if chosen test statistics fits with choice of comparison
  if(MultipleComparison==FALSE){
    STAT_pval_options <- c("t.test", "wilcox.test","chisq.test", "cor.test", "lmFit")
    if(STAT_pval %in% STAT_pval_options == FALSE){
      stop("Check input. The selected STAT_pval option for Hypothesis testing is not valid for multiple comparison (one-vs-all or all-vs-all). Please select one of the following: ",paste(STAT_pval_options,collapse = ", ")," or specify numerator and denumerator." )
    }
  }else{
    STAT_pval_options <- c("aov", "kruskal.test", "lmFit")
    if(STAT_pval %in% STAT_pval_options == FALSE){
      stop("Check input. The selected STAT_pval option for Hypothesis testing is not valid for one-vs-one comparsion. Multiple comparison is selected. Please select one of the following: ",paste(STAT_pval_options,collapse = ", ")," or change numerator and denumerator." )
    }
  }
  #if((STAT_pval =="wilcox-test" & nrow(Num)<5)|(STAT_pval =="wilcox-test" & nrow(Denom)<5)){# check number of samples for wilcoxons test
  #  warning("Number of samples measured per condition is <5 in at least one of the two conditions, which is small for using wilcox.test. Consider using another test.")
  #}

  ## ------------ Check General parameters ----------- ##
  #6. General parameters
  STAT_padj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(STAT_padj %in% STAT_padj_options == FALSE){
    stop("Check input. The selected STAT_padj option for multiple Hypothesis testing correction is not valid. Please select one of the folowing: ",paste(STAT_padj_options,collapse = ", "),"." )
  }
  if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for analysis of Consuption/Release experiment or =FALSE if not.")
  }
  if(is.logical(Plot) == FALSE){
    stop("Check input. The plot value should be either =TRUE if a Volcano plot presenting the DMA results is to be exported or =FALSE if not.")
  }
  Save_as_Plot_options <- c("svg","pdf","png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }
  Save_as_Results_options <- c("txt","csv", "xlsx" )
  if(Save_as_Results %in% Save_as_Results_options == FALSE){
    stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
  }

  #7. Are sample numbers enough?
  Num <- Input_data %>%
    filter(Input_SettingsFile_Sample$Conditions %in% numerator) %>%
    select_if(is.numeric)#only keep numeric columns with metabolite values
  Denom <- Input_data %>%
    filter(Input_SettingsFile_Sample$Conditions %in% denominator) %>%
    select_if(is.numeric)

  if(nrow(Num)==1){
    stop("There is only one sample available for ", numerator, ", so no statistical test can be performed.")
  } else if(nrow(Denom)==1){
    stop("There is only one sample available for ", denominator, ", so no statistical test can be performed.")
  }else if(nrow(Num)==0){
    stop("There is no sample available for ", numerator, ".")
  }else if(nrow(Denom)==0){
    stop("There is no sample available for ", denominator, ".")
  }

  ## ------------ Check Missingness ------------- ##
  #7.
  # If missing value imputation has not been performed the input data will most likely contain NA or 0 values for some metabolites, which will lead to Log2FC = NA.
  # Here we will check how many metabolites this affects in Num and Denom, and weather all replicates of a metabolite are affected.
  Num_Miss <- replace(Num, Num==0, NA)
  Num_Miss <- Num_Miss[, (colSums(is.na(Num_Miss)) > 0), drop = FALSE]

  Denom_Miss <- replace(Denom, Denom==0, NA)
  Denom_Miss <- Denom_Miss[, (colSums(is.na(Denom_Miss)) > 0), drop = FALSE]

  if((ncol(Num_Miss)>0 & ncol(Denom_Miss)==0)){
    message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(Num_Miss)
  } else if(ncol(Num_Miss)==0 & ncol(Denom_Miss)>0){
    message("In `denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(Denom_Miss)
  } else if(ncol(Num_Miss)>0 & ncol(Denom_Miss)>0){
    message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), " and in `denominator`",paste0(toString(denominator)), " ",ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "),". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  } else{
    message("There are no NA/0 values")
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  }


  ## ------------ Create Results output folder ----------- ##
  #8. Folders:
  if(is.null(Folder_Name)){
    name <- paste("MetaProViz_Results",Sys.Date(),sep = "_" )
  }else{
    if(grepl('[^[:alnum:]]', Folder_Name)){
      stop("The 'Folder_Name' must not contain any special character.")
    }else{
      name <- paste("MetaProViz_Results",Sys.Date(),Folder_Name,sep = "_" )
    }
  }
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  Results_folder_DMA_folder <- file.path(Results_folder,"DMA") # Make DMA results folder
  if (!dir.exists(Results_folder_DMA_folder)) {dir.create(Results_folder_DMA_folder)}
  Results_folder_DMA_folder_Shapiro_folder <- file.path(Results_folder_DMA_folder,"Shapiro") # Make DMA results folder
  if (!dir.exists(Results_folder_DMA_folder_Shapiro_folder)) {dir.create(Results_folder_DMA_folder_Shapiro_folder)}
  Results_folder_Conditions <- file.path(Results_folder_DMA_folder,paste0(toString(numerator),"_vs_",toString(denominator))) # Make comparison folder
  if (!dir.exists(Results_folder_Conditions)) {dir.create(Results_folder_Conditions)}

  # Prepare output names
  if(OutputName==""){
    OutputName <- OutputName
  }else{
    OutputName <- paste("_",OutputName, sep="")
  }
  if(CoRe==TRUE){
    OutputName <- paste(OutputName,"_CoRe", sep="")
  }


  # Here is the Shapiro
  Shapiro_output <-suppressWarnings(MetaProViz:::Shapiro(Input_data=Input_data,
                                                         Input_SettingsFile_Sample=Input_SettingsFile_Sample,
                                                         Input_SettingsInfo=Input_SettingsInfo,
                                                         STAT_pval=STAT_pval,
                                                         CoRe=CoRe,
                                                         OutputName=OutputName,
                                                         Save_as_Plot=Save_as_Plot,
                                                         QQplots=FALSE,
                                                         Save_as_Results=Save_as_Results,
                                                         Plot=FALSE,
                                                         Folder_Name=Results_folder_DMA_folder_Shapiro_folder))

  ###############################################################################################################################################################################################################
  #### Prepare the data ######
  #1. Metabolite names:
  savedMetaboliteNames <-  data.frame("InnputName"=colnames(Input_data))
  savedMetaboliteNames$Metabolite <- paste0("M", seq(1,length(colnames(Input_data))))
  colnames(Input_data) <- savedMetaboliteNames$Metabolite

  ################################################################################################################################################################################################
  ############### Calculate Log2FC, pval, padj, tval and add additional info ###############
  Log2FC_table <- list()# Create an empty list to store results data frames
  for(column in 1:dim(comparisons)[2]){
    C1 <- Input_data %>% # Numerator
      filter(Input_SettingsFile_Sample$Conditions %in% comparisons[1,column]) %>%
      select_if(is.numeric)#only keep numeric columns with metabolite values
    C2 <- Input_data %>% # Deniminator
      filter(Input_SettingsFile_Sample$Conditions %in%  comparisons[2,column]) %>%
      select_if(is.numeric)

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

    if(CoRe==TRUE){#Calculate absolute distance between the means. log2 transform and add sign (-/+):
      #CoRe values can be negative and positive, which can does not allow us to calculate a Log2FC.
      Mean_C1_t <- as.data.frame(t(Mean_C1))%>%
        rownames_to_column("Metabolite")
      Mean_C2_t <- as.data.frame(t(Mean_C2))%>%
        rownames_to_column("Metabolite")
      Mean_Merge <-merge(Mean_C1_t, Mean_C2_t, by="Metabolite", all=TRUE)%>%
        rename("C1"=2,
               "C2"=3)

      #Deal with NA/0s
      Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss#Column to enable the check if mean values of 0 are due to missing values (NA/0) and not by coincidence

      if(any((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))==TRUE){
        Mean_Merge <- Mean_Merge%>%
          mutate(C1 = case_when(C2 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                C1 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                C2 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                C1 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                TRUE ~ paste(C1)))%>%
          mutate(C2 = case_when(C1 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                C2 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                C1 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                C2 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                TRUE ~ paste(C2)))%>%
          mutate(C1 = as.numeric(C1))%>%
          mutate(C2 = as.numeric(C2))

        X <- Mean_Merge%>%
          subset((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))
        message("We added +1 to the mean value of metabolite(s) ", paste0(X$Metabolite, collapse = ", "), ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
      }

      #Add the distance column:
      Mean_Merge$`Log2(Distance)` <-log2(abs(Mean_Merge$C1 - Mean_Merge$C2))

      Mean_Merge <- Mean_Merge%>%#Now we can adapt the values to take into account the distance
        mutate(`Log2(Distance)` = case_when(C1 > C2 ~ paste(`Log2(Distance)`*+1),#If C1>C2 the distance stays positive to reflect that C1 > C2
                                            C1 < C2 ~ paste(`Log2(Distance)`*-1),#If C1<C2 the distance gets a negative sign to reflect that C1 < C2
                                            TRUE ~ 'NA'))%>%
        mutate(`Log2(Distance)` = as.numeric(`Log2(Distance)`))

      #Add additional information:
      temp1 <- Mean_C1
      temp2 <- Mean_C2
      #Add Info of CoRe:
      CoRe_info <- rbind(temp1, temp2,rep(0,length(temp1)))
      for (i in 1:length(temp1)){
        if (temp1[i]>0 & temp2[i]>0){
          CoRe_info[3,i] <- "Released"
        }else if (temp1[i]<0 & temp2[i]<0){
          CoRe_info[3,i] <- "Consumed"
        }else if(temp1[i]>0 & temp2[i]<0){
          CoRe_info[3,i] <- paste("Released in" ,comparisons[1,column] , "and Consumed",comparisons[2,column] , sep=" ")
        } else if(temp1[i]<0 & temp2[i]>0){
          CoRe_info[3,i] <- paste("Consumed in" ,comparisons[1,column] , " and Released",comparisons[2,column] , sep=" ")
        }else{
          CoRe_info[3,i] <- "No Change"
        }
      }

      CoRe_info <- t(CoRe_info) %>% as.data.frame()
      CoRe_info <- rownames_to_column(CoRe_info, "Metabolite")
      names(CoRe_info)[2] <- paste("Mean",  comparisons[1,column], sep="_")
      names(CoRe_info)[3] <- paste("Mean",  comparisons[2,column], sep="_")
      names(CoRe_info)[4] <- "CoRe_specific"

      CoRe_info <-CoRe_info%>%
        mutate(CoRe = case_when(CoRe_specific == "Released" ~ 'Released',
                                CoRe_specific == "Consumed" ~ 'Consumed',
                                TRUE ~ 'Released/Consumed'))

      Log2FC_C1vC2 <-merge(Mean_Merge[,c(1,5)], CoRe_info[,c(1,4:5,2:3)], by="Metabolite", all.x=TRUE)

      #Add info on Input:
      temp3 <- as.data.frame(t(C1))%>%rownames_to_column("Metabolite")
      temp4 <- as.data.frame(t(C2))%>%rownames_to_column("Metabolite")
      temp_3a4 <-merge(temp3, temp4, by="Metabolite", all=TRUE)
      Log2FC_C1vC2 <-merge(Log2FC_C1vC2, temp_3a4, by="Metabolite", all.x=TRUE)

      #Add info on Pathways:
      if(is.null(Input_SettingsFile_Metab)!=TRUE & 'Metabolite' %in% colnames(Input_SettingsFile_Metab)){
        Pathways <- merge(savedMetaboliteNames , Input_SettingsFile_Metab, by.x="InnputName", by.y="Metabolite", all.y=TRUE)
        Log2FC_C1vC2<- merge(Log2FC_C1vC2, Pathways[,-c(1)],by="Metabolite", all.x=T)
      }

      #Return DFs
      ##Make reverse DF
      Log2FC_C2vC1 <-Log2FC_C1vC2
      Log2FC_C2vC1$`Log2(Distance)` <- Log2FC_C2vC1$`Log2(Distance)` *-1

      ##Save them
      if(MultipleComparison == TRUE){
        logname <- paste(comparisons[1,column], comparisons[2,column],sep="_vs_")
        logname_reverse <- paste(comparisons[2,column], comparisons[1,column],sep="_vs_")

        # Store the data frame in the results list, named after the contrast
        Log2FC_table[[logname]] <- Log2FC_C1vC2
        Log2FC_table[[logname_reverse]] <- Log2FC_C2vC1
      }else{
         Log2FC_table <- Log2FC_C1vC2
      }
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

      Mean_Merge <- Mean_Merge%>%
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

      #Add info on Input:
      temp3 <- as.data.frame(t(C1))%>%rownames_to_column("Metabolite")
      temp4 <- as.data.frame(t(C2))%>%rownames_to_column("Metabolite")
      temp_3a4 <-merge(temp3, temp4, by="Metabolite", all=TRUE)
      Log2FC_C1vC2 <-merge(Mean_Merge[,c(1,8)], temp_3a4, by="Metabolite", all.x=TRUE)

      #Add info on Pathways:
      if(is.null(Input_SettingsFile_Metab)!=TRUE & 'Metabolite' %in% colnames(Input_SettingsFile_Metab)){
        Pathways <- merge(savedMetaboliteNames , Input_SettingsFile_Metab, by.x="InnputName", by.y="Metabolite", all.y=TRUE)
        Log2FC_C1vC2<- merge(Log2FC_C1vC2, Pathways[,-c(1)],by="Metabolite", all.x=T)
      }

      #Return DFs
      ##Make reverse DF
      Log2FC_C2vC1 <-Log2FC_C1vC2
      Log2FC_C2vC1$Log2FC <- Log2FC_C2vC1$Log2FC*-1

      if(MultipleComparison == TRUE){
        logname <- paste(comparisons[1,column], comparisons[2,column],sep="_vs_")
        logname_reverse <- paste(comparisons[2,column], comparisons[1,column],sep="_vs_")


        # Store the data frame in the results list, named after the contrast
        Log2FC_table[[logname]] <- Log2FC_C1vC2
        Log2FC_table[[logname_reverse]] <- Log2FC_C2vC1
      }else{
        Log2FC_table <- Log2FC_C1vC2
      }
    }else{
      stop("Please choose CoRe= TRUE or CoRe=FALSE.")
    }
  }

  ################################################################################################################################################################################################
  ############### Perform Hypothesis testing ###############
  if(MultipleComparison == FALSE){
    if(STAT_pval=="lmFit"){
      STAT_C1vC2 <- DMA_Stat_limma(Input_data=Input_data,
                                   Input_SettingsFile_Sample=Input_SettingsFile_Sample,
                                   Input_SettingsInfo=Input_SettingsInfo,
                                   STAT_padj=STAT_padj,
                                   Log2FC_table=Log2FC_table,
                                   CoRe=CoRe,
                                   all_vs_all=all_vs_all,
                                   MultipleComparison=MultipleComparison)
    }else{
      STAT_C1vC2 <-MetaProViz:::DMA_Stat_single(C1=C1, C2=C2, Log2FC_table=Log2FC_table, Metabolites_Miss=Metabolites_Miss, STAT_pval=STAT_pval, STAT_padj=STAT_padj)
    }
    }else{ # MultipleComparison = TRUE
      # conditions =as.factor(conditions)
      if(all_vs_all ==TRUE){
        message("No conditions were specified as numerator or denumerator. Performing multiple testing `all-vs-all` using ", paste(STAT_pval), ".")
        }else{# for 1 vs all
          message("No condition was specified as numerator and ",paste(denominator), " was selected as a denominator. Performing multiple testing `all-vs-one` using ", paste(STAT_pval), ".")
          # conditions=relevel(conditions, ref = denominator)
        }

      if(STAT_pval=="aov"){
        STAT_C1vC2 <- AOV(Input_data=Input_data,
                          Input_SettingsInfo=Input_SettingsInfo,
                          conditions=conditions,
                          STAT_padj=STAT_padj,
                          Log2FC_table=Log2FC_table,
                          all_vs_all=all_vs_all,
                          comparisons=comparisons)
        }else if(STAT_pval=="kruskal.test"){
          STAT_C1vC2 <-Kruskal(Input_data=Input_data,
                              conditions=conditions,
                              STAT_padj=STAT_padj,
                              Log2FC_table=Log2FC_table,
                              all_vs_all=all_vs_all,
                              comparisons=comparisons)
          }else if(STAT_pval=="lmFit"){
            STAT_C1vC2 <- DMA_Stat_limma(Input_data=Input_data,
                                         Input_SettingsFile_Sample=Input_SettingsFile_Sample,
                                         Input_SettingsInfo=Input_SettingsInfo,
                                         STAT_padj=STAT_padj,
                                         Log2FC_table=Log2FC_table,
                                         CoRe=CoRe,
                                         all_vs_all=all_vs_all,
                                         MultipleComparison=MultipleComparison)
          }
      }

  ################################################################################################################################################################################################
  ###############  Add the previous metabolite names back ###############
  if(MultipleComparison==FALSE){
    DMA_Output <- merge(savedMetaboliteNames, STAT_C1vC2, by="Metabolite")
    DMA_Output$Metabolite <- NULL
    colnames(DMA_Output)[1] <- "Metabolite"
  }else{
    DMA_Output <- lapply(STAT_C1vC2, function(df){
      merged_df <- merge(savedMetaboliteNames, df, by = "Metabolite", all.y = TRUE)
      merged_df <-merged_df[,-1]%>%#remove the names we used as part of the function and add back the input names.
        dplyr::rename("Metabolite"=1)
      return(merged_df)
    })
  }

  ################################################################################################################################################################################################
  ###############  Folder ###############
  if(MultipleComparison==FALSE){
    if (Save_as_Results == "xlsx"){
      xlsDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator), OutputName, ".xlsx"))   # Save the DMA results table
      writexl::write_xlsx(DMA_Output,xlsDMA, col_names = TRUE) # save the DMA result DF
      }else if (Save_as_Results == "csv"){
        csvDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator), OutputName, ".csv"))
        write.csv(DMA_Output,csvDMA) # save the DMA result DF
      }else if (Save_as_Results == "txt"){
        txtDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator), OutputName, ".txt"))
        write.table(DMA_Output,txtDMA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
      }
  }else{
    for(DF in names(DMA_Output)){
      DMA_Output_Save <- DMA_Output[[DF]]
      DF_save <- gsub("[^A-Za-z0-9._-]", "_", DF)## Remove special characters and replace spaces with underscores

      if (Save_as_Results == "xlsx"){
        xlsDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",DF_save, OutputName, ".xlsx"))   # Save the DMA results table
        writexl::write_xlsx(DMA_Output_Save,xlsDMA, col_names = TRUE) # save the DMA result DF
      }else if (Save_as_Results == "csv"){
        csvDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",DF_save, OutputName, ".csv"))
        write.csv(DMA_Output_Save,csvDMA) # save the DMA result DF
      }else if (Save_as_Results == "txt"){
        txtDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",DF_save, OutputName, ".txt"))
        write.table(DMA_Output_Save,txtDMA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
      }
    }
  }

  if(CoRe==TRUE){
    x <- "Log2(Distance)"
    VolPlot_SettingsInfo= c(color="CoRe")
    VOlPlot_SettingsFile = DMA_Output
  }else{
    x <- "Log2FC"
    VolPlot_SettingsInfo= NULL
    VOlPlot_SettingsFile = NULL
  }

  ################################################################################################################################################################################################
  ###############  Plots ###############
  volplotList = list()
  if(MultipleComparison==TRUE){
    for(DF in names(DMA_Output)){
      Volplotdata<- DMA_Output[[DF]]

      dev.new()
      VolcanoPlot <- invisible(MetaProViz::VizVolcano(Plot_Settings="Standard",
                                                      Input_data=Volplotdata,
                                                      Plot_SettingsInfo=VolPlot_SettingsInfo,
                                                      Plot_SettingsFile=VOlPlot_SettingsFile[[DF]],
                                                      y= "p.adj",
                                                      x= x,
                                                      AdditionalInput_data= NULL,
                                                      OutputPlotName= DF,
                                                      Comparison_name= c(Input_data="Cond1", AdditionalInput_data= "Cond2"),
                                                      xlab= NULL,#"~Log[2]~FC"
                                                      ylab= NULL,#"~-Log[10]~p.adj"
                                                      pCutoff= 0.05,
                                                      FCcutoff= 0.5,
                                                      color_palette= NULL,
                                                      shape_palette=NULL,
                                                      SelectLab= "",
                                                      Connectors=  FALSE,
                                                      Subtitle=  bquote(italic("Differential Metabolite Analysis")),
                                                      Theme= NULL,
                                                      Save_as_Plot= NULL))

      volplotList[[DF]]<- VolcanoPlot

      DF_save <- gsub("[^A-Za-z0-9._-]", "_", DF)## Remove special characters and replace spaces with underscores
      volcanoDMA <- file.path(Results_folder_Conditions,paste0( "Volcano_Plot_", DF_save, OutputName,".",Save_as_Plot))
      ggsave(volcanoDMA,plot=VolcanoPlot, width=10, height=8) # save the volcano plot
      dev.off()
    }
  }else{
  # Make a simple Volcano plot
  dev.new()
  VolcanoPlot <- invisible(MetaProViz::VizVolcano(Plot_Settings="Standard",
                                                  Input_data=DMA_Output,
                                                  Plot_SettingsInfo=VolPlot_SettingsInfo,
                                                  Plot_SettingsFile=VOlPlot_SettingsFile,
                                                  y= "p.adj",
                                                  x= x,
                                                  AdditionalInput_data= NULL,
                                                  OutputPlotName= paste0(toString(numerator)," versus ",toString(denominator)),
                                                  Comparison_name= c(Input_data="Cond1", AdditionalInput_data= "Cond2"),
                                                  xlab= NULL,#"~Log[2]~FC"
                                                  ylab= NULL,#"~-Log[10]~p.adj"
                                                  pCutoff= 0.05,
                                                  FCcutoff= 0.5,
                                                  color_palette= NULL,
                                                  shape_palette=NULL,
                                                  SelectLab= "",
                                                  Connectors=  FALSE,
                                                  Subtitle=  bquote(italic("Differential Metabolite Analysis")),
                                                  Theme= NULL,
                                                  Save_as_Plot= NULL))
  volplotList[[paste0(toString(numerator)," versus ",toString(denominator))]]<- VolcanoPlot
  dev.off()
  #plot(VolcanoPlot)

  volcanoDMA <- file.path(Results_folder_Conditions,paste0( "Volcano_Plot_",toString(numerator),"_versus_",toString(denominator),OutputName,".",Save_as_Plot))
  ggsave(volcanoDMA,plot=VolcanoPlot, width=10, height=8) # save the volcano plot
  }


  #Here we make a list in which we will save the output
  suppressWarnings(DMA_output_list <- list("DF" = list("Shapiro_result"=Shapiro_output$DF$Shapiro_result,"DMA_result"=DMA_Output),"Plot"=list( "Distributions"=Shapiro_output$Plot$Distributions, "Volcano"=volplotList)))


  if(Plot == TRUE){
    for (plot in volplotList){
      print(plot)
    }
  }
  return(invisible(DMA_output_list))
}



##########################################################################################
### ### ### DMA helper function: Internal Function to perform single comparison ### ### ###
##########################################################################################

#' @param C1 This is the C1 (Condition 1) DF generated within the DMA function.
#' @param C2 This is the C2 (Condition 2) DF generated within the DMA function.
#' @param Log2FC_table this is the Log2FC DF generated within the DMA function.
#' @param Metabolites_Miss these are the metabolites with missing values generated within the DMA function.
#' @param STAT_pval Passed to DMA
#' @param STAT_padj Passed to DMA
#'
#' @keywords DMA helper function
#' @noRd
#'

DMA_Stat_single <- function(C1, C2, Log2FC_table, Metabolites_Miss, STAT_pval, STAT_padj){
  ## ------------ Perform Hypothesis testing ----------- ##
  # For C1 and C2 we use 0, since otherwise we can not perform the statistical testing.
  C1[is.na(C1)] <- 0
  C2[is.na(C2)] <- 0

  #### 1. p.value and test statistics (=t.val)
  T_C1vC2 <-mapply(STAT_pval, x= as.data.frame(C2), y = as.data.frame(C1), SIMPLIFY = F)

  VecPVAL_C1vC2 <- c()
  VecTVAL_C1vC2 <- c()
  for(i in 1:length(T_C1vC2)){
    p_value <- unlist(T_C1vC2[[i]][3])
    t_value <- unlist(T_C1vC2[[i]])[1]   # Extract the t-value
    VecPVAL_C1vC2[i] <- p_value
    VecTVAL_C1vC2[i] <- t_value
  }
  Metabolite <- colnames(C2)
  PVal_C1vC2 <- data.frame(Metabolite, p.val = VecPVAL_C1vC2, t.val = VecTVAL_C1vC2)

  #we set p.val= NA, for metabolites that had 1 or more replicates with NA/0 values and remove them prior to p-value adjustment
  PVal_C1vC2$`NA/0` <- PVal_C1vC2$Metabolite %in% Metabolites_Miss
  PVal_C1vC2 <-PVal_C1vC2%>%
    mutate(p.val = case_when(`NA/0`== TRUE ~ NA,
                             TRUE ~ paste(VecPVAL_C1vC2)))
  PVal_C1vC2$p.val = as.numeric(as.character(PVal_C1vC2$p.val))

  #### 2. p.adjusted
  #Split data for p.value adjustment to exclude NA
  PVal_NA <- PVal_C1vC2[is.na(PVal_C1vC2$p.val), c(1:3)]
  PVal_C1vC2 <-PVal_C1vC2[!is.na(PVal_C1vC2$p.val), c(1:3)]

  #perform adjustment
  VecPADJ_C1vC2 <- p.adjust((PVal_C1vC2[,2]),method = STAT_padj, n = length((PVal_C1vC2[,2]))) #p-adjusted
  Metabolite <- PVal_C1vC2[,1]
  PADJ_C1vC2 <- data.frame(Metabolite, p.adj = VecPADJ_C1vC2)
  STAT_C1vC2 <- merge(PVal_C1vC2,PADJ_C1vC2, by="Metabolite")

  #Add Metabolites that have p.val=NA back into the DF for completeness.
  if(nrow(PVal_NA)>0){
    PVal_NA$p.adj <- NA
    STAT_C1vC2 <- rbind(STAT_C1vC2, PVal_NA)
  }

  #Add Log2FC
  STAT_C1vC2 <- merge(Log2FC_table,STAT_C1vC2[,c(1:2,4,3)], by="Metabolite")

  #order for t.value
  STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value

  Output <- STAT_C1vC2
}


# all-vs-all:
################################################################
### ### ### AOV: Internal Function to perform anova  ### ### ###
################################################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Input_SettingsInfo Passed to DMA
#' @param conditions Factor with sample group information.
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param Log2FC_table Table with Metabolites are rows and a Log2FC column \strong(Default = Log2FC_table)
#' @param all_vs_all  True for multiple comparison of all against all or FALSE for multiple comparison of all against one selected base.
#' @param comparisons Dataframe containing the comparison information for multiple comparison. The first row is the numerator, the second row is the denominator and each column is a different comparison.
#'
#' @keywords Kruskal test,Hypothesis testing, p.value
#' @export


AOV <-function(Input_data,
               Input_SettingsInfo,
               conditions,
               STAT_padj,
               Log2FC_table,
               all_vs_all,
               comparisons
){
  ## 1. Anova p.val
  aov.res= apply(Input_data,2,function(x) aov(x~conditions))

  ## 2. Tukey test p.adj
  posthoc.res = lapply(aov.res, TukeyHSD, conf.level=0.95)
  Tukey_res <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'p adj']))
  Tukey_res <- as.data.frame(Tukey_res)

  comps <-   paste(comparisons[1, ], comparisons[2, ], sep="-")# normal
  opp_comps <-  paste(comparisons[2, ], comparisons[1, ], sep="-")

  if(sum(opp_comps %in%  colnames(Tukey_res))>0){# if opposite comparisons is true
    for (comp in 1: length(opp_comps)){
      colnames(Tukey_res)[colnames(Tukey_res) %in% opp_comps[comp]] <-  comps[comp]
    }
  }

  ## 3. t.val
  Tukey_res_diff <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'diff'])) %>% as.data.frame()

  if (sum(opp_comps %in%  colnames(Tukey_res_diff))>0){# if oposite comparisons is true
    for (comp in 1: length(opp_comps)){
      colnames(Tukey_res_diff)[colnames(Tukey_res_diff) %in% opp_comps[comp]] <-  comps[comp]
    }
  }

  #Make output DFs:
  Pval_table <- Tukey_res
  Pval_table <- rownames_to_column(Pval_table,"Metabolite")

  Tval_table <- rownames_to_column(Tukey_res_diff,"Metabolite")

  common_col_names <- setdiff(names(Tukey_res_diff), "row.names")#Here we need to adapt for one_vs_all or all_vs_all

  results_list <- list()
  for(col_name in common_col_names){
    # Create a new data frame by merging the two data frames
    merged_df <- merge(Pval_table[,c("Metabolite",col_name)], Tval_table[,c("Metabolite",col_name)], by="Metabolite", all=TRUE)%>%
      dplyr::rename("p.adj"=2,
                    "t.val"=3)

    #We need to add _vs_ into the comparison col_name
    pattern <- paste(conditions, collapse = "|")
    conditions_present <- unique(unlist(regmatches(col_name, gregexpr(pattern, col_name))))
    modified_col_name <- paste(conditions_present[1], "vs", conditions_present[2], sep = "_")

    # Add the new data frame to the list with the column name as the list element name
    results_list[[modified_col_name]] <- merged_df
  }

  # Merge the data frames in list1 and list2 based on the "Metabolite" column
  list_names <-  names(results_list)

  merged_list <- list()
  for(name in list_names){
    # Check if the data frames exist in both lists
    if(name %in% names(results_list) && name %in% names(Log2FC_table)){
      merged_df <- merge(results_list[[name]], Log2FC_table[[name]], by = "Metabolite", all = TRUE)
      merged_df <- merged_df[,c(1,4,2:3,5:ncol(merged_df))]#reorder the columns
      merged_list[[name]] <- merged_df
    }
  }

  # Make sure the right comparisons are returned:
  if(all_vs_all==TRUE){
    STAT_C1vC2 <- merged_list
  }else if(all_vs_all==FALSE){
    #remove the comparisons that are not needed:
    modified_df_list <- list()
    for(df_name in names(merged_list)){
      if(endsWith(df_name, Input_SettingsInfo[["denominator"]])){
        modified_df_list[[df_name]] <- merged_list[[df_name]]
      }
    }
    STAT_C1vC2 <- modified_df_list
  }

  return(STAT_C1vC2)
}

###########################################################################
### ### ### Kruskal: Internal Function to perform Kruskal test  ### ### ###
###########################################################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param conditions Factor with sample group information.
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param Log2FC_table Table with Metabolites are rows and a Log2FC column \strong(Default = Log2FC_table)
#' @param all_vs_all  True for multiple comparison of all against all or FALSE for multiple comparison of all against one selected base.
#' @param comparisons Dataframe containing the comparison information for multiple comparison. The first row is the numerator, the second row is the denominator and each column is a different comparison.
#'
#' @keywords Kruskal test,Hypothesis testing, p.value
#' @export


Kruskal <-function(Input_data,
                   conditions,
                   STAT_padj,
                   Log2FC_table,
                   all_vs_all,
                   comparisons
){
  # Kruskal test (p.val)
  aov.res= apply(Input_data,2,function(x) kruskal.test(x~conditions))
  anova_res<-do.call('rbind', lapply(aov.res, function(x) {x["p.value"]}))
  anova_res <- as.matrix(mutate_all(as.data.frame(anova_res), function(x) as.numeric(as.character(x))))
  colnames(anova_res) = c("Kruskal_p.val")

  # Dunn test (p.adj)
  Dunndata <- Input_data %>%
    mutate(conditions = conditions) %>%
    select(conditions, everything())%>% as.data.frame()

  # Applying a loop to obtain p.adj and t.val:
  Dunn_Pres<- data.frame(comparisons = paste(comparisons[1,],    comparisons[2,], sep = "_vs_" ))
  Dunn_Tres<- Dunn_Pres
  for(col in 2:dim(Dunndata)[2]){
    data = Dunndata[,c(1,col)]
    colnames(data)[2] <- gsub("^\\d+", "", colnames(data)[2])

    ## If a metabolite starts with number remove it
    formula <- as.formula(paste(colnames(data)[2], "~ conditions"))
    posthoc.res= rstatix::dunn_test(data, formula, p.adjust.method = STAT_padj)

    res <- data.frame(comparisons =  paste(posthoc.res$group1, posthoc.res$group2, sep = "_vs_" ))
    res[[colnames(Dunndata)[col] ]] <-  posthoc.res$p.adj
    res <- res[res$comparisons %in% Dunn_Pres$comparisons ,] # take only the comparisons selected
    Dunn_Pres <- merge(Dunn_Pres,res,by="comparisons")

    tvalues <- data.frame(comparisons =  paste(posthoc.res$group1, posthoc.res$group2, sep = "_vs_" ))
    tvalues[[colnames(Dunndata)[col] ]] <-  posthoc.res$statistic
    tvalues <- tvalues[tvalues$comparisons %in% Dunn_Pres$comparisons ,]
    Dunn_Tres <- merge(Dunn_Tres,tvalues,by="comparisons")
  }

  #Make output DFs:
  Dunn_Pres <- column_to_rownames(Dunn_Pres, "comparisons")%>% t() %>% as.data.frame()
  Pval_table <- as.matrix(mutate_all(as.data.frame(Dunn_Pres), function(x) as.numeric(as.character(x))))
  Pval_table <- Pval_table %>% as.data.frame()
  Pval_table <- rownames_to_column(Pval_table, "Metabolite")

  Dunn_Tres <- column_to_rownames(Dunn_Tres, "comparisons")%>% t() %>% as.data.frame()
  Tval_table <- as.matrix(mutate_all(as.data.frame(Dunn_Tres), function(x) as.numeric(as.character(x))))
  Tval_table <- Tval_table %>% as.data.frame()
  Tval_table <- rownames_to_column(Tval_table, "Metabolite")

  common_col_names <- setdiff(names(Dunn_Pres), "row.names")

  results_list <- list()
  for(col_name in common_col_names){
    # Create a new data frame by merging the two data frames
    merged_df <- merge(Pval_table[,c("Metabolite",col_name)], Tval_table[,c("Metabolite",col_name)], by="Metabolite", all=TRUE)%>%
      dplyr::rename("p.adj"=2,
                    "t.val"=3)
    # Add the new data frame to the list with the column name as the list element name
    results_list[[col_name]] <- merged_df
  }

  # Merge the data frames in list1 and list2 based on the "Metabolite" column
  merged_list <- list()
  for(name in common_col_names){
    # Check if the data frames exist in both lists
    if(name %in% names(results_list) && name %in% names(Log2FC_table)){
      merged_df <- merge(results_list[[name]], Log2FC_table[[name]], by = "Metabolite", all = TRUE)
      merged_df <- merged_df[,c(1,4,2:3,5:ncol(merged_df))]#reorder the columns
      merged_list[[name]] <- merged_df
    }
  }

  STAT_C1vC2 <- merged_list

  return(STAT_C1vC2)
}

##########################################################################################
### ### ### DMA helper function: Internal Function to perform limma ### ### ###
##########################################################################################

#' @param Input_data Passed to DMA
#' @param Input_SettingsFile_Sample Passed to DMA
#' @param Input_SettingsInfo Passed to DMA
#' @param Log2FC_table this is the Log2FC DF generated within the DMA function.
#' @param STAT_padj Passed to DMA
#' @param CoRe Passed to DMA
#' @param all_vs_all generated within the DMA function
#' @param MultipleComparison generated within the DMA function
#'
#' @keywords DMA helper function
#' @noRd
#'

DMA_Stat_limma <- function(Input_data, Input_SettingsFile_Sample, Input_SettingsInfo, Log2FC_table, STAT_padj, CoRe, all_vs_all, MultipleComparison){
  ####------ Ensure that Input_data is ordered by conditions and sample names are the same as in Input_SettingsFile_Sample:
  targets <- Input_SettingsFile_Sample%>%
    rownames_to_column("sample")
  targets<- targets[,c("sample", Input_SettingsInfo[["conditions"]])]%>%
    dplyr::rename("condition"=2)%>%
    arrange(sample)#Order the column "sample" alphabetically
  targets$condition_limma_compatible <-make.names(targets$condition)#make appropriate condition names accepted by limma

  if(MultipleComparison==FALSE){
    #subset the data:
    targets<-targets%>%
      subset(condition==Input_SettingsInfo[["numerator"]] | condition==Input_SettingsInfo[["denominator"]])%>%
      arrange(sample)#Order the column "sample" alphabetically

    Limma_input <- Input_data%>%rownames_to_column("sample")
    Limma_input <-merge(targets[,1:2],  Limma_input, by="sample", all.x=TRUE)
    Limma_input <- Limma_input[,-2]%>%
      arrange(sample)#Order the column "sample" alphabetically
  }else if(MultipleComparison==TRUE){
    Limma_input <- Input_data%>%rownames_to_column("sample")%>%
    arrange(sample)#Order the column "sample" alphabetically
  }

  #Check if the order of the "sample" column is the same in both data frames
  if(identical(targets$sample, Limma_input$sample)==FALSE){
    stop("The order of the 'sample' column is different in both data frames. Please make sure that Input_SettingsFile_Sample and Input_data contain the same rownames and sample numbers.")
  }

  targets_limma <-targets[,-2]%>%
    dplyr::rename("condition"="condition_limma_compatible")

  #We need to transpose the df to run limma. Also, if the data is not log2 transformed, we will not calculate the Log2FC as limma just substracts one condition from the other
  Limma_input <- as.data.frame(t(Limma_input%>%column_to_rownames("sample")))
  Limma_input_log2 <- log2(Limma_input) # communicate the log2 transformation --> how does limma deals with NA when calculating the change?

  #### ------Run limma:
  ####  Make design matrix:
  fcond <- as.factor(targets_limma$condition)#all versus all

  design <- model.matrix(~0 + fcond)# Create the design matrix
  colnames(design) <- levels(fcond) # Give meaningful column names to the design matrix

  #### Fit the linear model
  fit <- limma::lmFit(Limma_input_log2, design)

  ####  Make contrast matrix:
  if(all_vs_all ==TRUE & MultipleComparison==TRUE){
    unique_conditions <- levels(fcond)# Get unique conditions

    # Create an empty contrast matrix
    num_conditions <- length(unique_conditions)
    num_comparisons <- num_conditions * (num_conditions - 1) / 2
    cont.matrix <- matrix(0, nrow = num_comparisons, ncol = num_conditions)

    # Initialize an index for the column in the contrast matrix
    i <- 1

    # Initialize column and row names
    colnames(cont.matrix) <- unique_conditions
    rownames(cont.matrix) <- character(num_comparisons)

    # Loop through all pairwise combinations of unique conditions
    for (condition1 in 1:(num_conditions - 1)) {
      for (condition2 in (condition1 + 1):num_conditions) {
        # Create the pairwise comparison vector
        comparison <- rep(0, num_conditions)

        comparison[condition2] <- -1
        # Add the comparison vector to the contrast matrix
        cont.matrix[i, ] <- comparison
        # Set row name
        rownames(cont.matrix)[i] <- paste(unique_conditions[condition1], "_vs_", unique_conditions[condition2], sep="")
        i <- i + 1
      }
    }
  }else if(all_vs_all ==FALSE & MultipleComparison==TRUE){
    unique_conditions <- levels(fcond)# Get unique conditions
    denominator  <- Input_SettingsInfo[["denominator"]]

    # Create an empty contrast matrix
    num_conditions <- length(unique_conditions)
    num_comparisons <- num_conditions - 1
    cont.matrix <- matrix(0, nrow = num_comparisons, ncol = num_conditions)


    # Initialize an index for the column in the contrast matrix
    i <- 1

    # Initialize column and row names
    colnames(cont.matrix) <- unique_conditions
    rownames(cont.matrix) <- character(num_comparisons)

    # Loop through all pairwise combinations of unique conditions
    for(condition in 2:num_conditions){
      # Create the pairwise comparison vector
      comparison <- rep(0, num_conditions)
      if(unique_conditions[1]==Input_SettingsInfo[["denominator"]]){
        comparison[1] <- -1
        comparison[condition] <- 1
        # Add the comparison vector to the contrast matrix
        cont.matrix[i, ] <- comparison
        # Set row name
        rownames(cont.matrix)[i] <- paste(unique_conditions[condition], "_vs_", unique_conditions[1], sep = "")
      }else{
        comparison[1] <- 1
        comparison[condition] <- -1
        # Add the comparison vector to the contrast matrix
        cont.matrix[i, ] <- comparison
        # Set row name
        rownames(cont.matrix)[i] <- paste(unique_conditions[1], "_vs_", unique_conditions[condition], sep = "")

        }
   i <- i + 1
    }
  }else if(all_vs_all ==FALSE & MultipleComparison==FALSE){
    unique_conditions <- levels(fcond)# Get unique conditions

    # Create an empty contrast matrix
    num_conditions <- length(unique_conditions)
    num_comparisons <- num_conditions * (num_conditions - 1) / 2
    cont.matrix <- matrix(0, nrow = num_comparisons, ncol = num_conditions)

    # Initialize an index for the column in the contrast matrix
    i <- 1

    # Initialize column and row names
    colnames(cont.matrix) <- unique_conditions
    rownames(cont.matrix) <- character(num_comparisons)

    # Loop through all pairwise combinations of unique conditions
    for (condition1 in 1:(num_conditions - 1)) {
      for (condition2 in (condition1 + 1):num_conditions) {
        # Create the pairwise comparison vector
        comparison <- rep(0, num_conditions)
        if(unique_conditions[condition1]==Input_SettingsInfo[["denominator"]]){
          comparison[condition1] <- -1
          comparison[condition2] <- 1
          # Add the comparison vector to the contrast matrix
          cont.matrix[i, ] <- comparison
          # Set row name
          rownames(cont.matrix)[i] <- paste(unique_conditions[condition2], "_vs_", unique_conditions[condition1], sep="")
        }else{
          comparison[condition1] <- 1
          comparison[condition2] <- -1
          # Add the comparison vector to the contrast matrix
          cont.matrix[i, ] <- comparison
          # Set row name
          rownames(cont.matrix)[i] <- paste(unique_conditions[condition1], "_vs_", unique_conditions[condition2], sep="")
        }
        i <- i + 1
      }
    }
  }

  cont.matrix<- t(cont.matrix)

  # Fit the linear model with contrasts
  #fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2)# Perform empirical Bayes moderation

  #### ------Extract results:
  contrast_names <- colnames(fit2$coefficients)  # Get all contrast names

  results_list <- list()# Create an empty list to store results data frames
  for (contrast_name in contrast_names) {
    # Extract results for the current contrast
    res.t <- limma::topTable(fit2, coef=contrast_name, n=Inf, sort.by="n", adjust.method = STAT_padj)%>% # coef= the comparison the test is done for!
      dplyr::rename("Log2FC"=1,
                    "t.val"=3,
                    "p.val"=4,
                    "p.adj"=5)

   res.t <- res.t%>%
      rownames_to_column("Metabolite")

    # Store the data frame in the results list, named after the contrast
    results_list[[contrast_name]] <- res.t
  }

  #If CoRe=TRUE, we need to exchange the Log2FC with the Distance and we need to combine the lists
  #Make the name_match_df
  name_match_df <- as.data.frame(names(results_list))%>%
    separate("names(results_list)", into=c("a", "b"), sep="_vs_", remove=FALSE)

  name_match_df <-merge(name_match_df, targets , by.x="a", by.y="condition_limma_compatible", all.x=TRUE)%>%
    dplyr::rename("Condition1"=5)
  name_match_df <- merge(name_match_df, targets , by.x="b", by.y="condition_limma_compatible", all.x=TRUE)%>%
    dplyr::rename("Condition2"=7)%>%
    unite("New", "Condition1", "Condition2", sep="_vs_", remove=FALSE)

  name_match_df<- name_match_df[,c(3,5)]%>%
    distinct(New, .keep_all = TRUE)

  #Match the lists using name_match_df
  for(i in 1:nrow(name_match_df)){
    old_name <- name_match_df$`names(results_list)`[i]
    new_name <- name_match_df$New[i]
    results_list[[new_name]] <- results_list[[old_name]]
    results_list[[old_name]] <- NULL
  }

  if(CoRe==TRUE){
   # Merge the data frames in list1 and list2 based on the "Metabolite" column
    merged_list <- list()
    for(i in 1:nrow(name_match_df)){
      list_dfs <- name_match_df$New[i]

      # Check if the data frames exist in both lists
      if(list_dfs %in% names(results_list) && list_dfs %in% names(Log2FC_table)){
        merged_df <- merge(results_list[[list_dfs]], Log2FC_table[[list_dfs]], by = "Metabolite", all = TRUE)
        merged_list[[list_dfs]] <- merged_df
      }
    }
    STAT_C1vC2 <- merged_list
  }else{
    STAT_C1vC2 <- results_list
  }

  if(MultipleComparison==FALSE){
    nameComp <- names(STAT_C1vC2)
    STAT_C1vC2 <-STAT_C1vC2[[nameComp]]
  }

  return(STAT_C1vC2)
}


#############################################################################################
### ### ### Shapiro function: Internal Function to perform Shapiro test and plots ### ### ###
#############################################################################################

#' @param Input_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Input_SettingsFile_Sample DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param Input_SettingsInfo \emph{Optional: } Named vector including the information about the conditions column c(conditions="ColumnName_Plot_SettingsFile"). Can additionally pass information on numerator or denominator c(numerator = "ColumnName_Plot_SettingsFile", denumerator = "ColumnName_Plot_SettingsFile") for specifying which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all). Using =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator and are presented in the right side on the Volcano plot (For CoRe the Log2Distance). \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test" or "cor.test", for one-vs-all or all-vs-all comparison choose aov (=annova), kruskal.test or lmFit (=limma) \strong{Default = "t-test"}
#' @param OutputName String which is added to the output files of the DMA.
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used \strong{Default = FALSE}
#' @param QQplots \emph {Optional: } TRUE or FALSE for whether QQ plots should be plotted  \strong{Default = TRUE}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{Default = "csv"}
#' @param Plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param Folder_Name {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @keywords Shapiro test,Normality testing, Density plot, QQplot
#' @export
#'

Shapiro <-function(Input_data,
                   Input_SettingsFile_Sample,
                   Input_SettingsInfo = c(conditions="Conditions", numerator = NULL, denumerator = NULL),
                   STAT_pval ="kruskal.test",
                   OutputName="",
                   CoRe=FALSE,
                   QQplots=TRUE,
                   Save_as_Plot="svg",
                   Save_as_Results="csv",
                   Plot=TRUE,
                   Folder_Name=NULL
){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "gtools")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }

  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(class(Input_data) != "data.frame"){
    stop("Input_data should be a data.frame. It's currently a ", paste(class(Input_data), ".",sep = ""))
  }
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Input_SettingsFile_Sample, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Input_SettingsFile_Sample"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Input_SettingsFile_Sample.")
      } else{
        Input_data <- Input_data
      }
    }
  }

  ## ------------ Check Input SettingsInfo ----------- ##
  #3. Input_SettingsInfo
  if(Input_SettingsInfo[["conditions"]] %in% Input_SettingsInfo==TRUE){
    if(Input_SettingsInfo[["conditions"]] %in% colnames(Input_SettingsFile_Sample)== FALSE){
      stop("The ",Input_SettingsInfo[["conditions"]], " column selected as Conditions in Input_SettingsInfo was not found in Input_SettingsFile_Sample. Please check your input.")
    }else{# if true rename to Conditions
      Input_SettingsFile_Sample<- Input_SettingsFile_Sample%>%
        dplyr::rename("Conditions"= paste(Input_SettingsInfo[["conditions"]]) )
    }
  }else{
    stop("You have to provide a Input_SettingsInfo for conditions.")
  }

  ##########################
  if("denominator" %in% names(Input_SettingsInfo)==TRUE){
    if(Input_SettingsInfo[["denominator"]] %in% Input_SettingsFile_Sample$Conditions==FALSE){
      stop("The ",Input_SettingsInfo[["denominator"]], " column selected as denominator in Input_SettingsInfo was not found in Input_SettingsFile_Sample. Please check your input.")
    }else{
      denominator <- Input_SettingsInfo[["denominator"]]
    }
  }
  if("numerator" %in% names(Input_SettingsInfo)==TRUE){
    if(Input_SettingsInfo[["numerator"]] %in% Input_SettingsFile_Sample$Conditions  == FALSE){
      stop("The ",Input_SettingsInfo[["numerator"]], " column selected as numerator in Input_SettingsInfo was not found in Input_SettingsFile_Sample. Please check your input.")
    }else{
      numerator <- Input_SettingsInfo[["numerator"]]
    }
  }
  if("denominator" %in% names(Input_SettingsInfo)==FALSE  & "numerator" %in% names(Input_SettingsInfo) ==TRUE){
    stop("Check input. The selected denominator option is empty while ",paste(Input_SettingsInfo[["numerator"]])," has been selected as a numerator. Please add a denminator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison." )
  }

  ## ------------ Check Denominator/numerator ----------- ##
  #4.  Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("denominator" %in% names(Input_SettingsInfo)==FALSE  & "numerator" %in% names(Input_SettingsInfo) ==FALSE){
    conditions = Input_SettingsFile_Sample$Conditions
    denominator <-unique(Input_SettingsFile_Sample$Conditions)
    numerator <-unique(Input_SettingsFile_Sample$Conditions)
  }else if("denominator" %in% names(Input_SettingsInfo)==TRUE  & "numerator" %in% names(Input_SettingsInfo)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = Input_SettingsFile_Sample$Conditions
    denominator <- Input_SettingsInfo[["denominator"]]
    numerator <-unique(Input_SettingsFile_Sample$Conditions)
  }

  ## ------------ Check General parameters ----------- ##
  #6. General parameters
  if(is.logical(Plot) == FALSE){
    stop("Check input. The plot value should be either =TRUE if a Volcano plot presenting the DMA results is to be exported or =FALSE if not.")
  }
  Save_as_Plot_options <- c("svg","pdf","png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }
  Save_as_Results_options <- c("txt","csv", "xlsx" )
  if(Save_as_Results %in% Save_as_Results_options == FALSE){
    stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
  }

  #7. Are sample numbers enough?
  Num <- Input_data %>%
    filter(Input_SettingsFile_Sample$Conditions %in% numerator) %>%
    select_if(is.numeric)#only keep numeric columns with metabolite values
  Denom <- Input_data %>%
    filter(Input_SettingsFile_Sample$Conditions %in% denominator) %>%
    select_if(is.numeric)

  if(nrow(Num)==1){
    stop("There is only one sample available for ", numerator, ", so no statistical test can be performed.")
  } else if(nrow(Denom)==1){
    stop("There is only one sample available for ", denominator, ", so no statistical test can be performed.")
  }else if(nrow(Num)==0){
    stop("There is no sample available for ", numerator, ".")
  }else if(nrow(Denom)==0){
    stop("There is no sample available for ", denominator, ".")
  }


  ## ------------ Create Results output folder ----------- ##
  if(is.null(Folder_Name)==TRUE){
    WorkD <- getwd()
    Results_folder <- file.path(WorkD, name)
    if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
    Results_folder_DMA_folder <- file.path(Results_folder,"DMA") # Make DMA results folder
    if (!dir.exists(Results_folder_DMA_folder)) {dir.create(Results_folder_DMA_folder)}
    Results_folder_DMA_folder_Shapiro_folder <- file.path(Results_folder_DMA_folder,"Shapiro") # Make DMA results folder
    if (!dir.exists(Results_folder_DMA_folder_Shapiro_folder)) {dir.create(Results_folder_DMA_folder_Shapiro_folder)}

  }else{
    Results_folder_DMA_folder_Shapiro_folder <- Folder_Name
  }

  ###############################################################################################################################################################################################################
  ## ------------ Check data normality and statistical test chosen and generate Output DF----------- ##
  # Before Hypothesis testing, we have to decide whether to use a parametric or a non parametric test. We can test the data normality using the Shapiro test.
  ##-------- First: Load the data and perform the shapiro.test on each metabolite across the samples of one condition. this needs to be repeated for each condition:
  #Prepare the input:
  Input_shaptest <- replace(Input_data, Input_data==0, NA) %>% #Shapiro test ignores NAs!
    filter(Input_SettingsFile_Sample$Conditions %in% numerator | Input_SettingsFile_Sample$Conditions %in% denominator)%>%
    select_if(is.numeric)
  temp<- as.vector(sapply(Input_shaptest, function(x) var(x)) == 0)#  we have to remove features with zero variance if there are any.
  Input_shaptest <- Input_shaptest[,!temp]
  Input_shaptest_Cond <-merge(data.frame(Conditions = Input_SettingsFile_Sample[, "Conditions", drop = FALSE]), Input_shaptest, by=0, all.y=TRUE)

  UniqueConditions <- Input_SettingsFile_Sample%>%
    subset(Input_SettingsFile_Sample$Conditions %in% numerator | Input_SettingsFile_Sample$Conditions %in% denominator, select = c(1))
  UniqueConditions <- unique(UniqueConditions$Conditions)

  #Generate the results
  shapiro_results <- list()
  for (i in UniqueConditions) {
    # Subset the data for the current condition
    subset_data <- Input_shaptest_Cond%>%
      column_to_rownames("Row.names")%>%
      subset(Conditions == i, select = -c(1))

    # Apply Shapiro-Wilk test to each feature in the subset
    shapiro_results[[i]] <- as.data.frame(sapply(subset_data, function(x) shapiro.test(x)))
  }

  #Make the output DF
  DF_shapiro_results <- as.data.frame(matrix(NA, nrow = length(UniqueConditions), ncol = ncol(Input_shaptest)))
  rownames(DF_shapiro_results) <- UniqueConditions
  colnames(DF_shapiro_results) <- colnames(Input_shaptest)
  for(k in 1:length(UniqueConditions)){
    for(l in 1:ncol(Input_shaptest)){
      DF_shapiro_results[k, l] <- shapiro_results[[UniqueConditions[k]]][[l]]$p.value
    }
  }
  colnames(DF_shapiro_results) <- paste("Shapiro p.val(", colnames(DF_shapiro_results),")", sep = "")

  ##------ Second: Give feedback to the user if the chosen test fits the data distribution. The data are normal if the p-value of the shapiro.test > 0.05.
  Density_plots <- list()
  if(QQplots==TRUE){
    QQ_plots <- list()
  }
  for(x in 1:nrow(DF_shapiro_results)){
    transpose <- as.data.frame(t(DF_shapiro_results[x,]))
    Norm <- format((round(sum(transpose[[1]] > 0.05)/nrow(transpose),4))*100, nsmall = 2) # Percentage of normally distributed metabolites across samples
    NotNorm <- format((round(sum(transpose[[1]] < 0.05)/nrow(transpose),4))*100, nsmall = 2) # Percentage of not-normally distributed metabolites across samples
    if(STAT_pval =="kruskal.test" | STAT_pval =="wilcox.test"){
      message("For the condition ", colnames(transpose) ," ", Norm, " % of the metabolites follow a normal distribution and ", NotNorm, " % of the metabolites are not-normally distributed according to the shapiro test. You have chosen ",paste(STAT_pval), ", which is for non parametric Hypothesis testing. `shapiro.test` ignores missing values in the calculation.")
    }else{
      message("For the condition ", colnames(transpose) ," ", Norm, " % of the metabolites follow a normal distribution and ", NotNorm, " % of the metabolites are not-normally distributed according to the shapiro test. You have chosen ",paste(STAT_pval), ", which is for parametric Hypothesis testing. `shapiro.test` ignores missing values in the calculation.")
    }

    # Assign the calculated values to the corresponding rows in result_df
    DF_shapiro_results$`Metabolites with normal distribution [%]`[x] <- Norm
    DF_shapiro_results$`Metabolites with not-normal distribution [%]`[x] <- NotNorm

    #reorder the DF:
    DF_shapiro_results<-DF_shapiro_results[,c(ncol(DF_shapiro_results)-1, ncol(DF_shapiro_results), 1:(ncol(DF_shapiro_results)-2))]

    DF_shapiro_results_out<- t(DF_shapiro_results)%>% as.data.frame()%>% rownames_to_column("Shapiro_p.val")
    DF_shapiro_results_out$Shapiro_p.val <-  str_replace_all(DF_shapiro_results_out$Shapiro_p.val, "Shapiro p.val", " ")
    DF_shapiro_results_out$Shapiro_p.val <-gsub("[[:punct:]]", " ", DF_shapiro_results_out$Shapiro_p.val)

    # Save the DF Shapiro
    if (Save_as_Results == "xlsx"){
      writexl::write_xlsx(DF_shapiro_results_out,paste(Results_folder_DMA_folder_Shapiro_folder,"/DF_shapiro_results_table",OutputName,".",Save_as_Results,sep =  "")) # save the DMA result DF
    }else if (Save_as_Results == "csv"){
      write.csv(DF_shapiro_results_out,paste(Results_folder_DMA_folder_Shapiro_folder,"/DF_shapiro_results_table",OutputName,".",Save_as_Results,sep =  ""),row.names =FALSE) # save the DMA result DF
    }else if (Save_as_Results == "txt"){
      write.table(DF_shapiro_results_out,paste(Results_folder_DMA_folder_Shapiro_folder,"/DF_shapiro_results_table",OutputName,".",Save_as_Results,sep =  ""), col.names = TRUE, row.names = FALSE) # save the DMA result DF
    }

    ## Make Group wise data distribution plot and QQ plots
    subset_data <- Input_shaptest_Cond%>%
      column_to_rownames("Row.names")%>%
      subset(Conditions ==  colnames(transpose), select = -c(1))
    all_data <- unlist(subset_data)

    plot <- ggplot(data.frame(x = all_data), aes(x = x)) +
      geom_histogram(aes(y=..density..), binwidth=.5, colour="black", fill="white")  +
      geom_density(alpha = 0.2, fill = "grey45")

    density_values <- ggplot_build(plot)$data[[2]]

    plot <- ggplot(data.frame(x = all_data), aes(x = x)) +
      geom_histogram(aes(y=..density..), binwidth=.5, colour="black", fill="white") +
      geom_density(alpha=.2, fill="grey45") +
      scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))]))

    density_values2 <- ggplot_build(plot)$data[[2]]

    suppressWarnings( sampleDist <- ggplot(data.frame(x = all_data), aes(x = x)) +
                        geom_histogram(aes(y=..density..), binwidth=.5, colour="black", fill="white") +
                        geom_density(alpha=.2, fill="grey45") +
                        scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))])) +
                        theme_minimal()+
                        # geom_vline(xintercept =median(all_data) , linetype = "dashed", color = "red")+
                        labs(title=paste("Data distribution ",  colnames(transpose)), subtitle = paste(NotNorm, " of metabolites not normally distributed based on Shapiro test"),x="Abundance", y = "Density")#+
                      # geom_text(aes(x = density_values2$x[which.max(density_values2$y)], y = 0, label = "Median"),  vjust = 0, hjust = -0.5, color = "red", size = 3.5)  # Add label for
    )

    plot(sampleDist)
    Density_plots[[paste(colnames(transpose))]] <- recordPlot()

    if(CoRe==TRUE){
      ggsave(filename = paste0(Results_folder_DMA_folder_Shapiro_folder, "/Density_plot", paste(colnames(transpose)),OutputName,".",Save_as_Plot), plot = sampleDist, width = 10,  height = 8)
    }else{
      ggsave(filename = paste0(Results_folder_DMA_folder_Shapiro_folder, "/Density_plot", paste(colnames(transpose)),OutputName,".",Save_as_Plot), plot = sampleDist, width = 10,  height = 8)

    }
    # QQ plots
    if(QQplots==TRUE){
      # Make folders !has to be moved on top!
      conds <- unique(c(numerator, denominator))
      for(x in conds){
        Results_folder_DMA_folder_Shapiro_folder_Condition <- file.path(Results_folder_DMA_folder_Shapiro_folder, paste(x)) # Make DMA results folder
        if (!dir.exists(Results_folder_DMA_folder_Shapiro_folder_Condition)) {dir.create(Results_folder_DMA_folder_Shapiro_folder_Condition)}
      }
      #QQ plots for each groups for each metabolite for normality visual check
      qq_plot_list <- list()
      for (col_name in colnames(subset_data)){
        qq_plot <- ggplot(data.frame(x = subset_data[[col_name]]), aes(sample = x)) +
          geom_qq() +
          geom_qq_line(color = "red") +
          labs(title = paste("QQPlot for", col_name),x = "Theoretical", y="Sample")+ theme_minimal()

        plot.new()
        plot(qq_plot)
        qq_plot_list[[col_name]] <-  recordPlot()

        col_name2 <- (gsub("/","_",col_name))#remove "/" cause this can not be safed in a PDF name
        col_name2 <- gsub("-", "", col_name2)
        col_name2 <- gsub("/", "", col_name2)
        col_name2 <- gsub(" ", "", col_name2)
        col_name2 <- gsub("\\*", "", col_name2)
        col_name2 <- gsub("\\+", "", col_name2)
        col_name2 <- gsub(",", "", col_name2)
        col_name2 <- gsub("\\(", "", col_name2)
        col_name2 <- gsub("\\)", "", col_name2)

        ggsave(paste0(Results_folder_DMA_folder_Shapiro_folder, "/", paste(colnames(transpose)),"/",paste(col_name2),".",Save_as_Plot), plot = qq_plot, device = Save_as_Plot, width = 10,  height = 8)

        dev.off()
      }

      QQ_plots[[paste(colnames(transpose))]] <- qq_plot_list
    }
  }

  #Here we make a list in which we will save the output
  if(QQplots==TRUE){
    Shapiro_output_list <- list("DF" = list("Shapiro_result"=DF_shapiro_results),"Plot"=list( "Distributions"=Density_plots, "QQ_plots" = QQ_plots))
  }else{
    Shapiro_output_list <- list("DF" = list("Shapiro_result"=DF_shapiro_results),"Plot"=list( "Distributions"=Density_plots))
  }
  suppressWarnings(invisible(return(Shapiro_output_list)))

  if(Plot == TRUE){
    Shapiro_output$Plot$Distributions
  }
}
