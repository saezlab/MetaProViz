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
#' @param Input_SettingsFile DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param Input_SettingsInfo \emph{Optional: } Named vector including the information about the conditions column c(conditions="ColumnName_Plot_SettingsFile"). Can additionally pass information on numerator or denominator c(numerator = "ColumnName_Plot_SettingsFile", denumerator = "ColumnName_Plot_SettingsFile") for specifying which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all). Using =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator and are presented in the right side on the Volcano plot (For CoRe the Log2Distance). \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test" or "cor.test", for one-vs-all or all-vs-all comparison choose aov (=annova), kruskal.test or lmFit (=limma) \strong{Default = "t-test"}
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param OutputName String which is added to the output files of the DMA.
#' @param Input_MetaFile_Metab \emph{Optional: } DF which contains the metadata information , i.e. pathway information, retention time,..., for each metabolite. \strong{Default = NULL}
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
               Input_SettingsFile,
               Input_SettingsInfo = c(conditions="Conditions", numerator = NULL, denumerator = NULL),
               STAT_pval ="kruskal.test",
               STAT_padj="fdr",
               Input_MetaFile_Metab = NULL,
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
      Test_match <- merge(Input_SettingsFile, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Input_SettingsFile"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Input_SettingsFile.")
      } else{
        Input_data <- Input_data
      }
    }
  }

  #2.  Input_MetaFile_Metab
  if(is.null(Input_MetaFile_Metab) == FALSE){
    if('Metabolite' %in% colnames(Input_MetaFile_Metab) == FALSE){
      warning("The provided file Input_MetaFile_Metab must have a columns named: `Metabolite`.")
    }
  }

  ## ------------ Check Input SettingsInfo ----------- ##
  #3. Input_SettingsInfo
  if(Input_SettingsInfo[["conditions"]] %in% Input_SettingsInfo==TRUE){
    if(Input_SettingsInfo[["conditions"]] %in% colnames(Input_SettingsFile)== FALSE){
      stop("The ",Input_SettingsInfo[["conditions"]], " column selected as Conditions in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
    }else{# if true rename to Conditions
      Input_SettingsFile<- Input_SettingsFile%>%
        dplyr::rename("Conditions"= paste(Input_SettingsInfo[["conditions"]]) )
    }
  }else{
    stop("You have to provide a Input_SettingsInfo for conditions.")
  }

  ##########################
  if("denominator" %in% names(Input_SettingsInfo)==TRUE){
    if(Input_SettingsInfo[["denominator"]] %in% Input_SettingsFile$Conditions==FALSE){
      stop("The ",Input_SettingsInfo[["denominator"]], " column selected as denominator in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
    }else{
      denominator <- Input_SettingsInfo[["denominator"]]
    }
  }
  if("numerator" %in% names(Input_SettingsInfo)==TRUE){
    if(Input_SettingsInfo[["numerator"]] %in% Input_SettingsFile$Conditions  == FALSE){
      stop("The ",Input_SettingsInfo[["numerator"]], " column selected as numerator in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
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
    conditions = Input_SettingsFile$Conditions
    denominator <-unique(Input_SettingsFile$Conditions)
    numerator <-unique(Input_SettingsFile$Conditions)
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("denominator" %in% names(Input_SettingsInfo)==TRUE  & "numerator" %in% names(Input_SettingsInfo)==FALSE){
    #all-vs-one: Generate the pairwise combinations


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
    STAT_pval_options <- c("t.test", "wilcox.test","chisq.test", "cor.test")
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
    stop("Check input. The selected STAT_padj option for multiple Hypothesis testing correction is not valid. Please select one of the folowwing: ",paste(STAT_padj_options,collapse = ", "),"." )
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
    filter(Input_SettingsFile$Conditions %in% numerator) %>%
    select_if(is.numeric)#only keep numeric columns with metabolite values
  Denom <- Input_data %>%
    filter(Input_SettingsFile$Conditions %in% denominator) %>%
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
    name <- paste("MetaProViz_Results",Sys.Date(),Folder_Name,sep = "_" )
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


  ###############################################################################################################################################################################################################
  ## ------------ Check data normality and statistical test chosen and generate Output DF----------- ##
  # Before Hypothesis testing, we have to decide whether to use a parametric or a non parametric test. We can test the data normality using the Shapiro test.
  ##-------- First: Load the data and perform the shapiro.test on each metabolite across the samples of one condition. this needs to be repeated for each condition:
  #Prepare the input:
  Input_shaptest <- replace(Input_data, Input_data==0, NA) %>% #Shapiro test ignores NAs!
    filter(Input_SettingsFile$Conditions %in% numerator | Input_SettingsFile$Conditions %in% denominator)%>%
    select_if(is.numeric)
  temp<- as.vector(sapply(Input_shaptest, function(x) var(x)) == 0)#  we have to remove features with zero variance if there are any.
  Input_shaptest <- Input_shaptest[,!temp]
  Input_shaptest_Cond <-merge(data.frame(Conditions = Input_SettingsFile[, "Conditions", drop = FALSE]), Input_shaptest, by=0, all.y=TRUE)

  UniqueConditions <- Input_SettingsFile%>%
    subset(Input_SettingsFile$Conditions %in% numerator | Input_SettingsFile$Conditions %in% denominator, select = c(1))
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
  QQ_plots <- list()
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
    # # QQ plots
    # # Make folders !has to be moved on top!
    # conds <- unique(c(numerator, denominator))
    # for(x in conds){
    #   Results_folder_DMA_folder_Shapiro_folder_Condition <- file.path(Results_folder_DMA_folder_Shapiro_folder, paste(x)) # Make DMA results folder
    #   if (!dir.exists(Results_folder_DMA_folder_Shapiro_folder_Condition)) {dir.create(Results_folder_DMA_folder_Shapiro_folder_Condition)}
    # }
    # QQ plots for each groups for each metabolite for normality visual check
    # qq_plot_list <- list()
    # for (col_name in colnames(subset_data)){
    #   qq_plot <- ggplot(data.frame(x = subset_data[[col_name]]), aes(sample = x)) +
    #     geom_qq() +
    #     geom_qq_line(color = "red") +
    #     labs(title = paste("QQPlot for", col_name),x = "Theoretical", y="Sample")+ theme_minimal()
    #
    #   plot.new()
    #   plot(qq_plot)
    #   qq_plot_list[[col_name]] <-  recordPlot()
    #
    #   col_name2 <- (gsub("/","_",col_name))#remove "/" cause this can not be safed in a PDF name
    #   col_name2 <- gsub("-", "", col_name2)
    #   col_name2 <- gsub("/", "", col_name2)
    #   col_name2 <- gsub(" ", "", col_name2)
    #   col_name2 <- gsub("\\*", "", col_name2)
    #   col_name2 <- gsub("\\+", "", col_name2)
    #   col_name2 <- gsub(",", "", col_name2)
    #   col_name2 <- gsub("\\(", "", col_name2)
    #   col_name2 <- gsub("\\)", "", col_name2)
    #
    #   ggsave(paste0(Results_folder_DMA_folder_Shapiro_folder, "/", paste(colnames(transpose)),"/",paste(col_name2),".",Save_as_Plot), plot = qq_plot, device = Save_as_Plot, width = 10,  height = 8)
    #
    #   dev.off()
    #}

    #QQ_plots[[paste(colnames(transpose))]] <- qq_plot_list
  }


   ###############################################################################################################################################################################################################
  #### Prepare the data ######
  #1. Metabolite names:
  # If anova is used change the column names
  if(MultipleComparison == TRUE){
    for (i in 1:length(colnames(Input_data))){
      metabolite_name <- colnames(Input_data)[i]
      metabolite_name <- gsub("-", "", metabolite_name)
      metabolite_name <- gsub("/", "", metabolite_name)
      metabolite_name <- gsub(" ", "", metabolite_name)
      metabolite_name <- gsub("\\*", "", metabolite_name)
      metabolite_name <- gsub("\\+", "", metabolite_name)
      metabolite_name <- gsub(",", "", metabolite_name)
      metabolite_name <- gsub("\\(", "", metabolite_name)
      metabolite_name <- gsub("\\)", "", metabolite_name)
      colnames(Input_data)[i] <- metabolite_name
    }
  }



  Log2FC_table <- data.frame(Metabolite = colnames(Input_data))
  for (column in 1:dim(comparisons)[2]){

    C1 <- Input_data %>% # Numerator
      filter(Input_SettingsFile$Conditions %in% comparisons[2,column]) %>%
      select_if(is.numeric)#only keep numeric columns with metabolite values
    C2 <- Input_data %>% # Deniminator
      filter(Input_SettingsFile$Conditions %in%  comparisons[1,column]) %>%
      select_if(is.numeric)
  }


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


     Log2FC_table <-Mean_Merge[,c(1,5)]
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
      Log2FC_C1vC2 <-Mean_Merge[,c(1,8)]

      if(MultipleComparison == TRUE){
        logname <- paste("Log2FC", paste( comparisons[2,column], comparisons[1,column],sep="-"), sep="_")
        names(Log2FC_C1vC2)[2] <- logname
      }
      Log2FC_table <- merge(Log2FC_table, Log2FC_C1vC2, by= "Metabolite")


    }else{
      stop("Please choose CoRe= TRUE or CoRe=FALSE.")
    }


  ## ------------ Perform Hypothesis testing ----------- ##
  if(MultipleComparison == FALSE){
    STAT_C1vC2 <-MetaProViz:::DMA_Stat_single(C1=C1, C2=C2, Log2FC_table=Log2FC_table, Metabolites_Miss=Metabolites_Miss, STAT_pval=STAT_pval, STAT_padj=STAT_padj)

  }else{ # ANOVA = TRUE

    # for 1 vs all
    conditions =as.factor(conditions)
    #conditions=relevel(conditions, ref = "HK2")

    if(STAT_pval=="anova"){
      ## 1. Anova
      aov.res= apply(Input_data,2,function(x) aov(x~conditions))
      anova.res = lapply(aov.res  , anova)
      anova_res<-do.call('rbind', lapply(anova.res, function(x) {x["Pr(>F)"][1,]}))
      anova_res = cbind(anova_res,p.adjust(anova_res,STAT_padj )) %>% as.data.frame()
      colnames(anova_res) = c('Anova_p.val','Anova_p.adj')
      anova_res$Anova_p.val <- NULL # remove p.val

      ## 2. Tukey test
      posthoc.res = lapply(aov.res, TukeyHSD, conf.level=0.95)
      Tukey_res <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'p adj']))
      Tukey_res <- as.data.frame(Tukey_res)
      names(Tukey_res) <- paste("p.adj", names(Tukey_res), sep = "_")

      Pval_table <- merge(anova_res, Tukey_res, by=0)
      names(Pval_table)[1] <- "Metabolite"

      ### put Log2FC and pval dataframes
      merged_table <- merge(Pval_table, Log2FC_table, by = "Metabolite", all = TRUE)

    }else{# STAT_pval = kruskal.test
      # Kruskal test
      aov.res= apply(Input_data,2,function(x) kruskal.test(x~conditions))
      anova_res<-do.call('rbind', lapply(aov.res, function(x) {x["p.value"]}))
      anova_res <- as.matrix(mutate_all(as.data.frame(anova_res), function(x) as.numeric(as.character(x))))
      colnames(anova_res) = c("Kruskal_p.val")
      # Dunn test
      Dunndata <- Input_data %>%
        mutate(conditions = conditions) %>%
        select(conditions, everything())%>% as.data.frame()
      # apply doesnt for so going for a loop
      Dunn_res<- data.frame( comparisons = paste(comparisons[2,],    comparisons[1,], sep = "_" ))
      for(col in 2:dim(Dunndata)[2]){
        data = Dunndata[,c(1,col)]
        colnames(data)[2] <- gsub("^\\d+", "", colnames(data)[2])

        ## If a metabolite starts with number remove it
        formula <- as.formula(paste(colnames(data)[2], "~ conditions"))
        posthoc.res= rstatix::dunn_test(data, formula, p.adjust.method = STAT_padj)

        res <- data.frame(comparisons =  paste(posthoc.res$group2, posthoc.res$group1, sep = "_" ))
        res[[colnames(Dunndata)[col] ]] <-  posthoc.res$p.adj
        Dunn_res <- merge(Dunn_res,res,by="comparisons")
      }
      Dunn_res <- column_to_rownames(Dunn_res, "comparisons")%>% t() %>% as.data.frame()
      colnames(Dunn_res) <- paste("p.adj", colnames(Dunn_res), sep = "_")

      Dunn_res <- as.matrix(mutate_all(as.data.frame(Dunn_res), function(x) as.numeric(as.character(x))))

      Pval_table <- merge(anova_res, Dunn_res, by=0)
      names(Pval_table)[1] <- "Metabolite"

      merged_table <- merge(Pval_table, Log2FC_table, by = "Metabolite", all = TRUE)


    }
    #### 3. t.value
    for (i in 1:(dim(Pval_table)[2]-2)){
      t.val_name <- gsub("p.adj", "t.val", colnames(merged_table[i+2]) )
      merged_table[[t.val_name]] <- qnorm((1 - merged_table[[i+2]] / 2)) * sign(merged_table[[i+2+dim(Pval_table)[2]-2]]) # calculate and add t-value
    }

    # re-order columns
    order<- c(1,2)
    for(i in 1: (dim(Pval_table)[2]-2)){
      order <- append(order, dim(Pval_table)[2]+i)
      order <- append(order, 2+i)
      order <- append(order, (dim(Pval_table)[2] +dim(Log2FC_table)[2]-1) +i)
    }
    merged_table <- merged_table[, order]
    STAT_C1vC2 <- merged_table

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
        CoRe_info[3,i] <- paste("Released in" ,numerator , "and Consumed",denominator , sep=" ")
      } else if(temp1[i]<0 & temp2[i]>0){
        CoRe_info[3,i] <- paste("Consumed in" ,numerator , " and Released",denominator , sep=" ")
      }else{
        CoRe_info[3,i] <- "No Change"
      }
    }

    CoRe_info <- t(CoRe_info) %>% as.data.frame()
    CoRe_info <- rownames_to_column(CoRe_info, "Metabolite")
    names(CoRe_info)[2] <- paste("Mean", "CoRe", numerator, sep="_")
    names(CoRe_info)[3] <- paste("Mean", "CoRe", denominator, sep="_")
    names(CoRe_info)[4] <- "CoRe"
    STAT_C1vC2 <- merge(STAT_C1vC2,CoRe_info,by= "Metabolite")
    STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value
  }

  ## ------------ Add pathway information to DMA results ----------- ##
  if(is.null(Input_MetaFile_Metab)!=TRUE & 'Metabolite' %in% colnames(Input_MetaFile_Metab)){
    STAT_C1vC2<- merge(STAT_C1vC2, Input_MetaFile_Metab,by="Metabolite", all.x=T)
  }
  DMA_Output <- STAT_C1vC2

  if(MultipleComparison==TRUE){
    if(All==TRUE){
      numerator <- "Anova_All"
      denominator<- "All"
    }
  }


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

  # Make a simple Volcano plot
  dev.new()
  if(CoRe==TRUE){
    x <- "Log2(Distance)"
    VolPlot_SettingsInfo= c(color="CoRe")
    VOlPlot_SettingsFile = DMA_Output

  }else{
    x <- "Log2FC"
    VolPlot_SettingsInfo= NULL
    VOlPlot_SettingsFile = NULL
  }
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
                                                  SelectLab= DMA_Output$Metabolite,
                                                  Connectors=  FALSE,
                                                  Subtitle=  bquote(italic("Differential Metabolite Analysis")),
                                                  Theme= NULL,
                                                  Save_as_Plot= NULL))
    dev.off()
    plot(VolcanoPlot)

    volcanoDMA <- file.path(Results_folder_Conditions,paste0( "Volcano_Plot_",toString(numerator),"_versus_",toString(denominator),OutputName,".",Save_as_Plot))
    ggsave(volcanoDMA,plot=VolcanoPlot, width=10, height=8) # save the volcano plot

    output_list <- list()  #Here we make a list in which we will save the output
    DMA_output_list <- list("DF" = list("Shapiro_result"=DF_shapiro_results,"DMA_result"=DMA_Output),"Plot"=list( "Distributions"=Density_plots, "QQ_plots" = QQ_plots, "Volcano"=VolcanoPlot))



    if(Plot == TRUE){
      VolcanoPlot
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
#message("No conditions were specified as numerator or denumerator. Performing multiple testing `all-vs-all` using", paste(STAT_pval), ".")

#all-vs-one:
#message("No conditions were specified as numerator. Performing multiple testing `one-vs-all` using", paste(STAT_pval), ".")

# one-vs-one:
#message("conditions were specified as numerator and denumerator. Performing testing `one-vs-one` using", paste(STAT_pval), ".")



