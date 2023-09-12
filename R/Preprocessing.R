## ---------------------------
##
## Script name: Preprocessing
##
## Purpose of script: Metabolomics (raw ion counts) pre-processing, normalisation and outlier detection
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


#' Applies 80%-filtering rule, total-ion count normalisation, missing value imputation and HotellingT2 outlier detection
#'
#' @param Input_data DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Input_SettingsFile DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param Input_SettingsInfo  NULL or Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "BiologicalReplicates" including numerical values. For CoRe = TRUE a CoRe_norm_factor = "Columnname_Input_SettingsFile" and CoRe_media = "Columnname_Input_SettingsFile", have to also be added. Column CoRe_norm_factor is used for normalization and CoRe_media is used to specify the name of the media controls in the Conditions.
#' @param Feature_Filtering \emph{Optional: }If set to "None" then no feature filtering is performed. If set to "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = Modified}
#' @param Feature_Filt_Value \emph{Optional: } Percentage of feature filtering. \strong{Default = 0.8}
#' @param TIC_Normalization \emph{Optional: } If TRUE, Total Ion Count normalization is performed. \strong{Default = TRUE}
#' @param MVI \emph{Optional: } If TRUE, Missing Value Imputation (MVI) based on half minimum is performed \strong{Default = TRUE}
#' @param MVI_Percentage \emph(Optional: ) Percentage (0 to 100)of imputed value based on the minimum value. \strong{Default = 50}
#' @param HotellinsConfidence \emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. Must be numeric.\strong{Default = 0.99}
#' @param ExportQCPlots \emph{Optional: } Select whether the quality control (QC) plots will be exported. \strong{Default = TRUE}
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed and the CoRe value will be calculated. Please consider providing a Normalisation factor column called "CoRe_norm_factor" in your "Input_SettingsFile" DF, where the column "Conditions" matches. The normalisation factor must be a numerical value obtained from growth rate that has been obtained from a growth curve or growth factor that was obtained by the ratio of cell count/protein quantification at the start point to cell count/protein quantification at the end point.. Additionally control media samples have to be available in the "Input" DF and defined as "CoRe_media" samples in the "Conditions" column in the "Input_SettingsFile" DF. \strong{Default = FALSE}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#'
#' @keywords 80% filtering rule, Missing Value Imputation, Total Ion Count normalization, PCA, HotellingT2, multivariate quality control charts,
#' @export


###################################################
### ### ### Metabolomics pre-processing ### ### ###
###################################################


Preprocessing <- function(Input_data,
                          Input_SettingsFile,
                          Input_SettingsInfo,
                          Feature_Filtering = "Modified",
                          Feature_Filt_Value = 0.8,
                          TIC_Normalization = TRUE,
                          MVI= TRUE,
                          MVI_Percentage=50,
                          HotellinsConfidence = 0.99,
                          ExportQCPlots = TRUE,
                          CoRe = FALSE,
                          Save_as_Plot = "svg"
){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", # general scripting
                        "factoextra", # visualize PCA
                        "qcc", # for hotelling plots
                        "ggplot2", # For visualization PCA
                        "hash", # Dictionary in R for making column of outliers
                        "reshape", # for melting df for anova
                        "patchwork", # for combining ggplots
                        "inflection"
                        )# For finding inflection point/ Elbow knee /PCA component selection # https://cran.r-project.org/web/packages/inflection/inflection.pdf # https://deliverypdf.ssrn.com/delivery.php?ID = 454026098004123081018105104090015093000085002012023032095093077109069092095000114006057018122039107109012089110120018031068078025094036037013095100070100076109026029024044005068010070117123085122016083112098002109001027028000024115096122101001083084026&EXT = pdf&INDEX = TRUE # https://arxiv.org/abs/1206.5478

  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------------ Run ------------------- ##

  #######################################################################################
  ### ### ### Check Input Information and add Input_SettingsFile information ### ### ###

  #1.  Input data
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Input_SettingsFile, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Input_SettingsFile"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Input_SettingsFile.")
      } else(
        Input_data <- Input_data
      )
    }
  }

  #2. The Input_settings: Input_SettingInfo and Input_SettingFile
  if(is.vector(Input_SettingsInfo)==TRUE & is.null(Input_SettingsFile)==TRUE){
    stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile.")
  }
  if(is.null(Input_SettingsInfo)==TRUE & is.null(Input_SettingsFile)==FALSE){
    warning("You have added a Plot_SettingsFile DF but the Plot_SettingsInfo option is empty. If you want to preprocess based some experimental condition please specify it in the Input_SettingsInfo.")
  }
  if(is.null(Input_SettingsInfo)==TRUE & is.null(Input_SettingsFile)==TRUE){
    message("No Input_Settings have been added.")
  }

  if(is.vector(Input_SettingsInfo)==TRUE){
    if ( "Conditions" %in% names(Input_SettingsInfo)){
      if(Input_SettingsInfo[["Conditions"]] %in% colnames(Input_SettingsFile)== FALSE){
        stop("The ",Input_SettingsInfo[["Conditions"]], " column selected as Conditions in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
      }else{# if true rename to Conditions
        Input_SettingsFile<- Input_SettingsFile%>%
          dplyr::rename("Conditions"= paste(Input_SettingsInfo[["Conditions"]]) )
      }
    }

    if ( "Biological_Replicates" %in% names(Input_SettingsInfo)){
      if(Input_SettingsInfo[["Biological_Replicates"]] %in% colnames(Input_SettingsFile)== FALSE){
        stop("The ",Input_SettingsInfo[["Biological_Replicates"]], " column selected as Biological_Replicates in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
      }else{# if true rename to Conditions
        Input_SettingsFile<- Input_SettingsFile%>%
          dplyr::rename("Biological_Replicates"= paste(Input_SettingsInfo[["Biological_Replicates"]]) )
      }
    }
  }

  #3. Core parameters
  if (CoRe ==  TRUE){   # parse CoRe normalisation factor
    message("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
    if ("CoRe_media" %in% names(Input_SettingsInfo)){
      Input_SettingsFile <- Input_SettingsFile %>%
        mutate(Conditions = ifelse(Conditions == paste(Input_SettingsInfo[["CoRe_media"]]), "CoRe_media", Conditions))
    }

    if(length(grep("CoRe_media", Input_SettingsFile$Conditions)) < 1){     # Check for CoRe_media samples
      stop("No CoRe_media samples were provided in the 'Conditions' in the Experimental design'. For a CoRe experiment control media samples without cells have to be measured and be added in the 'Conditions'
           column labeled as 'CoRe_media' (see @param section). Please make sure that you used the correct labelling or whether you need CoRE = FALSE for your analysis")
    }
    if ("CoRe_norm_factor" %in% names(Input_SettingsInfo)){
      Input_SettingsFile<- Input_SettingsFile%>%
        dplyr::rename("CoRe_norm_factor"= paste(Input_SettingsInfo[["CoRe_norm_factor"]]) )
      CoRe_norm_factor <-   Input_SettingsFile %>% filter(Conditions!="CoRe_media") %>% select(CoRe_norm_factor) %>%pull()

    }else{
      warning("No growth rate or growth factor provided for normalising the CoRe result, hence CoRe_norm_factor set to 1 for each sample")
      CoRe_norm_factor <- as.numeric(rep(1,dim(Input_SettingsFile %>% filter(Conditions!="CoRe_media"))[1]))
    }
  }

  #4. Parse Conditions
  if ( "Conditions" %in% colnames(Input_SettingsFile)){   # Parse Condition and Replicate information
    Conditions <- Input_SettingsFile$Conditions
  }else{
    Conditions <- NULL
  }

  #5. General parameters
  Feature_Filtering_options <- c("Standard","Modified", "none")
  if(Feature_Filtering %in% Feature_Filtering_options == FALSE ){
    stop("Check input. The selected Feature_Filtering option is not valid. Please select one of the folowwing: ",paste(Feature_Filtering_options,collapse = ", "),"." )
  }
  if( is.numeric(Feature_Filt_Value) == FALSE |Feature_Filt_Value > 1 | Feature_Filt_Value < 0){
    stop("Check input. The selected Filtering value should be numeric and between 0 and 1.")
  }
  if(is.logical(TIC_Normalization) == FALSE){
    stop("Check input. The TIC_Normalization value should be either =TRUE if TIC normalization is to be performed or =FALSE if no data normalization is to be applied.")
  }
  if(is.logical(MVI) == FALSE){
    stop("Check input. MVI value should be either =TRUE if mising value imputation should be performed or =FALSE if not.")
  }
  if( is.numeric(MVI_Percentage)== FALSE |HotellinsConfidence > 100 | HotellinsConfidence < 0){
    stop("Check input. The selected MVI_Percentage value should be numeric and between 0 and 100.")
  }
  if( is.numeric(HotellinsConfidence)== FALSE |HotellinsConfidence > 1 | HotellinsConfidence < 0){
    stop("Check input. The selected Filtering value should be numeric and between 0 and 1.")
  }
  if(is.logical(ExportQCPlots) == FALSE){
    stop("Check input. The ExportQCPlots value should be either =TRUE if QC plots are to be exported or =FALSE if not.")
  }
  if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for preprocessing of Consuption/Release experiment or =FALSE if not.")
  }
  Save_as_Plot_options <- c("svg","pdf", "png")
  if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
    stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
  }

  Input_data <- as.matrix(mutate_all(as.data.frame(Input_data), function(x) as.numeric(as.character(x))))


  #############################################
  ### ### ### Create output folders ### ### ###

  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_Preprocessing_folder = file.path(Results_folder, "Preprocessing")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  if (!dir.exists(Results_folder_Preprocessing_folder)) {dir.create(Results_folder_Preprocessing_folder)}  # check and create folder
  Results_folder_Preprocessing_Outlier_detection_folder = file.path(Results_folder_Preprocessing_folder, "Outlier_detection")   # Create Outlier_Detection directory
  if (!dir.exists(Results_folder_Preprocessing_Outlier_detection_folder)) {dir.create(Results_folder_Preprocessing_Outlier_detection_folder)}
  if (ExportQCPlots ==  TRUE){   # Create Quality_Control_PCA directory
    Results_folder_Preprocessing_folder_Quality_Control_PCA_folder = file.path(Results_folder_Preprocessing_folder, "Quality_Control_PCA")
    if (!dir.exists(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder)) {dir.create(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder)}
  }


  #########################################
  ### ### ### Feature filtering ### ### ###
  #message("Feature filtering is performed to reduce missing values that can bias the analysis and cause methods to underperform, which leads to low precision in the statistical analysis. REF: Steuer et. al. (2007), Methods Mol Biol. 358:105-26., doi:10.1007/978-1-59745-244-1_7.")
  #Prepare data:
  Input_data <- replace(Input_data, Input_data==0, NA)
  Original_Input_SettingsFile<-Input_SettingsFile

  #Perfrom filtering as selected
  if (Feature_Filtering ==  "Modified"){
    message("Here we apply the modified 80%-filtering rule that takes the class information (Column `Conditions`) into account, which additionally reduces the effect of missing values. REF: Yang et. al., (2015), doi: 10.3389/fmolb.2015.00004)")
    message(paste("filtering value selected:", Feature_Filt_Value))

    feat_filt_data <- as.data.frame(Input_data)
    feat_filt_Conditions <- Conditions[!Input_SettingsFile$Conditions=="CoRe_media"]

    if(CoRe== TRUE){ # remove CoRe_media samples for feature filtering
      feat_filt_data <- feat_filt_data %>% filter(!Input_SettingsFile$Conditions=="CoRe_media")
    }

    if(is.null(unique(feat_filt_Conditions)) ==  TRUE){
      stop("Condition information is missing from the Experimental design.")
    }
    if(length(unique(feat_filt_Conditions)) ==  1){
      stop("To perform the Modified feature filtering there have to be at least 2 different Conditions in the `Condition` column in the Experimental design. Consider using the Standard feature filtering option.")
    }

    miss <- c()
    # for (i in unique_conditions){
    split_Input <- split(feat_filt_data, feat_filt_Conditions) # splits data frame into a list of dataframes by condition

    for (m in split_Input){ # Select metabolites to be filtered for different conditions
      for(i in 1:ncol(m)) {
        if(length(which(is.na(m[,i]))) > (1-Feature_Filt_Value)*nrow(m))
          miss <- append(miss,i)
      }
    }
    #   }

    if(length(miss) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      filtered_matrix <- Input_data
      feat_file_res <- "There where no metabolites exluded"
      write.table(feat_file_res,row.names =  FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    } else {
      names<-unique(colnames(Input_data)[miss])
      message(length(unique(miss)) ," metabolites where removed: ", paste0(names, collapse = ", "))
      filtered_matrix <- Input_data[,-miss]
      write.table(unique(colnames(Input_data)[miss]),row.names = FALSE, file =  paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    }
  }
  if (Feature_Filtering ==  "Standard"){
    message("Here we apply the so-called 80%-filtering rule, which removes metabolites with missing values in more than 80% of samples. REF: Smilde et. al. (2005), Anal. Chem. 77, 6729–6736., doi:10.1021/ac051080y")
    message(paste("filtering value selected:", Feature_Filt_Value))

    feat_filt_data <- as.data.frame(Input_data)

    if(CoRe== TRUE){ # remove CoRe_media samples for feature filtering
      feat_filt_data <- feat_filt_data %>% filter(!Input_SettingsFile$Conditions=="CoRe_media")
    }

    split_Input <- feat_filt_data

    miss <- c()
    message("***Performing standard feature filtering***")
    for(i in 1:ncol(split_Input)) { # Select metabolites to be filtered for one condition
      if(length(which(is.na(split_Input[,i]))) > (1-Feature_Filt_Value)*nrow(split_Input))
        miss <- append(miss,i)
    }

    if(length(miss) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      filtered_matrix <- Input_data
      feat_file_res <- "There where no metabolites exluded"
      write.table(feat_file_res,row.names =  FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    } else {
      names<-unique(colnames(Input_data)[miss])
      message(length(unique(miss)) ," metabolites where removed: ", paste0(names, collapse = ", "))
      filtered_matrix <- Input_data[,-miss]
      write.table(unique(colnames(Input_data)[miss]),row.names =  FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    }
  }
  if (Feature_Filtering ==  "None"){
    warning("No feature filtering is selected.")
    filtered_matrix <- as.data.frame(Input_data)
  }

  filtered_matrix <- as.data.frame(mutate_all(as.data.frame(filtered_matrix), function(x) as.numeric(as.character(x))))



  ################################################
  ### ### ### Missing value Imputation ### ### ###

  # Do MVI for control media samples
  if(CoRe==TRUE){
    replaceNAdf <- filtered_matrix%>% filter( Input_SettingsFile$Conditions=="CoRe_media")

    # find metabolites with NA
    na_percentage <- colMeans(is.na(replaceNAdf)) * 100
    highNA_metabs <- na_percentage[na_percentage>20]

    # report metabolites with NA
    if(sum(na_percentage)>0){
      message("NA values were found in Control_media samples for metabolites.")
      if(sum(na_percentage>20)>0){
        message("Metabolites with high NA load in Control_media samples are: ",paste(names(highNA_metabs), collapse = ", "), ".")
      }
    }
    # if all values are NA set to 0
    replaceNAdf[,which(sapply(replaceNAdf, function(x)all(is.na(x))))]=0
    # If there is at least 1 value use the half minimum per feature
    replaceNAdf <- apply(replaceNAdf, 2,  function(x) {x[is.na(x)] <-  min(x, na.rm = TRUE)/2
    return(x)
    }) %>% as.data.frame()

    # replace the samples in the original dataframe
    filtered_matrix[rownames(filtered_matrix) %in% rownames(replaceNAdf), ] <- replaceNAdf
  }

  # Do MVI for the samples
  if (MVI ==  TRUE){
    message("Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0")

    NA_removed_matrix <- filtered_matrix %>% as.data.frame()
    for (feature  in colnames(NA_removed_matrix)){
      feature_data <- merge(NA_removed_matrix[feature] , Input_SettingsFile %>% select(Conditions), by= 0)
      feature_data <-column_to_rownames(feature_data, "Row.names")

      imputed_feature_data <- feature_data %>%
        group_by(Conditions) %>%
        mutate(across(all_of(feature), ~replace(., is.na(.), min(., na.rm = TRUE)*(MVI_Percentage/100))))

      NA_removed_matrix[[feature]] <- imputed_feature_data[[feature]]
    }

  } else if(MVI == FALSE){
    message("Missing value imputation is not performed.")
    stored_NA_positions <- which(is.na(filtered_matrix), arr.ind = TRUE)
    NA_removed_matrix <- replace(filtered_matrix, is.na(filtered_matrix), 0)
  }

  #######################################################
  ### ### ### Total Ion Current Normalization ### ### ###

  if (TIC_Normalization ==  TRUE){
    message("Total Ion Count (TIC) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022")
    RowSums <- rowSums(NA_removed_matrix)
    Median_RowSums <- median(RowSums) #This will built the median
    Data_TIC_Pre <- apply(NA_removed_matrix, 2, function(i) i/RowSums) #This is dividing the ion intensity by the total ion count
    Data_TIC <- Data_TIC_Pre*Median_RowSums #Multiplies with the median metabolite intensity
    Data_TIC <- as.data.frame(Data_TIC)
  } else {  # TIC_Normalization == FALSE
    Data_TIC <- as.data.frame(NA_removed_matrix)
    warning("***Total Ion Count normalization is not performed***")
    message("Total Ion Count (TIC) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022")
  }

  # Normalization QC plots
  # pre normalization
  log_NA_removed_matrix <- log(NA_removed_matrix) %>% t() %>% as.data.frame() # log tranforms the data
  medians <- apply(log_NA_removed_matrix, 2, median) # get median
  RLA_data_raw <- log_NA_removed_matrix - medians   # Subtract the medians from each column
  RLA_data_long <- pivot_longer(RLA_data_raw, cols = everything(), names_to = "Group")
  names(RLA_data_long)<- c("Samples", "Intensity")

  # Create the ggplot boxplot
  RLA_data_raw <- ggplot(RLA_data_long, aes(x = Samples, y = Intensity)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, color = "red", linetype = "solid") +
    labs(title = "Before Normalization")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


  # after normalization
  log_Data_TIC <- log(Data_TIC) %>% t() %>% as.data.frame()
  medians <- apply(log_Data_TIC, 2, median)
  RLA_data_norm <- log_Data_TIC - medians   # Subtract the medians from each column
  RLA_data_long <- pivot_longer(RLA_data_norm, cols = everything(), names_to = "Group")
  names(RLA_data_long)<- c("Samples", "Intensity")

  # Create the ggplot boxplot
  RLA_data_norm <- ggplot(RLA_data_long, aes(x = Samples, y = Intensity)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, color = "red", linetype = "solid") +
    labs(title = "Aftre Normalization")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


  norm_plots <- RLA_data_raw + RLA_data_norm

  if (ExportQCPlots == TRUE){
    ggsave(filename = paste0(Results_folder_Preprocessing_folder, "/Normalization.",Save_as_Plot), plot = norm_plots, width = 10,  height = 8)
  }


  # Start QC plot list
  qc_plot_list <- list()
  qc_plot_list_counter = 1

  if (CoRe ==  TRUE){
    CoRe_medias <-  Data_TIC[grep("CoRe_media", Conditions),]
    if(dim(CoRe_medias)[1]==1){
      warning("Only 1 CoRe_media sample was found. Thus, the consistency of the CoRe_media samples cannot be checked. It is assumed that the CoRe_media samples are already summed.")
      CoRe_media_df <- CoRe_medias %>% t() %>% as.data.frame()
      colnames(CoRe_medias) <- "CoRe_mediaMeans"

    }else{

      ## Media_control PCA
      media_pca_data <- merge(Input_SettingsFile %>% select(Conditions), Data_TIC, by=0) %>%
        column_to_rownames("Row.names") %>%
        mutate(Sample_type = case_when(Conditions == "CoRe_media" ~ "CoRe_media",
                                       TRUE ~ "Sample"))

      pca_QC_media <-invisible(MetaProViz::VizPCA(Input_data=media_pca_data %>%select(-Conditions, -Sample_type), Plot_SettingsInfo= c(color="Sample_type"),
                                                 Plot_SettingsFile= media_pca_data, OutputPlotName = "QC Media_samples",
                                                 Save_as_Plot =  NULL))

      qc_plot_list[[qc_plot_list_counter]] <- pca_QC_media
      qc_plot_list_counter = qc_plot_list_counter+1

      if (ExportQCPlots == TRUE){
        ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/PCA_Media_samples.",Save_as_Plot), plot = pca_QC_media, width = 10,  height = 8)
      }


      ## Check metabolite variance
      # Thresholds
      Therhold_cv = 1
      data_cv <- CoRe_medias

      ## Coefficient of Variation
      result_df <- apply(data_cv, 2,   function(x) { sd(x, na.rm =T)/  mean(x, na.rm =T) } ) %>% t()%>% as.data.frame()
      result_df[1, is.na(result_df[1,])]<- 0
      rownames(result_df)[1] <- "CV"

      result_df <- result_df %>% t()%>%as.data.frame() %>% rowwise() %>%
        mutate(High_var = CV > Therhold_cv) %>% as.data.frame()
      rownames(result_df)<- colnames(data_cv)

      high_var_metabs <- sum(result_df$High_var == TRUE)
      if(high_var_metabs>0){
        message(paste0(high_var_metabs, " of variables have high variability in the CoRe_media samples. Consider checking the pooled samples to decide whether to remove these metabolites or not."))
      }

      # Export/Save CV table ?
      # Export some QC plots?

      # Do sample outlier testing
      if(dim(CoRe_medias)[1]>3){


        Outlier_data <- CoRe_medias
        Outlier_data <- Outlier_data %>% mutate_all(.funs = ~ FALSE)

        while(high_var_metabs>0){

          #remove the furthest value from the mean
          if(high_var_metabs>1){
            max_var_pos <-  data_cv[,result_df$High_var == TRUE]  %>%
              as.data.frame() %>%
              mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
              summarise_all(.funs = ~ which.max(abs(.)))
          }else{
            max_var_pos <-  data_cv[,result_df$High_var == TRUE]  %>%
              as.data.frame() %>%
              mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
              summarise_all(.funs = ~ which.max(abs(.)))
            colnames(max_var_pos)<- colnames( data_cv)[result_df$High_var == TRUE]

          }

          # Remove rows based on positions
          for(i in 1:length(max_var_pos)){
            data_cv[max_var_pos[[i]],names(max_var_pos)[i]] <- NA
            Outlier_data[max_var_pos[[i]],names(max_var_pos)[i]] <- TRUE
          }

          # ReCalculate coefficient of variation for each column in the filtered data
          result_df <- apply(data_cv, 2,   function(x) { sd(x, na.rm =T)/  mean(x, na.rm =T) } ) %>% t()%>% as.data.frame()
          result_df[1, is.na(result_df[1,])]<- 0
          rownames(result_df)[1] <- "CV"

          result_df <- result_df %>% t()%>%as.data.frame() %>% rowwise() %>%
            mutate(High_var = CV > Therhold_cv) %>% as.data.frame()
          rownames(result_df)<- colnames(data_cv)

          high_var_metabs <- sum(result_df$High_var == TRUE)
        }


        data_cont <- Outlier_data %>% t() %>% as.data.frame()

        # List to store results
        fisher_test_results <- list()
        large_contingency_table <- matrix(0, nrow = 2, ncol = ncol(data_cont))

        for (i in 1:length(colnames(data_cont))) {
          sample = colnames(data_cont)[i]
          current_sample <- data_cont[, sample]

          contingency_table <- matrix(0, nrow = 2, ncol = 2)
          contingency_table[1, 1] <- sum(current_sample)
          contingency_table[2, 1] <- sum(!current_sample)
          contingency_table[1, 2] <- sum(rowSums(data_cont) - current_sample)
          contingency_table[2, 2] <- dim(data_cont %>% select(!all_of(sample)))[1]*dim(data_cont %>% select(!all_of(sample)))[2] -sum( rowSums(data_cont) - current_sample)

          # Fisher's exact test
          fisher_test_result <- fisher.test(contingency_table)
          fisher_test_results[[sample]] <- fisher_test_result

          # Calculate the sum of "TRUE" and "FALSE" for the current sample
          large_contingency_table[1, i] <- sum(current_sample)  # Sum of "TRUE"
          large_contingency_table[2, i] <- sum(!current_sample) # Sum of "FALSE"

        }

        # Convert the matrix into a data_contframe for better readability
        contingency_data_contframe <- as.data.frame(large_contingency_table)
        colnames(contingency_data_contframe) <- colnames(data_cont)
        rownames(contingency_data_contframe) <- c("High_var", "Low_var")

        contingency_data_contframe <- contingency_data_contframe %>% mutate(Total = rowSums(contingency_data_contframe))
        contingency_data_contframe <- rbind(contingency_data_contframe, Total= colSums(contingency_data_contframe))

        if (ExportQCPlots == TRUE){
          assign("CoRe_media_contingency_table", contingency_data_contframe, envir=.GlobalEnv)
        }


        different_samples <- c()
        for (sample in colnames(data_cont)) {
          p_value <- fisher_test_results[[sample]]$p.value
          if (p_value < 0.05) {  # Adjust the significance level as needed
            different_samples <- c(different_samples, sample)
          }
        }

        if(is.null(different_samples)==FALSE){
          warning("The CoRe_media samples ", paste(different_samples, collapse = ", "), " were found to be different from the rest. They will not be included in the sum of the CoRe_media samples.")
        }
        # Filter the CoRe_media samples
        CoRe_medias <- CoRe_medias %>% filter(!rownames(CoRe_medias) %in% different_samples)
      }
      CoRe_media_df <- as.data.frame(data.frame("CoRe_mediaMeans"=  colMeans( CoRe_medias, na.rm = TRUE)))

    }

    # Remove CoRe_media samples from the data
    Data_TIC <- Data_TIC[-grep("CoRe_media", Conditions),]

    Data_TIC_CoRe <- as.data.frame(t( apply(t(Data_TIC),2, function(i) i-CoRe_media_df$CoRe_mediaMeans)))  #Subtract from each sample the CoRe_media mean
    message("CoRe data are normalised using CoRe_norm_factor")
    Data_TIC <- apply(Data_TIC_CoRe, 2, function(i) i*CoRe_norm_factor)

    if (var(CoRe_norm_factor) ==  0){
      warning("The growth rate or growth factor for normalising the CoRe result, is the same for all samples")
    }
    # Remove CoRe_media samples from the data
    Input_SettingsFile <- Input_SettingsFile[Input_SettingsFile$Conditions!="CoRe_media",]
    Conditions <- Conditions[!Conditions=="CoRe_media"]
  }

  data_norm <- Data_TIC %>% as.data.frame()


  #####################################################
  ### ### ### Sample outlier identification ### ### ###

  message("Identification of outlier samples is performed using Hotellin's T2 test to define sample outliers in a mathematical way (Confidence = 0.99 ~ p.val < 0.01) REF: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3), 360–378, doi:https://doi.org/10.1214/aoms/1177732979.")
  message(paste("HotellinsConfidence value selected:", HotellinsConfidence))

  Outlier_filtering_loop = 10
  sample_outliers <- list()
  scree_plot_list <- list()
  outlier_plot_list <- list()
  metabolite_zero_var_total_list <- list()
  zero_var_metab_warning = FALSE
  k =  1
  a =  1
  for (loop in 1:Outlier_filtering_loop){   # here we do 10 rounds of hotelling filtering

    #################################################
    ### ### ### Zero variance metabolites ### ### ###
    metabolite_var <- as.data.frame( apply(data_norm, 2, var) %>% t()) # calculate each metabolites variance
    metabolite_zero_var_list <- list( colnames(metabolite_var)[which(metabolite_var[1,]==0)]) # takes the names of metabollites with zero variance and puts them in list

    if(sum(metabolite_var[1,]==0)==0){
      metabolite_zero_var_total_list[loop] <- 0
    } else if(sum(metabolite_var[1,]==0)>0){
      metabolite_zero_var_total_list[loop] <- metabolite_zero_var_list
      zero_var_metab_warning = TRUE # This is used later to print and save the zero variance metabolites if any are found.
    }

    for (metab in metabolite_zero_var_list){  # Remove the metabolites with zero variance from the data to do PCA
      data_norm <- data_norm %>% select(-all_of(metab))
    }

    ### ### PCA  ### ###
    PCA.res <- prcomp(data_norm, center =  TRUE, scale. =  TRUE)
    outlier_PCA_data <- data_norm
    outlier_PCA_data$Conditions <- Conditions

    pca_outlier <-invisible(MetaProViz::VizPCA(Input_data=data_norm, Plot_SettingsInfo= c(color="Conditions"),
                                               Plot_SettingsFile= outlier_PCA_data, OutputPlotName = paste("PCA outlier test filtering round ",loop),
                                               Save_as_Plot =  NULL))

    plot.new()
    plot(pca_outlier)
    outlier_plot_list[[k]] <- recordPlot()
    dev.off()
    k = k+1

    ### ### Scree plot ### ###
    inflect_df <- as.data.frame(c(1:length(PCA.res$sdev))) # get Scree plot values for inflection point calculation
    colnames(inflect_df) <- "x"
    inflect_df$y <- summary(PCA.res)$importance[2,]
    inflect_df$Cumulative <- summary(PCA.res)$importance[3,]
    screeplot_cumul <- format(round(inflect_df$Cumulative[1:20]*100, 1), nsmall = 1) #make cumulative variation labels for plot
    knee = inflection::uik(inflect_df$x,inflect_df$y) # Calculate the knee and select optimal number of components
    npcs = knee -1 #Note: we subtract 1 components from the knee cause the root of the knee is the PC that does not add something. npcs = 30

    # Make a scree plot with the selected component cut-off for HotellingT2 test
    screeplot <- factoextra::fviz_screeplot(PCA.res, main = paste("PCA Explained variance plot filtering round ",loop, sep = ""),
                                            addlabels = TRUE,
                                            ncp = 20,
                                            geom = c("bar", "line"),
                                            barfill = "grey",
                                            barcolor = "grey",
                                            linecolor = "black",linetype = 1) + theme_classic()+ geom_vline(xintercept = npcs+0.5, linetype = 2, color = "red") +
      annotate("text", x = c(1:20),y = -0.8,label = screeplot_cumul,col = "black", size = 3)

    plot.new()
    plot(screeplot)
    outlier_plot_list[[k]] <- recordPlot() # save plot
    dev.off()
    k = k+1

    ### ### HotellingT2 test for outliers ### ###
    data_hot <- as.matrix(PCA.res$x[,1:npcs])
    message("***Checking for outliers***")
    hotelling_qcc <- qcc::mqcc(data_hot, type = "T2.single",labels = rownames(data_hot),confidence.level = HotellinsConfidence, title = paste("Outlier filtering via HotellingT2 test filtering round ",loop,", with ",HotellinsConfidence, "% Confidence",  sep = ""), plot = FALSE)
    HotellingT2plot_data <- as.data.frame(hotelling_qcc$statistics)
    HotellingT2plot_data <- rownames_to_column(HotellingT2plot_data, "Samples")
    colnames(HotellingT2plot_data) <- c("Samples", "Group summary statisctics")
    outlier <- HotellingT2plot_data %>% filter(HotellingT2plot_data$`Group summary statisctics`>hotelling_qcc$limits[2])
    limits <- as.data.frame(hotelling_qcc$limits)
    legend <- colnames(HotellingT2plot_data[2])
    LegendTitle = "Limits"

    HotellingT2plot <- ggplot(HotellingT2plot_data, aes(x = Samples, y = `Group summary statisctics`, group = 1, fill = ))
    HotellingT2plot <- HotellingT2plot +
      geom_point(aes(x = Samples,y = `Group summary statisctics`), color = 'blue', size = 2) +
      geom_point(data = outlier, aes(x = Samples,y = `Group summary statisctics`), color = 'red',size = 3) +
      geom_line(linetype = 2)

    #draw the horizontal lines corresponding to the LCL,UCL
    HotellingT2plot <- HotellingT2plot + geom_hline(aes(yintercept = limits[,1]), color = "black", data = limits,  show.legend = F) +
      geom_hline(aes(yintercept = limits[,2], linetype = "UCL"), color = "red", data = limits, show.legend = T) +
      #only the LCl and UCL to be shown in y axis
      scale_y_continuous(breaks = sort(c(ggplot_build(HotellingT2plot)$layout$panel_ranges[[1]]$y.major_source, c(limits[,1],limits[,2]))))

    HotellingT2plot <- HotellingT2plot + theme_classic()
    HotellingT2plot <- HotellingT2plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    HotellingT2plot <- HotellingT2plot + ggtitle(paste("Outlier filtering via Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"))
    #plot1 <- plot1 + scale_linetype_manual(name = LegendTitle,values = "dashed") # this line instead of the next for bashed red line
    HotellingT2plot <- HotellingT2plot + scale_linetype_discrete(name = LegendTitle,)
    HotellingT2plot <- HotellingT2plot + theme(plot.title = element_text(size = 10, face = "bold")) +
      theme(axis.text = element_text(size = 12))

    plot.new()
    plot(HotellingT2plot)
    outlier_plot_list[[k]] <- recordPlot()
    dev.off()
    k = k+1

    ### Save the outlier detection plots in the outlier detection folder
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "/PCA_OD_round_" ,a ,".", Save_as_Plot, sep = ""),
           plot = pca_outlier, width = 10,height = 8)
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "//Scree_plot_OD_round_" ,a ,".",Save_as_Plot, sep = ""),
           plot = screeplot, width = 10,height = 8)
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "/Hotelling_OD_round_" ,a ,".", Save_as_Plot, sep = ""),
           plot = HotellingT2plot, width = 10,height = 8)
    a = a+1

    #Return the plots to environment:
    #The `outlier_plot_list` contains all the different plots that are saved as part of the Hotellins T2 test rounds. For each round three plots are recorded.
    assign("Outlier_Plots",  outlier_plot_list, envir=.GlobalEnv)

    # Here for the outliers we use confidence of 0.999 and p.val < 0.01.
    if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 0){ # loop for outliers until no outlier is detected
      data_norm <- data_norm  # filter the selected outliers from the data
      break
    } else if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 1){
      data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
      Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

      # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
      hotelling_qcc[["violations"]][["beyond.limits"]][1] <-   rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][1]]
      sample_outliers[loop] <- list(hotelling_qcc[["violations"]][["beyond.limits"]])

    } else {
      data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
      Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

      # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
      sm_out <- c() # list of outliers samples
      for (i in 1:length(hotelling_qcc[["violations"]][["beyond.limits"]])){
        sm_out <-  append(sm_out, rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][i]]  )
      }
      sample_outliers[loop] <- list(sm_out )
    }
  }

  pdf(file = paste(Results_folder_Preprocessing_Outlier_detection_folder, "/Outlier_testing.pdf", sep = ""), onefile = TRUE ) # or Outlier detection related plots
  for (plot in outlier_plot_list) {
    replayPlot(plot)
  }
  dev.off()


  # Print Outlier detection results about samples and metabolites
  if(length(sample_outliers) > 0){   # Print outlier samples
    message("There are possible outlier samples in the data") #This was a warning
    for (i in 1:length(sample_outliers)  ){
      message("Filtering round ",i ," Outlier Samples: ", paste( head(sample_outliers[[i]]) ," "))
    }
  }else{message("No sample outliers were found")}


  #######################################################
  ### ### ### Zero variance metabolites part2.### ### ###
  zero_var_metab_export_df <- data.frame(1,2)
  names(zero_var_metab_export_df) <- c("Filtering round","Metabolite")

  # Print zero variance metabolites
  if (zero_var_metab_warning==TRUE){
    warning("Metabolites with zero variance have been identified in the data. As scaling in PCA cannot be applied when features have zero variace, these metabolites are not taken into account for the outlier detection and the PCA plots.")
  }
  count = 1
  for (i in 1:length(metabolite_zero_var_total_list)){
    if (metabolite_zero_var_total_list[[i]] != 0){
      message("Filtering round ",i ,". Zero variance metabolites identified: ", paste( metabolite_zero_var_total_list[[i]] ," "))
      zero_var_metab_export_df[count,"Filtering round"] <- paste(i)
      zero_var_metab_export_df[count,"Metabolite"] <- paste(metabolite_zero_var_total_list[[i]])
      count = count +1
    }
  }
  if (zero_var_metab_warning==TRUE){
    write.table(zero_var_metab_export_df, row.names = FALSE, file =  paste(Results_folder_Preprocessing_folder,"/Zero_variance_metabolites",".csv",sep =  "")) # save zero var metabolite list
  }

  #############################################
  ### ### ### Make Output Dataframe ### ### ###

  total_outliers <- hash::hash() # make a dictionary
  if(length(sample_outliers) > 0){ # Create columns with outliers to merge to output dataframe
    for (i in 1:length(sample_outliers)  ){
      total_outliers[[paste("Outlier_filtering_round_",i, sep = "")]] <- sample_outliers[i]
    }
  }

  data_norm_filtered_full <- as.data.frame(Data_TIC)

  if(MVI == FALSE){ # if no MVI is selected put back the NAs in the dataframe
    data_norm_filtered_full[stored_NA_positions] <- NA
  }

  if(length(total_outliers) > 0){  # add outlier information to the full output dataframe
    data_norm_filtered_full$Outliers <- "no"
    for (i in 1:length(total_outliers)){
      for (k in 1:length( hash::values(total_outliers)[i] ) ){
        data_norm_filtered_full[as.character(hash::values(total_outliers)[[i]]) , "Outliers"] <- hash::keys(total_outliers)[i]
      }
    }
  }else{
    data_norm_filtered_full$Outliers <- "no"
  }

  data_norm_filtered_full <- data_norm_filtered_full %>% relocate(Outliers) #Put Outlier columns in the front
  data_norm_filtered_full <- merge(Input_SettingsFile, data_norm_filtered_full,  by = 0) # add the design in the output df (merge by rownames/sample names)
  rownames(data_norm_filtered_full) <- data_norm_filtered_full$Row.names
  data_norm_filtered_full$Row.names <- c()

  ################################################
  ### ### ### Quality Control (QC) PCA ### ### ###


  dtp <- data_norm_filtered_full %>%
    select(Conditions,Biological_Replicates, Outliers) %>%
    mutate(Outliers = case_when(Outliers == "no" ~ 'no',
                                Outliers == "Outlier_filtering_round_1" ~ ' Outlier_filtering_round = 1',
                                Outliers == "Outlier_filtering_round_2" ~ ' Outlier_filtering_round = 2',
                                Outliers == "Outlier_filtering_round_3" ~ ' Outlier_filtering_round = 3',
                                Outliers == "Outlier_filtering_round_4" ~ ' Outlier_filtering_round = 4',
                                TRUE ~ 'Outlier_filtering_round = or > 5'))
  dtp$Outliers <- relevel( as.factor(dtp$Outliers), ref="no")

  pca_QC <-invisible(MetaProViz::VizPCA(Input_data=as.data.frame(Data_TIC), Plot_SettingsInfo= c(color="Conditions", shape = "Outliers"),
                                        Plot_SettingsFile= dtp,OutputPlotName = "Quality Control PCA Condition clustering and Outlier check",
                                        Save_as_Plot =  NULL))


  qc_plot_list[[qc_plot_list_counter]] <- pca_QC
  qc_plot_list_counter = qc_plot_list_counter+1

  if (ExportQCPlots == TRUE){
    ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/PCA_Condition_Clustering.",Save_as_Plot), plot = pca_QC, width = 10,  height = 8)
  }


  if(is.null(Input_SettingsFile$Biological_Replicates)!= TRUE){
    pca_QC_repl <-invisible(MetaProViz::VizPCA(Input_data=as.data.frame(Data_TIC), Plot_SettingsInfo= c(color="Conditions", shape = "Biological_Replicates"),
                                               Plot_SettingsFile= dtp,OutputPlotName =  "Quality Control PCA replicate spread check",
                                               Save_as_Plot =  NULL))

    qc_plot_list[[qc_plot_list_counter]] <- pca_QC_repl

    if (ExportQCPlots == TRUE){
      ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/PCA_replicate_distribution.",Save_as_Plot), plot = pca_QC_repl, width = 10,  height = 8)
    }
  }

  #Return the QC plots to environment:
  assign("QC_Plots",  qc_plot_list, envir=.GlobalEnv)

  #########################################################
  ### ### ###  Make list with output dataframes ### ### ###

  output_list <- list()  #Here we make a list in which we will save the output
  preprocessing_output_list <- list(Input_SettingsFile = Original_Input_SettingsFile, Raw_data = as.data.frame(Input_data), Processed_data = data_norm_filtered_full)

  ##Write to file
  preprocessing_output_list_out <- lapply(preprocessing_output_list, function(x) rownames_to_column(x, "Sample_ID")) #  # use this line to make a sample_ID column in each dataframe
  writexl::write_xlsx(preprocessing_output_list_out, paste(Results_folder_Preprocessing_folder, "/Preprocessing_output.xlsx", sep = ""))#,showNA = TRUE)

  # Return the result
  assign("PreProcessing_res",  preprocessing_output_list, envir=.GlobalEnv)
}






############################################################
### ### ### Merge analytical replicates function ### ### ###
############################################################

#' Merges the analytical replicates of an experiment
#'
#' @param Input Dataframe which contains unique sample identifiers as row names the Experimental design and the metabolite numerical values in columns with metabolite identifiers as column names. Needs to have Conditions, Biological_Replicates and Analytical_Replicate columns
#'
#' @keywords Analyrical Replicate Merge
#' @export


ReplicateSum <- function(Input_data){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ###############################################
  ### ### ### Check Input Information ### ### ###

  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  }

  # Parse Condition and Replicate information
  if ( "Conditions" %in% colnames(Input_data)){
    Conditions <- Input_data$Conditions
  }else{
    stop("Column `Condition` is required.")
  }
  if ( "Biological_Replicates" %in% colnames(Input_data)){
    Biological_Replicates <- Input_data$Biological_Replicates
  }else{
    stop("Column `Biological_Replicates` is required.")
  }
  if ( "Analytical_Replicates" %in% colnames(Input_data)){
    Analytical_Replicates <- Input_data$Analytical_Replicates
  }else{
    stop("Column `Analytical_Replicates` is required.")
  }

  ############################################
  ### ### ### Create Output folder ### ### ###

  # This searches for a folder called "Results" within the current working directory and if its not found it creates one
  Results_folder = paste(getwd(), "/MetaProViz_Results_",Sys.Date(),  sep =  "")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  Results_folder_Preprocessing_folder = paste(Results_folder, "/Preprocessing", sep = "")
  if (!dir.exists(Results_folder_Preprocessing_folder)) {dir.create(Results_folder_Preprocessing_folder)}  # check and create folder

  ##############################################
  ### ### ### Load data and process ### ### ###

  #Load the data
  Input_data_numeric <-  select_if(Input_data, is.numeric) # take only the numeric values. This includes the replicate information
  Input_data_numeric <- merge(Input_data %>% select(Conditions), Input_data_numeric, by = 0)
  Input_data_numeric <- column_to_rownames(Input_data_numeric, "Row.names")

  # Make the replicate Sums
  Input_data_numeric_summed <- as.data.frame( Input_data_numeric %>%
                                                group_by(Biological_Replicates, Conditions) %>%
                                                summarise_all("mean") %>% select(-Analytical_Replicates))

  # Make a number of merged replicates column
  nReplicates <-  Input_data_numeric %>%
    group_by(Biological_Replicates, Conditions) %>%
    summarise_all("max") %>% ungroup() %>% select(Analytical_Replicates, Biological_Replicates, Conditions) %>% rename(n_Replicates_Summed = Analytical_Replicates)

  Input_data_numeric_summed <- merge(nReplicates,Input_data_numeric_summed, by = c("Conditions","Biological_Replicates"))%>%
    unite(UniqueID, c("Conditions","Biological_Replicates"), sep="_", remove=FALSE)%>% # Create a uniqueID
    column_to_rownames("UniqueID")# set UniqueID to rownames

  # Export result
  writexl::write_xlsx(Input_data_numeric_summed, paste(Results_folder_Preprocessing_folder, "/ReplicateSum_output.xlsx", sep = ""))#,showNA = TRUE)

  # Return the result

  return(Input_data_numeric_summed)
}




##########################################################################
### ### ### Metabolite detection estimation using pool samples ### ### ###
##########################################################################

#' Description
#'
#' @param Input_data DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. Can be either a full dataset or a dataset with just the pooled samples.
#' @param Input_SettingsFile  \emph{Optional: } DF which contains information about the samples when a full dataset is inserted as Input_data. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), has to exist.\strong{Default = NULL}
#' @param Input_SettingsInfo  \emph{Optional: } NULL or Named vector including the Pooled_Sample information (Name of the pooled samples in the Conditions in the Input_SettingsFile)  : c(Conditions="Pooled_Samples). \strong{Default = NULL}
#' @param Unstable_feature_remove  \emph{Optional: }  Parameter to automatically remove unstable(high variance) metabolites from the dataset. Used only when a full dataset is used as Input_data. \strong{Default = FALSE}
#' @param threshold_cv \emph{Optional: } Filtering threshold for high variance metabolites using the Coefficient of Variation. \strong{Default = 1}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt", ot NULL \strong{default: "csv"}
#'
#' @keywords Coefficient of Variation, high variance metabolites
#' @export


Pool_Estimation <- function(Input_data,
                            Input_SettingsFile = NULL,
                            Input_SettingsInfo = NULL,
                            Unstable_feature_remove = FALSE,
                            Therhold_cv = 1,
                            Save_as_Plot = "svg",
                            Save_as_Results = "csv" # txt or csv
){


  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Input_data <- Input_data
    }
  }
  # Check if the next lines work correctly in case of duplicated metabolties (colnames)
  if(sum( duplicated(colnames(Input_data))) > 0){
    doublons <- as.character(colnames(Input_data)[duplicated(colnames(Input_data))])#number of duplications
    data <-data[!duplicated(colnames(Input_data)),]#remove duplications
    warning("Input_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates.")
  }

  if(is.null(Input_SettingsFile)==FALSE){
    Test_match <- merge(Input_SettingsFile, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Input_SettingsFile"?
    if(nrow(Test_match) ==  0){
      stop("row.names Input_data need to match row.names Input_SettingsFile")
    } else{
      Input_data <- Input_data
    }

    # 2. Input_Settings
    if("Conditions" %in% names(Input_SettingsInfo)==TRUE){
      if(Input_SettingsInfo[["Conditions"]] %in% Input_SettingsFile[["Conditions"]]== FALSE ){
        stop("You have chosen Conditions = ",paste(Input_SettingsInfo[["Conditions"]]), ", ", paste(Input_SettingsInfo[["Conditions"]])," was not found in Plot_SettingsFile as sample Condition. Please insert the name of the pooled samples as stated in the Conditions column of the Input_SettingsFile."   )
      }
    }else{
      stop("The Conditions column was not found in the Input_SettingsFile. Either input a correct Input_SettingsFile or use a DF containing only yhe pooled samples as Input_data.")
    }
  }

  # 3. General parameters
  if(is_bare_logical(Unstable_feature_remove)==FALSE){
    stop("Check input. The Unstable_feature_remove value should be either =TRUE if metabolites with high variance are to be removed or =FALSE if not.")
  }
  if( is.numeric(threshold_cv)== FALSE | threshold_cv < 0){
    stop("Check input. The selected threshold_cv value should be a positive numeric value.")
  }
  if(is.null(Save_as_Plot)==FALSE){
    Save_as_Plot_options <- c("svg","pdf","png")
    if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
    }
  }
  if(is.null(Save_as_Results)==FALSE){
    Save_as_Results_options <- c("txt","csv", "xlsx" )
    if(Save_as_Results %in% Save_as_Results_options == FALSE){
      stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
    }
  }

  # Start QC plot list
  pool_plot_list <- list()
  pool_plot_list_counter = 1


  if(is.null(Input_SettingsFile)==TRUE){
    Input_numeric <- Input_data
  }else{
    pool_data <- Input_data[Input_SettingsFile[["Conditions"]]== Input_SettingsInfo[["Conditions"]],]
    Input_numeric <- pool_data
  }


  if(is.null(Save_as_Plot)==FALSE |is.null(Save_as_Results)==FALSE ){
    ## ------------ Create Results output folder ----------- ##
    name <- paste0("MetaProViz_Results_",Sys.Date())
    WorkD <- getwd()
    Results_folder <- file.path(WorkD, name)
    if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
    Results_folder_Preprocessing_folder = file.path(Results_folder, "Preprocessing")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
    if (!dir.exists(Results_folder_Preprocessing_folder)) {dir.create(Results_folder_Preprocessing_folder)}  # check and create folder
    Results_folder_Preprocessing_folder_Pool_Estimation = file.path(Results_folder_Preprocessing_folder, "Pool_Estimation")   # Create Outlier_Detection directory
    if (!dir.exists(Results_folder_Preprocessing_folder_Pool_Estimation)) {dir.create(Results_folder_Preprocessing_folder_Pool_Estimation)}

  }


  # Pool sample PCA
  if(is.null(Input_SettingsFile)==TRUE){
    pca_data <- Input_numeric
    pca_QC_pool <-invisible(MetaProViz::VizPCA(Input_data=pca_data, OutputPlotName = "QC Pool samples",Save_as_Plot =  NULL))
  }else{
    pca_data <- merge(Input_SettingsFile %>% select(Conditions), Input_data, by=0) %>%
      column_to_rownames("Row.names") %>%
      mutate(Sample_type = case_when(Conditions == Input_SettingsInfo[["Conditions"]] ~ "Pool",
                                     TRUE ~ "Sample"))

    pca_QC_pool <-invisible(MetaProViz::VizPCA(Input_data=pca_data %>%select(-Conditions, -Sample_type), Plot_SettingsInfo= c(color="Sample_type"),
                                               Plot_SettingsFile= pca_data, OutputPlotName = "QC Pool samples",
                                               Save_as_Plot =  NULL))
  }

  pool_plot_list[[pool_plot_list_counter]] <- pca_QC_pool
  pool_plot_list_counter = pool_plot_list_counter+1

  if (is.null(Save_as_Plot) == FALSE){
    ggsave(filename = paste0(Results_folder_Preprocessing_folder_Pool_Estimation, "/PCA_Pool_samples.",Save_as_Plot), plot = pca_QC_pool, width = 10,  height = 8)
  }

  ## Coefficient of Variation
  CV_data <- Input_numeric
  result_df <- apply(CV_data, 2,  function(x) { sd(x, na.rm =T)/  mean(x, na.rm =T) }  ) %>% t()%>% as.data.frame()
  rownames(result_df)[1] <- "CV"
  # Since we expect the pool samples to have very small variation since they "are the same sample" we go for CV=1 as threshold


  result_df_final <- result_df %>%
    t()%>%
    as.data.frame() %>%
    rowwise() %>%
    mutate(High_var = CV > threshold_cv) %>% as.data.frame()

  rownames(result_df_final)<- colnames(Input_numeric)
  result_df_final_out <- rownames_to_column(result_df_final,"Metabolite" )

  # Save results
  if(is.null(Save_as_Results)==FALSE){
    if (Save_as_Results == "xlsx"){
      xlsDMA <- file.path(Results_folder_Preprocessing_folder_Pool_Estimation,paste0("Pool_Estimation_CV", ".xlsx"))   # Save the DMA results table
      writexl::write_xlsx(result_df_final_out,xlsDMA, col_names = TRUE) # save the DMA result DF
    }else if (Save_as_Results == "csv"){
      csvDMA <- file.path(Results_folder_Preprocessing_folder_Pool_Estimation,paste0("Pool_Estimation_CV", ".csv"))
      write.csv(result_df_final_out,csvDMA,row.names = FALSE) # save the DMA result DF
    }else if (Save_as_Results == "txt"){
      txtDMA <- file.path(Results_folder_Preprocessing_folder_Pool_Estimation,paste0("Pool_Estimation_CV", ".txt"))
      write.table(result_df_final_out,txtDMA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
    }
  }

  #Make histogram of CVs
  HistCV <- invisible(ggplot(Pool_Estimation_result, aes(CV)) +
                        geom_histogram(aes(y=..density..), color="black", fill="white")+
                        geom_vline(aes(xintercept=Therhold_cv),
                                   color="darkred", linetype="dashed", size=1)+
                        geom_density(alpha=.2, fill="#FF6666") +
                        labs(title="Coefficient of Variation for metabolites of Pool samples",x="Coefficient of variation (CV)", y = "Frequency")+
                        theme_classic())

  pool_plot_list[[pool_plot_list_counter]] <- HistCV
  pool_plot_list_counter = pool_plot_list_counter+1

  if (is.null(Save_as_Plot) == FALSE){
    suppressMessages(ggsave(filename = paste0(Results_folder_Preprocessing_folder_Pool_Estimation, "/Pool_CV_Histogram",".",Save_as_Plot), plot = invisible(HistCV), width = 8,  height = 8))
  }

  ViolinCV <- invisible(ggplot(Pool_Estimation_result, aes(y=CV, x=High_var, label=row.names(Pool_Estimation_result)))+
                          geom_violin(alpha = 0.5 , fill="#FF6666")+
                          geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
                          #geom_point(position = position_jitter(seed = 1, width = 0.2))+
                          geom_text(aes(label=ifelse(CV>Therhold_cv,as.character(row.names(Pool_Estimation_result)),'')), hjust=0, vjust=0)+
                          labs(title="Coefficient of Variation for metabolites of Pool samples",x="Coefficient of variation (CV)", y = "Frequency")+
                          theme_classic())

  pool_plot_list[[pool_plot_list_counter]] <- ViolinCV
  pool_plot_list_counter = pool_plot_list_counter+1

  if (is.null(Save_as_Plot) == FALSE){
    suppressMessages(suppressWarnings(ggsave(filename = paste0(Results_folder_Preprocessing_folder_Pool_Estimation, "/Pool_CV_Violin",".",Save_as_Plot), plot = ViolinCV, width = 8,  height = 8)))
  }
  # Remove unstable metabolites
  filtered_Input_data <- Input_data # in case only pool samples in Input_data
  if(is.null(Input_SettingsFile)==FALSE){
    if(Unstable_feature_remove ==TRUE){
      unstable_metabs <- rownames(result_df_final)[result_df_final[["High_var_Metabs"]]]
      if(length(unstable_metabs)>0){
        filtered_Input_data <- Input_data %>% select(!unstable_metabs)
      }else{
        filtered_Input_data <- Input_data
      }
    }
  }


  DF_list <- list("Filtered_Input_data" = filtered_Input_data, "CV_result" = result_df_final_out )
  Pool_Estimation_res_list <- list("DFs"= DF_list,"Plots"=pool_plot_list )
  # Return the result
  assign("PoolEstimation_res",  Pool_Estimation_res_list, envir=.GlobalEnv)

}

