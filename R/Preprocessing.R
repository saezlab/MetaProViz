## ---------------------------
##
## Script name: PreProcessing
##
## Purpose of script: Metabolomics (raw ion counts) pre-processing, normalization, outlier detection and QC plots
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


#' Applies 80%-filtering rule, total-ion count normalization, missing value imputation and HotellingT2 outlier detection
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  NULL or Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "BiologicalReplicates" including numerical values. For CoRe = TRUE a CoRe_norm_factor = "Columnname_Input_SettingsFile" and CoRe_media = "Columnname_Input_SettingsFile", have to also be added. Column CoRe_norm_factor is used for normalization and CoRe_media is used to specify the name of the media controls in the Conditions.
#' @param FeatureFilt \emph{Optional: }If NULL, no feature filtering is performed. If set to "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = Modified}
#' @param FeatureFilt_Value \emph{Optional: } Percentage of feature filtering. \strong{Default = 0.8}
#' @param TIC \emph{Optional: } If TRUE, Total Ion Count normalization is performed. \strong{Default = TRUE}
#' @param MVI \emph{Optional: } If TRUE, Missing Value Imputation (MVI) based on half minimum is performed \strong{Default = TRUE}
#' @param MVI_Percentage \emph(Optional: ) Percentage (0 to 100)of imputed value based on the minimum value. \strong{Default = 50}
#' @param HotellinsConfidence \emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. Must be numeric.\strong{Default = 0.99}
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed and the CoRe value will be calculated. Please consider providing a Normalisation factor column called "CoRe_norm_factor" in your "Input_SettingsFile" DF, where the column "Conditions" matches. The normalisation factor must be a numerical value obtained from growth rate that has been obtained from a growth curve or growth factor that was obtained by the ratio of cell count/protein quantification at the start point to cell count/protein quantification at the end point.. Additionally control media samples have to be available in the "Input" DF and defined as "CoRe_media" samples in the "Conditions" column in the "Input_SettingsFile" DF. \strong{Default = FALSE}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. If set to NULL, plots are not saved. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } Select the file type of output table. Options are "csv", "xlsx", "txt". If set to NULL, plots are not saved. \strong{Default = "csv"}
#' @param PrintPlot  \emph{Optional: } If TRUE prints an overview of resulting plots. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong(Default = NULL)
#'
#' @keywords 80% filtering rule, Missing Value Imputation, Total Ion Count normalization, PCA, HotellingT2, multivariate quality control charts,
#' @export


###################################################
### ### ### Metabolomics pre-processing ### ### ###
###################################################

PreProcessing <- function(InputData,
                          SettingsFile_Sample,
                          SettingsInfo,
                          FeatureFilt = "Modified",
                          FeatureFilt_Value = 0.8,
                          TIC = TRUE,
                          MVI= TRUE,
                          MVI_Percentage=50,
                          HotellinsConfidence = 0.99,
                          CoRe = FALSE,
                          SaveAs_Plot = "svg",
                          SaveAs_Table = "csv",
                          PrintPlot = TRUE,
                          FolderPath = NULL
){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", # general scripting
                        "factoextra", # visualize PCA
                        "qcc", # for hotelling plots
                        "ggplot2", # For visualization PCA
                        "hash", # Dictionary in R for making column of outliers
                        "reshape", # for melting df for anova
                        "gridExtra",
                        "inflection",
                        "patchwork" # for output plot grid
  )
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          SettingsFile_Sample=SettingsFile_Sample,
                          SettingsFile_Metab=NULL,
                          SettingsInfo= SettingsInfo,
                          SaveAs_Plot=SaveAs_Plot,
                          SaveAs_Table=SaveAs_Table,
                          CoRe=CoRe,
                          PrintPlot= PrintPlot)


  # HelperFunction `CheckInput` Specific
  MetaProViz:::CheckInput_PreProcessing(InputData=InputData,
                                        SettingsFile_Sample=SettingsFile_Sample,
                                        SettingsInfo=SettingsInfo,
                                        CoRe=CoRe,
                                        FeatureFilt=FeatureFilt,
                                        FeatureFilt_Value=FeatureFilt_Value,
                                        TIC=TIC,
                                        MVI=MVI,
                                        MVI_Percentage=MVI_Percentage,
                                        HotellinsConfidence=HotellinsConfidence)

  ## ------------------  Create output folders  and path ------------------- ##
  if(is.null(SaveAs_Plot)==FALSE |is.null(SaveAs_Table)==FALSE ){
    Folder <- MetaProViz:::SavePath(FolderName= "Processing",
                                    FolderPath=FolderPath)

    SubFolder_P <- file.path(Folder, "PreProcessing")
    if (!dir.exists(SubFolder_P)) {dir.create(SubFolder_P)}
  }


  ## ------------------ Prepare the data ------------------- ##
  #InputData files:
  InputData <-as.data.frame(InputData)%>%
      mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .))#Make sure all 0 are changed to NAs

  InputData <- as.data.frame(mutate_all(as.data.frame(InputData), function(x) as.numeric(as.character(x))))

  ###################################################################################################################################
  ## ------------------ 1. Feature filtering ------------------- ##
  if(is.null(FeatureFilt)==FALSE){
    InputData_Filtered <- MetaProViz:::FeatureFiltering(InputData=InputData,
                                                    FeatureFilt=FeatureFilt,
                                                    FeatureFilt_Value=FeatureFilt_Value,
                                                    SettingsFile_Sample=SettingsFile_Sample,
                                                    SettingsInfo=SettingsInfo,
                                                    CoRe=CoRe)

    InputData_Filt <- InputData_Filtered[["DF"]]
  }else{
    InputData_Filt <- InputData
  }

  ## ------------------ 2. Missing value Imputation ------------------- ##
  if(MVI==TRUE){
    MVIRes<- MetaProViz:::MVImputation(InputData=InputData_Filt,
                                       SettingsFile_Sample=SettingsFile_Sample,
                                       SettingsInfo=SettingsInfo,
                                       CoRe=CoRe,
                                       MVI_Percentage=MVI_Percentage)
  }else{
    MVIRes<- InputData_Filt
  }

  ## ------------------  3. Total Ion Current Normalization ------------------- ##
  if(TIC==TRUE){
    #Perform TIC
    TICRes_List <- MetaProViz:::TICNorm(InputData=MVIRes,
                                        SettingsFile_Sample=SettingsFile_Sample,
                                        TIC=TIC)
    TICRes <- TICRes_List[["DF"]][["Data_TIC"]]

    #Add plots to PlotList
    PlotList <- list()
    PlotList[["RLAPlot"]] <- TICRes_List[["Plot"]][["RLA_BeforeTICNorm"]]
    PlotList[["RLAPlot_TICnorm"]] <- TICRes_List[["Plot"]][["RLA_AfterTICNorm"]]
    PlotList[["RLAPlot_BeforeAfter_TICnorm"]] <- TICRes_List[["Plot"]][["norm_plots"]]
  }else{
    TICRes <- MVIRes

    #Add plots to PlotList
    RLAPlot_List <- MetaProViz:::TICNorm(InputData=MVIRes,
                                         SettingsFile_Sample=SettingsFile_Sample,
                                         TIC=TIC)
    PlotList[["RLAPlot"]] <- RLAPlot_List[["Plot"]][["RLA_BeforeTICNorm"]]
  }

  ## ------------------ 4. CoRe media QC (blank) and normalization ------------------- ##
  if(CoRe ==TRUE){
   data_CoReNorm <- MetaProViz:::CoReNorm(InputData= TICRes,
                                               SettingsFile_Sample=SettingsFile_Sample,
                                               SettingsInfo=SettingsInfo)

    TICRes <- data_CoReNorm[["DF"]][["Core_Norm"]]
  }

  # ------------------ Final Output:
  data_norm <- TICRes %>% as.data.frame()

  ###################################################################################################################################
  ## ------------------ Sample outlier identification ------------------- ##
  OutlierRes <-  MetaProViz:::OutlierDetection(InputData= data_norm,
                                               SettingsFile_Sample=SettingsFile_Sample,
                                               SettingsInfo=SettingsInfo,
                                               CoRe=CoRe,
                                               HotellinsConfidence=HotellinsConfidence)

  ###################################################################################################################################
  ## ------------------ Return ------------------- ##
  ## ---- DFs
  if(is.null(FeatureFilt)==FALSE){#Add metabolites that where removed as part of the feature filtering
    if(length(InputData_Filtered[["RemovedMetabolites"]])==0){
      DFList <- list("InputData_RawData"= merge(as.data.frame(SettingsFile_Sample), as.data.frame(InputData), by="row.names")%>% column_to_rownames("Row.names"),
                     "Filtered_metabolites"= as.data.frame(list(FeatureFiltering = c(FeatureFilt),
                                                            FeatureFilt_Value = c(FeatureFilt_Value),
                                                            RemovedMetabolites = c("None"))),
                     "Preprocessing_output"=OutlierRes[["DF"]][["data_outliers"]])
    }else{
      DFList <- list("InputData_RawData"= merge(as.data.frame(SettingsFile_Sample), as.data.frame(InputData), by="row.names")%>% column_to_rownames("Row.names"),
                     "Filtered_metabolites"= as.data.frame(list(FeatureFiltering = rep(FeatureFilt, length(InputData_Filtered[["RemovedMetabolites"]])),
                                                           FeatureFilt_Value = rep(FeatureFilt_Value, length(InputData_Filtered[["RemovedMetabolites"]])),
                                                           RemovedMetabolites = InputData_Filtered[["RemovedMetabolites"]])),
                     "Preprocessing_output"=OutlierRes[["DF"]][["data_outliers"]])
    }
  }else{
    DFList <- list("InputData_RawData"= merge(as.data.frame(SettingsFile_Sample), as.data.frame(InputData), by="row.names")%>% column_to_rownames("Row.names"), "Preprocessing_output"=OutlierRes[["DF"]][["data_outliers"]])
  }

  if(CoRe ==TRUE){
     DFList_CoRe <- list( "CV_CoRe_blank"= data_CoReNorm[["DF"]][["CV_CoRe_blank"]],"Variation_ContigencyTable_CoRe_blank"=data_CoReNorm[["DF"]][["Contigency_table_CoRe_blank"]])
     DFList <- c(DFList, DFList_CoRe)
  }

  ## ---- Plots
  if(is.null(TIC)==FALSE){
    PlotList <- c(TICRes_List[["Plot"]], OutlierRes[["Plot"]])
  }else{
    PlotList <- c(RLAPlot_List[["Plot"]])
  }

  if(CoRe ==TRUE){
    PlotList <- c(PlotList , data_CoReNorm[["Plot"]])
  }

  Res_List <- list("DF"= DFList ,"Plot" =PlotList)

  # Save Plots and DFs
  #As row names are not saved we need to make row.names to column for the DFs that needs this:
  DFList[["InputData_RawData"]] <- DFList[["InputData_RawData"]]%>%rownames_to_column("Code")
  DFList[["Preprocessing_output"]] <- DFList[["Preprocessing_output"]]%>%rownames_to_column("Code")

  suppressMessages(suppressWarnings(
    MetaProViz:::SaveRes(InputList_DF=DFList,
                         InputList_Plot= PlotList,
                         SaveAs_Table=SaveAs_Table,
                         SaveAs_Plot=SaveAs_Plot,
                         FolderPath= SubFolder_P,
                         FileName= "PreProcessing",
                         CoRe=CoRe,
                         PrintPlot=PrintPlot)))


  #Return
  invisible(return(Res_List))
}




############################################################
### ### ### Merge analytical replicates function ### ### ###
############################################################

#' Merges the analytical replicates of an experiment
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#'#@param SettingsFile_Sample DF which contains information about the samples Column "Conditions", "Biological_replicates" and "Analytical_Replicates has to exist.
#' @param SettingsInfo  \emph{Optional: } Named vector including the Conditions and Replicates information: c(Conditions="ColumnNameConditions", Biological_Replicates="ColumnName_SettingsFile_Sample", Analytical_Replicates="ColumnName_SettingsFile_Sample").\strong{Default = NULL}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt", ot NULL \strong{default: "csv"}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong(Default = NULL)
#'
#' @keywords Analytical Replicate Merge
#' @export


ReplicateSum <- function(InputData,
                         SettingsFile_Sample,
                         SettingsInfo = c(Conditions="Conditions", Biological_Replicates="Biological_Replicates", Analytical_Replicates="Analytical_Replicates"),
                         SaveAs_Table = "csv",
                         FolderPath = NULL){
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          SettingsFile_Sample=SettingsFile_Sample,
                          SettingsFile_Metab=NULL,
                          SettingsInfo = SettingsInfo,
                          SaveAs_Plot=NULL,
                          SaveAs_Table=SaveAs_Table,
                          CoRe=FALSE,
                          PrintPlot=FALSE)

  # `CheckInput` Specific
  if(SettingsInfo[["Conditions"]] %in% colnames(SettingsFile_Sample)){
   # Conditions <- InputData[[SettingsInfo[["Conditions"]] ]]
  }else{
    stop("Column `Conditions` is required.")
  }
  if(SettingsInfo[["Biological_Replicates"]] %in% colnames(SettingsFile_Sample)){
    #Biological_Replicates <- InputData[[SettingsInfo[["Biological_Replicates"]]]]
  }else{
    stop("Column `Biological_Replicates` is required.")
  }
  if(SettingsInfo[["Analytical_Replicates"]] %in% colnames(SettingsFile_Sample)){
    #Analytical_Replicates <- InputData[[SettingsInfo[["Analytical_Replicates"]]]]
  }else{
    stop("Column `Analytical_Replicates` is required.")
  }

  ## ------------ Create Results output folder ----------- ##
  if(is.null(SaveAs_Table)==FALSE ){
    Folder <- MetaProViz:::SavePath(FolderName= "Processing",
                                    FolderPath=FolderPath)
    SubFolder <- file.path(Folder, "ReplicateSum")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }

  ## ------------  Load data and process  ----------- ##
  Input <- merge(x= SettingsFile_Sample%>% select(!!SettingsInfo[["Conditions"]], !!SettingsInfo[["Biological_Replicates"]], !!SettingsInfo[["Analytical_Replicates"]]),
                 y= InputData,
                 by="row.names")%>%
    column_to_rownames("Row.names")%>%
    dplyr::rename("Conditions"=SettingsInfo[["Conditions"]],
                  "Biological_Replicates"=SettingsInfo[["Biological_Replicates"]],
                  "Analytical_Replicates"=SettingsInfo[["Analytical_Replicates"]])

  # Make the replicate Sums
  Input_data_numeric_summed <- as.data.frame(Input %>%
                                               group_by(Biological_Replicates, Conditions) %>%
                                               summarise_all("mean") %>% select(-Analytical_Replicates))

  # Make a number of merged replicates column
  nReplicates <-  Input %>%
    group_by(Biological_Replicates, Conditions) %>%
    summarise_all("max") %>% ungroup() %>% select(Analytical_Replicates, Biological_Replicates, Conditions) %>%
    dplyr::rename("n_AnalyticalReplicates_Summed "= "Analytical_Replicates")

  Input_data_numeric_summed <- merge(nReplicates,Input_data_numeric_summed, by = c("Conditions","Biological_Replicates"))%>%
    unite(UniqueID, c("Conditions","Biological_Replicates"), sep="_", remove=FALSE)%>% # Create a uniqueID
    column_to_rownames("UniqueID")# set UniqueID to rownames

  #--------------- return ------------------##
  MetaProViz:::SaveRes(InputList_DF=list("Sum_AnalyticalReplicates"=Input_data_numeric_summed%>%rownames_to_column("Code")),
                       InputList_Plot = NULL,
                       SaveAs_Table=SaveAs_Table,
                       SaveAs_Plot=NULL,
                       FolderPath= SubFolder,
                       FileName= "Sum_AnalyticalReplicates",
                       CoRe=FALSE,
                       PrintPlot=FALSE)

  #Return
  invisible(return(Input_data_numeric_summed))
}




##########################################################################
### ### ### Metabolite detection estimation using pool samples ### ### ###
##########################################################################

#' Description
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. Can be either a full dataset or a dataset with only the pool samples.
#' @param SettingsFile_Sample  \emph{Optional: } DF which contains information about the samples when a full dataset is inserted as Input_data. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), has to exist.\strong{Default = NULL}
#' @param SettingsInfo  \emph{Optional: } NULL or Named vector including the Conditions and PoolSample information (Name of the Conditions column and Name of the pooled samples in the Conditions in the Input_SettingsFile)  : c(Conditions="ColumnNameConditions, PoolSamples=NamePoolCondition. If no Conditions is added in the Input_SettingsInfo, it is assumed that the conditions column is named 'Conditions' in the Input_SettingsFile. ). \strong{Default = NULL}
#' @param CutoffCV \emph{Optional: } Filtering cutoff for high variance metabolites using the Coefficient of Variation. \strong{Default = 1}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt", ot NULL \strong{default: "csv"}
#' @param PrintPlot \emph{Optional: } If TRUE prints an overview of resulting plots. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong(Default = NULL)
#'
#' @keywords Coefficient of Variation, high variance metabolites
#' @export


PoolEstimation <- function(InputData,
                           SettingsFile_Sample = NULL,
                           SettingsInfo = NULL,
                           CutoffCV = 100,
                           SaveAs_Plot = "svg",
                           SaveAs_Table = "csv", # txt or csv
                           PrintPlot=TRUE,
                           FolderPath = NULL){

  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse")

  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`
  MetaProViz:::CheckInput(InputData=InputData,
                          SettingsFile_Sample=SettingsFile_Sample,
                          SettingsFile_Metab=NULL,
                          SettingsInfo=SettingsInfo,
                          SaveAs_Plot=SaveAs_Plot,
                          SaveAs_Table=SaveAs_Table,
                          CoRe=FALSE,
                          PrintPlot = PrintPlot)

  # `CheckInput` Specific
  if(is.null(SettingsFile_Sample)==FALSE){
    if("Conditions" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["Conditions"]] %in% colnames(SettingsFile_Sample)== FALSE ){
        stop("You have chosen Conditions = ",paste(SettingsInfo[["Conditions"]]), ", ", paste(SettingsInfo[["Conditions"]])," was not found in SettingsFile_Sample as column. Please insert the name of the experimental conditions as stated in the SettingsFile_Sample."   )
      }
    }
    if("PoolSamples" %in% names(SettingsInfo)==TRUE){
      if(SettingsInfo[["PoolSamples"]] %in% SettingsFile_Sample[["Conditions"]] == FALSE ){
        stop("You have chosen PoolSamples = ",paste(SettingsInfo[["PoolSamples"]] ), ", ", paste(SettingsInfo[["PoolSamples"]] )," was not found in SettingsFile_Sample as sample condition. Please insert the name of the pool samples as stated in the Conditions column of the SettingsFile_Sample."   )
      }
    }
  }

  if(is.numeric(CutoffCV)== FALSE | CutoffCV < 0){
    stop("Check input. The selected CutoffCV value should be a positive numeric value.")
  }

  ## ------------------  Create output folders  and path ------------------- ##
  if(is.null(SaveAs_Plot)==FALSE |is.null(SaveAs_Table)==FALSE ){
    Folder <- MetaProViz:::SavePath(FolderName= "Processing",
                                    FolderPath=FolderPath)

    SubFolder <- file.path(Folder, "PoolEstimation")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }



  ## ------------------ Prepare the data ------------------- ##
  #InputData files:
  if(is.null(SettingsFile_Sample)==TRUE){
    PoolData <- InputData
    PoolData[PoolData == 0] <- NA
  }else{
    PoolData <- InputData[SettingsFile_Sample[["Conditions"]] == SettingsInfo[["PoolSamples"]],]
    PoolData[PoolData == 0] <- NA
  }


  ###################################################################################################################################
  ## ------------------ Coefficient of Variation ------------------- ##
  result_df <- apply(PoolData, 2,  function(x) { (sd(x, na.rm =T)/  mean(x, na.rm =T))*100 }  ) %>% t()%>% as.data.frame()
  rownames(result_df)[1] <- "CV"

  NAvector <- apply(PoolData, 2,  function(x) {(sum(is.na(x))/length(x))*100 })# Calculate the NAs

  # Create Output DF
  result_df_final <- result_df %>%
    t()%>% as.data.frame() %>% rowwise() %>%
    mutate(HighVar = CV > CutoffCV) %>% as.data.frame()

  result_df_final$MissingValuePercentage <- NAvector

  rownames(result_df_final)<- colnames(InputData)
  result_df_final_out <- rownames_to_column(result_df_final,"Metabolite" )

  # Remove Metabolites from InputData based on CutoffCV
  if(is.null(SettingsFile_Sample)==FALSE){
      unstable_metabs <- rownames(result_df_final)[result_df_final[["HighVar_Metabs"]]]
      if(length(unstable_metabs)>0){
        filtered_Input_data <- InputData %>% select(!unstable_metabs)
      }else{
        filtered_Input_data <- NULL
  }
      }else{
    filtered_Input_data <- NULL
  }

  ## ------------------ QC plots ------------------- ##
  # Start QC plot list
  PlotList <- list()

  # 1. Pool Sample PCA
  dev.new()
  if(is.null(SettingsFile_Sample)==TRUE){
    pca_data <- PoolData
    pca_QC_pool <-invisible(MetaProViz::VizPCA(InputData=pca_data,
                                               PlotName = "QC Pool samples",
                                               SaveAs_Plot =  NULL))
  }else{
    pca_data <- merge(SettingsFile_Sample %>% select(Conditions), InputData, by=0) %>%
      column_to_rownames("Row.names") %>%
      mutate(Sample_type = case_when(Conditions == SettingsInfo[["PoolSamples"]] ~ "Pool",
                                     TRUE ~ "Sample"))

    pca_QC_pool <-invisible(MetaProViz::VizPCA(InputData=pca_data %>%select(-Conditions, -Sample_type),
                                               SettingsInfo= c(color="Sample_type"),
                                               SettingsFile_Sample= pca_data,
                                               PlotName = "QC Pool samples",
                                               SaveAs_Plot =  NULL))
  }
  dev.off()
  PlotList [["PCAPlot_PoolSamples"]] <- pca_QC_pool[["Plot_Sized"]][["QC Pool samples"]]


  # 2. Histogram of CVs
  HistCV <-suppressWarnings(invisible(ggplot(result_df_final_out, aes(CV)) +
                        geom_histogram(aes(y=after_stat(density)), color="black", fill="white")+
                        geom_vline(aes(xintercept=CutoffCV),
                                   color="darkred", linetype="dashed", size=1)+
                        geom_density(alpha=.2, fill="#FF6666") +
                        labs(title="Coefficient of Variation for metabolites of Pool samples",x="Coefficient of variation (CV%)", y = "Frequency")+
                        theme_classic()))

  PlotList [["Histogram_CV-PoolSamples"]] <- HistCV

  # 2. ViolinPlot of CVs
  ViolinCV <- invisible(ggplot(result_df_final_out, aes(y=CV, x=HighVar, label=row.names(result_df_final_out)))+
                          geom_violin(alpha = 0.5 , fill="#FF6666")+
                          geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
                          #geom_point(position = position_jitter(seed = 1, width = 0.2))+
                          geom_text(aes(label=ifelse(CV>CutoffCV,as.character(row.names(result_df_final_out)),'')), hjust=0, vjust=0)+
                          labs(title="Coefficient of Variation for metabolites of Pool samples",x="Metabolites", y = "Coefficient of variation (CV%)")+
                          theme_classic())

  PlotList [["ViolinPlot_CV-PoolSamples"]] <- ViolinCV

  ###################################################################################################################################
  ## ------------------ Return and Save ------------------- ##
  #Save
  if(is.null(filtered_Input_data)==FALSE){
    DF_list <- list("InputData" = InputData, "Filtered_InputData" = filtered_Input_data, "CV" = result_df_final_out )
  }else{
    DF_list <- list("InputData" = InputData, "CV" = result_df_final_out)
  }
  ResList <- list("DF"= DF_list,"Plot"=PlotList)

  #Save
  DF_list[["InputData"]]<-  DF_list[["InputData"]]%>%rownames_to_column("Code")

  MetaProViz:::SaveRes(InputList_DF=DF_list,
                      InputList_Plot = PlotList,
                      SaveAs_Table=SaveAs_Table,
                      SaveAs_Plot=SaveAs_Plot,
                      FolderPath= SubFolder,
                      FileName= "PoolEstimation",
                      CoRe=FALSE,
                      PrintPlot=PrintPlot)

  #Return
  invisible(return(ResList))
}

################################################################################################
### ### ### PreProcessing helper function: Internal Function to check function input ### ### ###
################################################################################################

#' Check input parameters
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Sample Passed to main function MetaProViz::PreProcessing()
#' @param SettingsInfo Passed to main function MetaProViz::PreProcessing()
#' @param CoRe Passed to main function MetaProViz::PreProcessing()
#' @param FeatureFilt Passed to main function MetaProViz::PreProcessing()
#' @param FeatureFilt_Value Passed to main function MetaProViz::PreProcessing()
#' @param TIC Passed to main function MetaProViz::PreProcessing()
#' @param MVI Passed to main function MetaProViz::PreProcessing()
#' @param MVI_Percentage Passed to main function MetaProViz::PreProcessing()
#' @param HotellinsConfidence Passed to main function MetaProViz::PreProcessing()
#'
#' @keywords Input check
#' @noRd
#'
#'

CheckInput_PreProcessing <- function(InputData,
                                     SettingsFile_Sample,
                                     SettingsInfo,
                                     CoRe,
                                     FeatureFilt,
                                     FeatureFilt_Value,
                                     TIC,
                                     MVI,
                                     MVI_Percentage,
                                     HotellinsConfidence){
  if(is.vector(SettingsInfo)==TRUE){
    #-------------SettingsInfo
    #CoRe
    if(CoRe == TRUE){   # parse CoRe normalisation factor
      message("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
      if("CoRe_media" %in% names(SettingsInfo)){
        if(length(grep(SettingsInfo[["CoRe_media"]], SettingsFile_Sample[[SettingsInfo[["Conditions"]]]])) < 1){     # Check for CoRe_media samples
          stop("No CoRe_media samples were provided in the 'Conditions' in the SettingsFile_Sample. For a CoRe experiment control media samples without cells have to be measured and be added in the 'Conditions'
           column labeled as 'CoRe_media' (see @param section). Please make sure that you used the correct labelling or whether you need CoRe = FALSE for your analysis")
        }
      }

      if ("CoRe_norm_factor" %in% names(SettingsInfo)==FALSE){
        warning("No growth rate or growth factor provided for normalising the CoRe result, hence CoRe_norm_factor set to 1 for each sample")
      }
    }
  }

  #-------------General parameters
  Feature_Filtering_options <- c("Standard","Modified")
  if(FeatureFilt %in% Feature_Filtering_options == FALSE & is.null(FeatureFilt)==FALSE){
    stop("Check input. The selected FeatureFilt option is not valid. Please set to NULL or select one of the folowwing: ",paste(Feature_Filtering_options,collapse = ", "),"." )
  }
  if(is.numeric(FeatureFilt_Value) == FALSE |FeatureFilt_Value > 1 | FeatureFilt_Value < 0){
    stop("Check input. The selected FeatureFilt_Value should be numeric and between 0 and 1.")
  }
  if(is.logical(TIC) == FALSE){
    stop("Check input. The TIC should be either `TRUE` if TIC normalization is to be performed or `FALSE` if no data normalization is to be applied.")
  }
  if(is.logical(MVI) == FALSE){
    stop("Check input. MVI value should be either `TRUE` if mising value imputation should be performed or `FALSE` if not.")
  }
  if(is.numeric(MVI_Percentage)== FALSE |HotellinsConfidence > 100 | HotellinsConfidence < 0){
    stop("Check input. The selected MVI_Percentage value should be numeric and between 0 and 100.")
  }
  if( is.numeric(HotellinsConfidence)== FALSE |HotellinsConfidence > 1 | HotellinsConfidence < 0){
    stop("Check input. The selected Filtering value should be numeric and between 0 and 1.")
  }

}


################################################################################################
### ### ### PreProcessing helper function: FeatureFiltering ### ### ###
################################################################################################

#' FeatureFiltering
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Sample Passed to main function MetaProViz::PreProcessing()
#' @param SettingsInfo Passed to main function MetaProViz::PreProcessing()
#' @param CoRe Passed to main function MetaProViz::PreProcessing()
#' @param FeatureFilt Passed to main function MetaProViz::PreProcessing()
#' @param FeatureFilt_Value Passed to main function MetaProViz::PreProcessing()
#'
#' @keywords feature filtering
#' @noRd
#'

FeatureFiltering <-function(InputData, FeatureFilt, FeatureFilt_Value, SettingsFile_Sample, SettingsInfo, CoRe){
  ## ------------------ Prepare the data ------------------- ##
  feat_filt_data <- as.data.frame(replace(InputData, InputData==0, NA))

  if(CoRe== TRUE){ # remove CoRe_media samples for feature filtering
    feat_filt_data <- feat_filt_data %>% filter(!SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] ==SettingsInfo[["CoRe_media"]])
    Feature_Filtering <- paste0(FeatureFilt, "_CoRe")
  }

  ## ------------------ Perform filtering ------------------ ##
  if(FeatureFilt ==  "Modified"){
    message("Here we apply the modified 80%-filtering rule that takes the class information (Column `Conditions`) into account, which additionally reduces the effect of missing values. REF: Yang et. al., (2015), doi: 10.3389/fmolb.2015.00004)")
    message(paste("filtering value selected:", FeatureFilt_Value))

    if(CoRe== TRUE){
      feat_filt_Conditions <- SettingsFile_Sample[[SettingsInfo[["Conditions"]]]][!SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]]]
    }else{
      feat_filt_Conditions <- SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]
    }

    if(is.null(unique(feat_filt_Conditions)) ==  TRUE){
      stop("Conditions information is missing.")
    }
    if(length(unique(feat_filt_Conditions)) ==  1){
      stop("To perform the Modified feature filtering there have to be at least 2 different Conditions in the `Condition` column in the Experimental design. Consider using the Standard feature filtering option.")
    }

    miss <- c()
    split_Input <- split(feat_filt_data, feat_filt_Conditions) # split data frame into a list of dataframes by condition

    for (m in split_Input){ # Select metabolites to be filtered for different conditions
      for(i in 1:ncol(m)) {
        if(length(which(is.na(m[,i]))) > (1-FeatureFilt_Value)*nrow(m))
          miss <- append(miss,i)
      }
    }

    if(length(miss) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      filtered_matrix <- InputData
      feat_file_res <- "There where no metabolites exluded"
    }else{
      names<-unique(colnames(InputData)[miss])
      message(length(unique(miss)) ," metabolites where removed: ", paste0(names, collapse = ", "))
      filtered_matrix <- InputData[,-miss]
    }
  }else if(Feature_Filtering ==  "Standard"){
    message("Here we apply the so-called 80%-filtering rule, which removes metabolites with missing values in more than 80% of samples. REF: Smilde et. al. (2005), Anal. Chem. 77, 6729–6736., doi:10.1021/ac051080y")
    message(paste("filtering value selected:", FeatureFilt_Value))

    split_Input <- feat_filt_data

    miss <- c()
    for(i in 1:ncol(split_Input)) { # Select metabolites to be filtered for one condition
      if(length(which(is.na(split_Input[,i]))) > (1-FeatureFilt_Value)*nrow(split_Input))
        miss <- append(miss,i)
    }

    if(length(miss) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      filtered_matrix <- InputData
      feat_file_res <- "There where no metabolites exluded"
    }else{
      names<-unique(colnames(InputData)[miss])
      message(length(unique(miss)) ," metabolites where removed: ", paste0(names, collapse = ", "))
      filtered_matrix <- InputData[,-miss]
    }
  }

  ## ------------------ Return ------------------ ##
  features_filtered <- unique(colnames(InputData)[miss]) %>% as.vector()
  filtered_matrix <- as.data.frame(mutate_all(as.data.frame(filtered_matrix), function(x) as.numeric(as.character(x))))

  Filtered_results <- list("DF"= filtered_matrix , "RemovedMetabolites" = features_filtered)
  invisible(return(Filtered_results))
}



################################################################################################
### ### ### PreProcessing helper function: Missing Value imputation ### ### ###
################################################################################################

#' MVI
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Sample Passed to main function MetaProViz::PreProcessing()
#' @param SettingsInfo Passed to main function MetaProViz::PreProcessing()
#' @param CoRe Passed to main function MetaProViz::PreProcessing()
#' @param MVI_Percentage Passed to main function MetaProViz::PreProcessing()
#'
#' @keywords Half minimum missing value imputation
#' @noRd
#'

MVImputation <-function(InputData, SettingsFile_Sample, SettingsInfo, CoRe, MVI_Percentage){
  ## ------------------ Prepare the data ------------------- ##
  filtered_matrix <- InputData
  filtered_matrix[filtered_matrix == 0] <- NA

  ## ------------------ Perform MVI ------------------ ##
  if(CoRe==TRUE){
    replaceNAdf <- filtered_matrix%>% filter(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]])

    # find metabolites with NA
    na_percentage <- colMeans(is.na(replaceNAdf)) * 100
    highNA_metabs <- na_percentage[na_percentage>20 & na_percentage<100]
    OnlyNA_metabs <- na_percentage[na_percentage==100]

    # report metabolites with NA
    if(sum(na_percentage)>0){
      message("NA values were found in Control_media samples for metabolites. For metabolites including NAs MVI is performed unless all samples of a metabolite are NA.")
      if(sum(na_percentage>20 & na_percentage<100)>0){
        message("Metabolites with high NA load (>20%) in Control_media samples are: ",paste(names(highNA_metabs), collapse = ", "), ".")
      }
      if(sum(na_percentage==100)>0){
        message("Metabolites with only NAs (=100%) in Control_media samples are: ",paste(names(highNA_metabs), collapse = ", "), ". Those NAs are set zero as we consider them true zeros")
      }
    }

    # if all values are NA set to 0
    replaceNAdf_zero <- as.data.frame(lapply(replaceNAdf, function(x) if(all(is.na(x))) replace(x, is.na(x), 0) else x))
    colnames( replaceNAdf_zero) <-  colnames(replaceNAdf)

    # If there is at least 1 value use the half minimum per feature
    replaceNAdf_Zero_MVI <- apply( replaceNAdf_zero, 2,  function(x) {x[is.na(x)] <-  min(x, na.rm = TRUE)/2
    return(x)
    }) %>% as.data.frame()

    # replace the samples in the original dataframe
    filtered_matrix[rownames(filtered_matrix) %in% rownames(replaceNAdf_Zero_MVI), ] <- replaceNAdf_Zero_MVI
  }

  # Do MVI for the samples
  message("Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0")

  NA_removed_matrix <- filtered_matrix %>% as.data.frame()
  for (feature  in colnames(NA_removed_matrix)){
    feature_data <- merge(NA_removed_matrix[feature] , SettingsFile_Sample %>% select(Conditions), by= 0)
    feature_data <-column_to_rownames(feature_data, "Row.names")

    imputed_feature_data <- feature_data %>%
      group_by(Conditions) %>%
      mutate(across(all_of(feature), ~replace(., is.na(.), min(., na.rm = TRUE)*(MVI_Percentage/100))))

    NA_removed_matrix[[feature]] <- imputed_feature_data[[feature]]
  }

  ## ------------------ Return ------------------ ##
  invisible(return(NA_removed_matrix))
}


################################################################################################
### ### ### PreProcessing helper function: Total ion Count Normalization ### ### ###
################################################################################################

#' TIC
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Sample Passed to main function MetaProViz::PreProcessing()
#' @param TIC Passed to main function MetaProViz::PreProcessing()
#'
#' @keywords total ion count normalisation
#' @noRd
#'

TICNorm <-function(InputData, SettingsFile_Sample, TIC){
  ## ------------------ Prepare the data ------------------- ##
  NA_removed_matrix <- InputData
  NA_removed_matrix[is.na(NA_removed_matrix)] <- 0#replace NA with 0

  ## ------------------ QC plot ------------------- ##
  ### Before TIC Normalization
  log_NA_removed_matrix <- log(NA_removed_matrix) %>% t() %>% as.data.frame() # log tranforms the data
  medians <- apply(log_NA_removed_matrix, 2, median) # get median
  RLA_data_raw <- log_NA_removed_matrix - medians   # Subtract the medians from each column
  RLA_data_long <- pivot_longer(RLA_data_raw, cols = everything(), names_to = "Group")
  names(RLA_data_long)<- c("Samples", "Intensity")
  RLA_data_long <- as.data.frame(RLA_data_long)
  for (row in 1:nrow(RLA_data_long)){ # add conditions
    RLA_data_long[row,"Conditions"] <- SettingsFile_Sample[rownames(SettingsFile_Sample) %in%RLA_data_long[row,1],"Conditions"]
  }

  # Create the ggplot boxplot
  RLA_data_raw <- ggplot(RLA_data_long, aes(x = Samples, y = Intensity, color = Conditions)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, color = "red", linetype = "solid") +
    labs(title = "Before TIC Normalization")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position = "none")

  if(TIC==TRUE){
    ## ------------------ Perform TIC ------------------- ##
    message("Total Ion Count (TIC) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022")
    RowSums <- rowSums(NA_removed_matrix)
    Median_RowSums <- median(RowSums) #This will built the median
    Data_TIC_Pre <- apply(NA_removed_matrix, 2, function(i) i/RowSums) #This is dividing the ion intensity by the total ion count
    Data_TIC <- Data_TIC_Pre*Median_RowSums #Multiplies with the median metabolite intensity
    Data_TIC <- as.data.frame(Data_TIC)

    ## ------------------ QC plot ------------------- ##
    ### After TIC normalization
    log_Data_TIC <- log(Data_TIC) %>% t() %>% as.data.frame()
    medians <- apply(log_Data_TIC, 2, median)
    RLA_data_norm <- log_Data_TIC - medians   # Subtract the medians from each column
    RLA_data_long <- pivot_longer(RLA_data_norm, cols = everything(), names_to = "Group")
    names(RLA_data_long)<- c("Samples", "Intensity")
    for (row in 1:nrow(RLA_data_long)){ # add conditions
      RLA_data_long[row,"Conditions"] <- SettingsFile_Sample[rownames(SettingsFile_Sample) %in%RLA_data_long[row,1],"Conditions"]
    }

    # Create the ggplot boxplot
    RLA_data_norm <- ggplot(RLA_data_long, aes(x = Samples, y = Intensity, color = Conditions)) +
      geom_boxplot() +
      geom_hline(yintercept = 0, color = "red", linetype = "solid") +
      labs(title = "After TIC Normalization")+
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position = "none")

    #Combine Plots
    dev.new()
    norm_plots <- suppressWarnings(gridExtra::grid.arrange(RLA_data_raw+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position = "none"),
                                                           RLA_data_norm+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position = "none"),
                                                           ncol = 2))
    dev.off()
    norm_plots <- ggplot2::ggplot() +theme_minimal()+ annotation_custom(norm_plots)

    ## ------------------ Return ------------------ ##
    Output_list <- list("DF" = list("Data_TIC"=as.data.frame(Data_TIC)),"Plot"=list( "norm_plots"=norm_plots, "RLA_AfterTICNorm"=RLA_data_norm,  "RLA_BeforeTICNorm" = RLA_data_raw ))
    invisible(return(Output_list))
  }else{
    ## ------------------ Return ------------------ ##
    Output_list <- list("Plot"=list("RLA_BeforeTICNorm" = RLA_data_raw ))
    invisible(return(Output_list))
  }
}




################################################################################################
### ### ### PreProcessing helper function: CoRe nomalisation ### ### ###
################################################################################################

#' CoReNorm
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Sample Passed to main function MetaProViz::PreProcessing()
#' @param SettingsInfo Passed to main function MetaProViz::PreProcessing()
#'
#' @keywords Consumption Release Normalisation
#' @noRd
#'

CoReNorm <-function(InputData, SettingsFile_Sample, SettingsInfo){
  ## ------------------ Prepare the data ------------------- ##
  Data_TIC <- InputData
  Data_TIC[is.na(Data_TIC)] <- 0

  ## ------------------ Perform QC ------------------- ##
  Conditions <- SettingsFile_Sample[["Conditions"]]
  CoRe_medias <-  Data_TIC[grep(SettingsInfo[["CoRe_media"]], Conditions),]

  if(dim(CoRe_medias)[1]==1){
    warning("Only 1 CoRe_media sample was found. Thus, the consistency of the CoRe_media samples cannot be checked. It is assumed that the CoRe_media samples are already summed.")
    CoRe_media_df <- CoRe_medias %>% t() %>% as.data.frame()
    colnames(CoRe_medias) <- "CoRe_mediaMeans"
  }else{
    ######################################################################################
    ## ------------------ QC Plots
    PlotList <- list()
    ##-- 1. PCA Media_control
    media_pca_data <- merge(x= SettingsFile_Sample %>% select(SettingsInfo[["Conditions"]]), y= Data_TIC, by=0) %>%
      column_to_rownames("Row.names") %>%
      mutate(Sample_type = case_when(Conditions == SettingsInfo[["CoRe_media"]] ~ "CoRe_media",
                                     TRUE ~ "Sample"))

    media_pca_data[is.na( media_pca_data)] <- 0

    dev.new()
    pca_QC_media <-invisible(MetaProViz::VizPCA(InputData=media_pca_data %>%select(-SettingsInfo[["Conditions"]], -Sample_type),
                                                SettingsInfo= c(color="Sample_type"),
                                                SettingsFile_Sample= media_pca_data,
                                                PlotName = "QC Media_samples",
                                                SaveAs_Plot =  NULL))
    dev.off()

    PlotList[["PCA_CoReMediaSamples"]] <- pca_QC_media[["Plot_Sized"]][["QC Media_samples"]]

    ##-- 2. Metabolite Variance Histogram
    # Coefficient of Variation
    result_df <- apply(CoRe_medias, 2,   function(x) { (sd(x, na.rm =T)/  mean(x, na.rm =T))*100 } ) %>% t()%>% as.data.frame()
    result_df[1, is.na(result_df[1,])]<- 0
    rownames(result_df)[1] <- "CV"

    CutoffCV <- 100
    result_df <- result_df %>% t()%>%as.data.frame() %>% rowwise() %>%
      mutate(HighVar = CV > CutoffCV) %>% as.data.frame()
    rownames(result_df)<- colnames(CoRe_medias)

    # calculate the NAs
    NAvector <- apply(CoRe_medias, 2,  function(x) { (sum(is.na(x))/length(x))*100 })
    result_df$MissingValuePercentage <- NAvector

    cv_result_df <- result_df

    HighVar_metabs <- sum(result_df$HighVar == TRUE)
    if(HighVar_metabs>0){
      message(paste0(HighVar_metabs, " of variables have high variability in the CoRe_media control samples. Consider checking the pooled samples to decide whether to remove these metabolites or not."))
    }

    #Make histogram of CVs
    HistCV <- invisible(ggplot(cv_result_df, aes(CV)) +
                          geom_histogram(aes(y=after_stat(density)), color="black", fill="white")+
                          geom_vline(aes(xintercept=CutoffCV),
                                     color="darkred", linetype="dashed", linewidth=1)+
                          geom_density(alpha=.2, fill="#FF6666") +
                          labs(title="Coefficient of Variation for metabolites of control media samples (no cells)",x="Coefficient of variation (CV)", y = "Frequency")+
                          theme_classic())

    PlotList[["Histogram_CoReMediaCV"]] <- HistCV

    #Make Violin of CVs
    ViolinCV <- invisible(ggplot(cv_result_df, aes(y=CV, x=HighVar, label=row.names(cv_result_df)))+
                            geom_violin(alpha = 0.5 , fill="#FF6666")+
                            geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
                            #geom_point(position = position_jitter(seed = 1, width = 0.2))+
                            geom_text(aes(label=ifelse(cv_result_df$CV>CutoffCV,as.character(row.names(cv_result_df)),'')), hjust=0, vjust=0)+
                            labs(title="Coefficient of Variation for metabolites of control media samples (no cells)",x="Coefficient of variation (CV)", y = "Frequency")+
                            theme_classic())

    PlotList[["CoRe_Media_CV_Violin"]] <- ViolinCV

    ######################################################################################
    ## ------------------ Outlier testing
    if(dim(CoRe_medias)[1]>=3){
      Outlier_data <- CoRe_medias
      Outlier_data <- Outlier_data %>% mutate_all(.funs = ~ FALSE)

      while(HighVar_metabs>0){
        #remove the furthest value from the mean
        if(HighVar_metabs>1){
          max_var_pos <-  CoRe_medias[,result_df$HighVar == TRUE]  %>%
            as.data.frame() %>%
            mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
            summarise_all(.funs = ~ which.max(abs(.)))
        }else{
          max_var_pos <-  CoRe_medias[,result_df$HighVar == TRUE]  %>%
            as.data.frame() %>%
            mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
            summarise_all(.funs = ~ which.max(abs(.)))
          colnames(max_var_pos)<- colnames(CoRe_medias)[result_df$HighVar == TRUE]
        }

        # Remove rows based on positions
        for(i in 1:length(max_var_pos)){
          CoRe_medias[max_var_pos[[i]],names(max_var_pos)[i]] <- NA
          Outlier_data[max_var_pos[[i]],names(max_var_pos)[i]] <- TRUE
        }

        # ReCalculate coefficient of variation for each column in the filtered data
        result_df <- apply(CoRe_medias, 2,   function(x) { sd(x, na.rm =T)/  mean(x, na.rm =T) } ) %>% t()%>% as.data.frame()
        result_df[1, is.na(result_df[1,])]<- 0
        rownames(result_df)[1] <- "CV"

        result_df <- result_df %>% t()%>%as.data.frame() %>% rowwise() %>%
          mutate(HighVar = CV > CutoffCV) %>% as.data.frame()
        rownames(result_df)<- colnames(CoRe_medias)

        HighVar_metabs <- sum(result_df$HighVar == TRUE)
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
      rownames(contingency_data_contframe) <- c("HighVar", "Low_var")

      contingency_data_contframe <- contingency_data_contframe %>% mutate(Total = rowSums(contingency_data_contframe))
      contingency_data_contframe <- rbind(contingency_data_contframe, Total= colSums(contingency_data_contframe))

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
    CoRe_media_df <- as.data.frame(data.frame("CoRe_mediaMeans"=  colMeans(CoRe_medias, na.rm = TRUE)))
  }

  cv_result_df <- rownames_to_column(cv_result_df, "Metabolite")

  ######################################################################################
  ##------------------------ Substract mean (media control) from samples
  message("CoRe data are normalised by substracting mean (blank) from each sample and multiplying with the CoRe_norm_factor")
  ##-- Check CoRe_norm_factor
  if(("CoRe_norm_factor" %in% names(SettingsInfo))){
    CoRe_norm_factor <-   SettingsFile_Sample %>% filter(!!as.name(SettingsInfo[["Conditions"]])!=SettingsInfo[["CoRe_media"]]) %>% select(SettingsInfo[["CoRe_norm_factor"]]) %>%pull()
    if(var(CoRe_norm_factor) ==  0){
      warning("The growth rate or growth factor for normalising the CoRe result, is the same for all samples")
    }
  }else{
    CoRe_norm_factor <- as.numeric(rep(1,dim(SettingsFile_Sample %>% filter(!!as.name(SettingsInfo[["Conditions"]])!=SettingsInfo[["CoRe_media"]]))[1]))
  }

  # Remove CoRe_media samples from the data
  Data_TIC <- merge(SettingsFile_Sample, Data_TIC, by="row.names")%>%
    filter(!!as.name(SettingsInfo[["Conditions"]])!=SettingsInfo[["CoRe_media"]])%>%
    column_to_rownames("Row.names")%>%
    select(-1:-ncol(SettingsFile_Sample))

  Data_TIC_CoReNorm_Media <- as.data.frame(t( apply(t(Data_TIC),2, function(i) i-CoRe_media_df$CoRe_mediaMeans)))  #Subtract from each sample the CoRe_media mean
  Data_TIC_CoReNorm <- as.data.frame(apply(Data_TIC_CoReNorm_Media, 2, function(i) i*CoRe_norm_factor))

  #Remove CoRe_media samples from the data
  #Input_SettingsFile <- Input_SettingsFile[Input_SettingsFile$Conditions!="CoRe_media",]
  #Conditions <- Conditions[!Conditions=="CoRe_media"]

  ######################################################################################
  ##------------------------ Return Plots and Data
  DF_list <- list("CV_CoRe_blank" = cv_result_df, "Contigency_table_CoRe_blank" = contingency_data_contframe, "Core_Norm" = Data_TIC_CoReNorm)

  #Return
  Output_list <- list("DF"= DF_list,"Plot"=PlotList)
  invisible(return(Output_list))
}



################################################################################################
### ### ### PreProcessing helper function: Outlier detection ### ### ###
################################################################################################

#' OutlierDetection
#'
#' @param InputData Passed to main function MetaProViz::PreProcessing()
#' @param SettingsFile_Sample Passed to main function MetaProViz::PreProcessing()
#' @param SettingsInfo Passed to main function MetaProViz::PreProcessing()
#' @param CoRe Passed to main function MetaProViz::PreProcessing()
#' @param HotellinsConfidence Passed to main function MetaProViz::PreProcessing()
#'
#' @keywords Hotellins T2 outlier detection
#' @noRd
#'

OutlierDetection <-function(InputData, SettingsFile_Sample, SettingsInfo, CoRe, HotellinsConfidence){
  # Message:
  message("Identification of outlier samples is performed using Hotellin's T2 test to define sample outliers in a mathematical way (Confidence = 0.99 ~ p.val < 0.01) REF: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3), 360–378, doi:https://doi.org/10.1214/aoms/1177732979.")
  message(paste("HotellinsConfidence value selected:", HotellinsConfidence))

  # Load the data:
  data_norm <- InputData%>%
    mutate_all(~ replace(., is.nan(.), 0))
  data_norm[is.na(data_norm)] <- 0 #replace NA with 0

  if(CoRe==TRUE){
    Conditions <- SettingsFile_Sample[[SettingsInfo[["Conditions"]]]][!SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]]]
  }else{
    Conditions <- SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]
  }


  # Prepare the lists to store the results:
  Outlier_filtering_loop = 10 #Here we do 10 rounds of hotelling filtering
  sample_outliers <- list()
  scree_plot_list <- list()
  outlier_plot_list <- list()
  metabolite_zero_var_total_list <- list()
  zero_var_metab_warning = FALSE


  #################################################
  ##--------- Perform Outlier testing:
  for(loop in 1:Outlier_filtering_loop){
    ##--- Zero variance metabolites
    metabolite_var <- as.data.frame(apply(data_norm, 2, var) %>% t()) # calculate each metabolites variance
    metabolite_zero_var_list <- list(colnames(metabolite_var)[which(metabolite_var[1,]==0)]) # takes the names of metabolites with zero variance and puts them in list

    if(sum(metabolite_var[1,]==0)==0){
      metabolite_zero_var_total_list[loop] <- 0
    }else if(sum(metabolite_var[1,]==0)>0){
      metabolite_zero_var_total_list[loop] <- metabolite_zero_var_list
      zero_var_metab_warning = TRUE # This is used later to print and save the zero variance metabolites if any are found.
    }

    for(metab in metabolite_zero_var_list){  # Remove the metabolites with zero variance from the data to do PCA
      data_norm <- data_norm %>% select(-all_of(metab))
    }

    ##---  PCA
    PCA.res <- prcomp(data_norm, center =  TRUE, scale. =  TRUE)
    outlier_PCA_data <- data_norm
    outlier_PCA_data$Conditions <- Conditions

    dev.new()
    pca_outlier <-invisible(MetaProViz::VizPCA(InputData=data_norm,
                                               SettingsInfo= c(color=SettingsInfo[["Conditions"]]),
                                               SettingsFile_Sample= outlier_PCA_data,
                                               PlotName = paste("PCA outlier test filtering round ",loop),
                                               SaveAs_Plot =  NULL))

    if(loop==1){
      pca_outlierloop1 <- pca_outlier[["Plot_Sized"]][[1]]
    }
    #suppressMessages(plot(pca_outlier[["Plot_Sized"]][[1]]))
    outlier_plot_list[[paste("PCA_round",loop,sep="")]] <- pca_outlier[["Plot_Sized"]][[1]]
    dev.off()

    ##--- Scree plot
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

    if(loop==1){
      scree_outlierloop1 <-screeplot
    }
    dev.new()
    #plot(screeplot)
    outlier_plot_list[[paste("ScreePlot_round",loop,sep="")]] <- screeplot # save plot
    dev.off()

    ##--- HotellingT2 test for outliers
    data_hot <- as.matrix(PCA.res$x[,1:npcs])
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
      scale_y_continuous(breaks = sort(c(ggplot_build(HotellingT2plot)$layout$panel_ranges[[1]]$y.major_source, c(limits[,1],limits[,2]))))

    HotellingT2plot <- HotellingT2plot + theme_classic()
    HotellingT2plot <- HotellingT2plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    HotellingT2plot <- HotellingT2plot + ggtitle(paste("Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"))
    HotellingT2plot <- HotellingT2plot + scale_linetype_discrete(name = LegendTitle,)
    HotellingT2plot <- HotellingT2plot + theme(plot.title = element_text(size = 13))+#, face = "bold")) +
      theme(axis.text = element_text(size = 12))

    if(loop==1){
      hotel_outlierloop1 <- HotellingT2plot
    }
    dev.new()
    #plot(HotellingT2plot)
    outlier_plot_list[[paste("HotellingsPlot_round",loop,sep="")]] <- HotellingT2plot
    dev.off()

    a<- loop
    if(CoRe==TRUE){
      a<- paste0(a,"_CoRe")
    }

    if(length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 0){ # loop for outliers until no outlier is detected
      data_norm <- data_norm
      break
    }else if(length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 1){
      data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]# filter the selected outliers from the data
      Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

      # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
      hotelling_qcc[["violations"]][["beyond.limits"]][1] <-   rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][1]]
      sample_outliers[loop] <- list(hotelling_qcc[["violations"]][["beyond.limits"]])
    }else{
      data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
      Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

      # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
      sm_out <- c() # list of outliers samples
      for (i in 1:length(hotelling_qcc[["violations"]][["beyond.limits"]])){
        sm_out <-  append(sm_out, rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][i]])
      }
      sample_outliers[loop] <- list(sm_out )
    }
  }

  #################################################
  ##-- Print Outlier detection results about samples and metabolites
  if(length(sample_outliers) > 0){   # Print outlier samples
    message("There are possible outlier samples in the data") #This was a warning
    for (i in 1:length(sample_outliers)  ){
      message("Filtering round ",i ," Outlier Samples: ", paste( head(sample_outliers[[i]]) ," "))
    }
  }else{message("No sample outliers were found")}


  ##--  Print Zero variance metabolites
  zero_var_metab_export_df <- data.frame(1,2)
  names(zero_var_metab_export_df) <- c("Filtering round","Metabolite")

  if(zero_var_metab_warning==TRUE){
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

  #############################################
  ##---- 1. Make Output DF
  total_outliers <- hash::hash() # make a dictionary
  if(length(sample_outliers) > 0){ # Create columns with outliers to merge to output dataframe
    for (i in 1:length(sample_outliers)  ){
      total_outliers[[paste("Outlier_filtering_round_",i, sep = "")]] <- sample_outliers[i]
    }
  }

  data_norm_filtered_full <- as.data.frame(replace(InputData, InputData==0, NA))

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
  data_norm_filtered_full <- merge(SettingsFile_Sample, data_norm_filtered_full,  by = 0) # add the design in the output df (merge by rownames/sample names)
  rownames(data_norm_filtered_full) <- data_norm_filtered_full$Row.names
  data_norm_filtered_full$Row.names <- c()

  ##-- 2.  Quality Control (QC) PCA
  MetaData_Sample <- data_norm_filtered_full %>%
    mutate(Outliers = case_when(Outliers == "no" ~ 'no',
                                Outliers == "Outlier_filtering_round_1" ~ ' Outlier_filtering_round = 1',
                                Outliers == "Outlier_filtering_round_2" ~ ' Outlier_filtering_round = 2',
                                Outliers == "Outlier_filtering_round_3" ~ ' Outlier_filtering_round = 3',
                                Outliers == "Outlier_filtering_round_4" ~ ' Outlier_filtering_round = 4',
                                TRUE ~ 'Outlier_filtering_round = or > 5'))
  MetaData_Sample$Outliers <- relevel(as.factor(MetaData_Sample$Outliers), ref="no")

  # 1. Shape Outliers
  if(length(sample_outliers)>0){
    dev.new()
    pca_QC <-invisible(MetaProViz::VizPCA(InputData=as.data.frame(InputData)%>%select(-zero_var_metab_export_df$Metabolite),
                                          SettingsInfo= c(color=SettingsInfo[["Conditions"]], shape = "Outliers"),
                                          SettingsFile_Sample= MetaData_Sample ,
                                          PlotName = "Quality Control PCA Condition clustering and outlier check",
                                          SaveAs_Plot =  NULL))
    dev.off()
    outlier_plot_list[["QC_PCA_and_Outliers"]] <- pca_QC[["Plot_Sized"]][[1]]
  }

  # 2. Shape Biological replicates
  if(SettingsInfo[["Biological_Replicates"]] %in% colnames(SettingsFile_Sample)){
    dev.new()
    pca_QC_repl <-invisible(MetaProViz::VizPCA(InputData=as.data.frame(InputData)%>%select(-zero_var_metab_export_df$Metabolite),
                                               SettingsInfo= c(color=SettingsInfo[["Conditions"]], shape = SettingsInfo[["Biological_Replicates"]]),
                                               SettingsFile_Sample= MetaData_Sample,
                                               PlotName =  "Quality Control PCA replicate spread check",
                                               SaveAs_Plot =  NULL))
    dev.off()

    outlier_plot_list[["QC_PCA_Replicates"]] <- pca_QC_repl[["Plot_Sized"]][[1]]
  }



  #############################################
  ##--- Save and Return plots and DFs
  DF_list <- list("Zero_variance_metabolites_CoRe" = zero_var_metab_export_df, "data_outliers" = data_norm_filtered_full)

  #Return
  Output_list <- list("DF"= DF_list,"Plot"=outlier_plot_list)
  invisible(return(Output_list))
}
