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



###################################################
### ### ### Metabolomics pre-processing ### ### ###
###################################################

#' Modularised Normalization: 80%-filtering rule, total-ion count normalization, missing value imputation and Outlier Detection: HotellingT2.
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  or Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "BiologicalReplicates" including numerical values. For CoRe = TRUE a CoRe_norm_factor = "Columnname_Input_SettingsFile" and CoRe_media = "Columnname_Input_SettingsFile", have to also be added. Column CoRe_norm_factor is used for normalization and CoRe_media is used to specify the name of the media controls in the Conditions.
#' @param FeatureFilt \emph{Optional: }If NULL, no feature filtering is performed. If set to "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = "Standard"}
#' @param FeatureFilt_Value \emph{Optional: } Percentage of feature filtering. \strong{Default = 0.8}
#' @param TIC \emph{Optional: } If TRUE, Total Ion Count normalization is performed. \strong{Default = TRUE}
#' @param MVI \emph{Optional: } If TRUE, Missing Value Imputation (MVI) based on half minimum is performed \strong{Default = TRUE}
#' @param MVI_Percentage \emph{Optional: } Percentage 0-100 of imputed value based on the minimum value. \strong{Default = 50}
#' @param HotellinsConfidence \emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. Must be numeric.\strong{Default = 0.99}
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed and the CoRe value will be calculated. Please consider providing a Normalisation factor column called "CoRe_norm_factor" in your "Input_SettingsFile" DF, where the column "Conditions" matches. The normalisation factor must be a numerical value obtained from growth rate that has been obtained from a growth curve or growth factor that was obtained by the ratio of cell count/protein quantification at the start point to cell count/protein quantification at the end point.. Additionally control media samples have to be available in the "Input" DF and defined as "CoRe_media" samples in the "Conditions" column in the "Input_SettingsFile" DF. \strong{Default = FALSE}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. If set to NULL, plots are not saved. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } Select the file type of output table. Options are "csv", "xlsx", "txt". If set to NULL, plots are not saved. \strong{Default = "csv"}
#' @param PrintPlot  \emph{Optional: } If TRUE prints an overview of resulting plots. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List with two elements: DF (including all output tables generated) and Plot (including all plots generated)
#'
#' @examples
#' Intra <- MetaProViz::ToyData("IntraCells_Raw")
#' Res <- MetaProViz::PreProcessing(InputData=Intra[-c(49:58) ,-c(1:3)],
#'                                  SettingsFile_Sample=Intra[-c(49:58) , c(1:3)],
#'                                  SettingsInfo = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates"))
#'
#' @keywords 80  percent filtering rule, Missing Value Imputation, Total Ion Count normalization, PCA, HotellingT2, multivariate quality control charts
#'
#' @importFrom dplyr mutate_all
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
#'
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
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`
  CheckInput(InputData=InputData,
                          SettingsFile_Sample=SettingsFile_Sample,
                          SettingsFile_Metab=NULL,
                          SettingsInfo= SettingsInfo,
                          SaveAs_Plot=SaveAs_Plot,
                          SaveAs_Table=SaveAs_Table,
                          CoRe=CoRe,
                          PrintPlot= PrintPlot)

  # HelperFunction `CheckInput` Specific
  CheckInput_PreProcessing(InputData=InputData,
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
    Folder <- SavePath(FolderName= "Processing",
                                    FolderPath=FolderPath)

    SubFolder_P <- file.path(Folder, "PreProcessing")
    if (!dir.exists(SubFolder_P)) {dir.create(SubFolder_P)}
  }

  ## ------------------ Prepare the data ------------------- ##
  #InputData files:
  InputData <-as.data.frame(InputData)%>%
    dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .))#Make sure all 0 are changed to NAs

  InputData <- as.data.frame(dplyr::mutate_all(as.data.frame(InputData), function(x) as.numeric(as.character(x))))

  ###################################################################################################################################
  ## ------------------ 1. Feature filtering ------------------- ##
  if(is.null(FeatureFilt)==FALSE){
    InputData_Filtered <- FeatureFiltering(InputData=InputData,
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
    MVIRes<- MVImputation(InputData=InputData_Filt,
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
    TICRes_List <- TICNorm(InputData=MVIRes,
                                        SettingsFile_Sample=SettingsFile_Sample,
                                        SettingsInfo=SettingsInfo,
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
    RLAPlot_List <- TICNorm(InputData=MVIRes,
                                         SettingsFile_Sample=SettingsFile_Sample,
                                         SettingsInfo=SettingsInfo,
                                         TIC=TIC)
    PlotList <- list()
    PlotList[["RLAPlot"]] <- RLAPlot_List[["Plot"]][["RLA_BeforeTICNorm"]]
  }

  ## ------------------ 4. CoRe media QC (blank) and normalization ------------------- ##
  if(CoRe ==TRUE){
   data_CoReNorm <- CoReNorm(InputData= TICRes,
                                          SettingsFile_Sample=SettingsFile_Sample,
                                          SettingsInfo=SettingsInfo)

    TICRes <- data_CoReNorm[["DF"]][["Core_Norm"]]
  }

  # ------------------ Final Output:
  data_norm <- TICRes %>% as.data.frame()

  ###################################################################################################################################
  ## ------------------ Sample outlier identification ------------------- ##
  OutlierRes <-  OutlierDetection(InputData= data_norm,
                                               SettingsFile_Sample=SettingsFile_Sample,
                                               SettingsInfo=SettingsInfo,
                                               CoRe=CoRe,
                                               HotellinsConfidence=HotellinsConfidence)

  ###################################################################################################################################
  ## ------------------ Return ------------------- ##
  ## ---- DFs
  if(is.null(FeatureFilt)==FALSE){#Add metabolites that where removed as part of the feature filtering
    if(length(InputData_Filtered[["RemovedMetabolites"]])==0){
      DFList <- list("InputData_RawData"= merge(as.data.frame(SettingsFile_Sample), as.data.frame(InputData), by="row.names")%>% tibble::column_to_rownames("Row.names"),
                     "Filtered_metabolites"= as.data.frame(list(FeatureFiltering = c(FeatureFilt),
                                                            FeatureFilt_Value = c(FeatureFilt_Value),
                                                            RemovedMetabolites = c("None"))),
                     "Preprocessing_output"=OutlierRes[["DF"]][["data_outliers"]])
    }else{
      DFList <- list("InputData_RawData"= merge(as.data.frame(SettingsFile_Sample), as.data.frame(InputData), by="row.names")%>% tibble::column_to_rownames("Row.names"),
                     "Filtered_metabolites"= as.data.frame(list(FeatureFiltering = rep(FeatureFilt, length(InputData_Filtered[["RemovedMetabolites"]])),
                                                           FeatureFilt_Value = rep(FeatureFilt_Value, length(InputData_Filtered[["RemovedMetabolites"]])),
                                                           RemovedMetabolites = InputData_Filtered[["RemovedMetabolites"]])),
                     "Preprocessing_output"=OutlierRes[["DF"]][["data_outliers"]])
    }
  }else{
    DFList <- list("InputData_RawData"= merge(as.data.frame(SettingsFile_Sample), as.data.frame(InputData), by="row.names")%>% tibble::column_to_rownames("Row.names"), "Preprocessing_output"=OutlierRes[["DF"]][["data_outliers"]])
  }

  if(CoRe ==TRUE){
     DFList_CoRe <- list( "CV_CoRe_blank"= data_CoReNorm[["DF"]][["CV_CoRe_blank"]],"Variation_ContigencyTable_CoRe_blank"=data_CoReNorm[["DF"]][["Contigency_table_CoRe_blank"]])
     DFList <- c(DFList, DFList_CoRe)
  }

  ## ---- Plots
  if(TIC==TRUE){
    PlotList <- c(TICRes_List[["Plot"]], OutlierRes[["Plot"]])
  }else{
    PlotList <- c(RLAPlot_List[["Plot"]], OutlierRes[["Plot"]])
  }

  if(CoRe ==TRUE){
    PlotList <- c(PlotList , data_CoReNorm[["Plot"]])
  }

  Res_List <- list("DF"= DFList ,"Plot" =PlotList)

  # Save Plots and DFs
  #As row names are not saved we need to make row.names to column for the DFs that needs this:
  DFList[["InputData_RawData"]] <- DFList[["InputData_RawData"]]%>%tibble::rownames_to_column("Code")
  DFList[["Preprocessing_output"]] <- DFList[["Preprocessing_output"]]%>%tibble::rownames_to_column("Code")

  suppressMessages(suppressWarnings(
    SaveRes(InputList_DF=DFList,
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
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return DF with the merged analytical replicates
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")
#' Res <- ReplicateSum(InputData=Intra[-c(49:58) ,-c(1:3)],
#'                                 SettingsFile_Sample=Intra[-c(49:58) , c(1:3)],
#'                                 SettingsInfo = c(Conditions="Conditions", Biological_Replicates="Biological_Replicates", Analytical_Replicates="Analytical_Replicates"))
#'
#' @keywords Analytical Replicate Merge
#'
#' @importFrom dplyr mutate_all summarise_all select rename ungroup group_by
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom rlang !! :=
#' @importFrom tidyr unite
#'
#' @export
#'
ReplicateSum <- function(InputData,
                         SettingsFile_Sample,
                         SettingsInfo = c(Conditions="Conditions", Biological_Replicates="Biological_Replicates", Analytical_Replicates="Analytical_Replicates"),
                         SaveAs_Table = "csv",
                         FolderPath = NULL){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`
  CheckInput(InputData=InputData,
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
    Folder <- SavePath(FolderName= "Processing",
                                    FolderPath=FolderPath)
    SubFolder <- file.path(Folder, "ReplicateSum")
    if (!dir.exists(SubFolder)) {dir.create(SubFolder)}
  }

  ## ------------  Load data and process  ----------- ##
  Input <- merge(x= SettingsFile_Sample%>% dplyr::select(!!SettingsInfo[["Conditions"]], !!SettingsInfo[["Biological_Replicates"]], !!SettingsInfo[["Analytical_Replicates"]]),
                 y= InputData,
                 by="row.names")%>%
    tibble::column_to_rownames("Row.names")%>%
    dplyr::rename("Conditions"=SettingsInfo[["Conditions"]],
                  "Biological_Replicates"=SettingsInfo[["Biological_Replicates"]],
                  "Analytical_Replicates"=SettingsInfo[["Analytical_Replicates"]])

  # Make the replicate Sums
  Input_data_numeric_summed <- as.data.frame(Input %>%
                                             dplyr::group_by(Biological_Replicates, Conditions) %>%
                                             dplyr::summarise_all("mean") %>% dplyr::select(-Analytical_Replicates))

  # Make a number of merged replicates column
  nReplicates <-  Input %>%
    dplyr::group_by(Biological_Replicates, Conditions) %>%
    dplyr::summarise_all("max") %>%
    dplyr::ungroup() %>%
    dplyr::select(Analytical_Replicates, Biological_Replicates, Conditions) %>%
    dplyr::rename("n_AnalyticalReplicates_Summed "= "Analytical_Replicates")

  Input_data_numeric_summed <- merge(nReplicates,Input_data_numeric_summed, by = c("Conditions","Biological_Replicates"))%>%
    tidyr::unite(UniqueID, c("Conditions","Biological_Replicates"), sep="_", remove=FALSE)%>% # Create a uniqueID
    tibble::column_to_rownames("UniqueID")# set UniqueID to rownames

  #--------------- return ------------------##
  SaveRes(InputList_DF=list("Sum_AnalyticalReplicates"=Input_data_numeric_summed%>%tibble::rownames_to_column("Code")),
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

#' Find metabolites with high variqability across total pool samples
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected. Can be either a full dataset or a dataset with only the pool samples.
#' @param SettingsFile_Sample  \emph{Optional: } DF which contains information about the samples when a full dataset is inserted as Input_data. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), has to exist.\strong{Default = NULL}
#' @param SettingsInfo  \emph{Optional: } NULL or Named vector including the Conditions and PoolSample information (Name of the Conditions column and Name of the pooled samples in the Conditions in the Input_SettingsFile)  : c(Conditions="ColumnNameConditions, PoolSamples=NamePoolCondition. If no Conditions is added in the Input_SettingsInfo, it is assumed that the conditions column is named 'Conditions' in the Input_SettingsFile. ). \strong{Default = NULL}
#' @param CutoffCV \emph{Optional: } Filtering cutoff for high variance metabolites using the Coefficient of Variation. \strong{Default = 30}
#' @param SaveAs_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf or NULL. \strong{Default = svg}
#' @param SaveAs_Table \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt", ot NULL \strong{default: "csv"}
#' @param PrintPlot \emph{Optional: } If TRUE prints an overview of resulting plots. \strong{Default = TRUE}
#' @param FolderPath \emph{Optional:} Path to the folder the results should be saved at. \strong{default: NULL}
#'
#' @return List with two elements: DF (including input and output table) and Plot (including all plots generated)
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")
#' Res <- PoolEstimation(InputData=Intra[ ,-c(1:3)],
#'                                 SettingsFile_Sample=Intra[ , c(1:3)],
#'                                 SettingsInfo = c(PoolSamples = "Pool", Conditions="Conditions"))
#'
#' @keywords Coefficient of Variation, high variance metabolites
#'
#' @importFrom dplyr case_when select rowwise mutate ungroup
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom logger log_info log_trace
#'
#' @export
#'
PoolEstimation <- function(InputData,
                           SettingsFile_Sample = NULL,
                           SettingsInfo = NULL,
                           CutoffCV = 30,
                           SaveAs_Plot = "svg",
                           SaveAs_Table = "csv",
                           PrintPlot=TRUE,
                           FolderPath = NULL){

  ## ------------ Create log file ----------- ##
  MetaProViz_Init()

  logger::log_info('Starting pool estimation.')
  ## ------------------ Check Input ------------------- ##
  # HelperFunction `CheckInput`
  CheckInput(InputData=InputData,
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
      if(SettingsInfo[["PoolSamples"]] %in% SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] == FALSE ){
        stop("You have chosen PoolSamples = ",paste(SettingsInfo[["PoolSamples"]] ), ", ", paste(SettingsInfo[["PoolSamples"]] )," was not found in SettingsFile_Sample as sample condition. Please insert the name of the pool samples as stated in the Conditions column of the SettingsFile_Sample."   )
      }
    }
  }

  if(is.numeric(CutoffCV)== FALSE | CutoffCV < 0){
    stop("Check input. The selected CutoffCV value should be a positive numeric value.")
  }

  ## ------------------  Create output folders  and path ------------------- ##
  if(is.null(SaveAs_Plot)==FALSE |is.null(SaveAs_Table)==FALSE ){
    Folder <- SavePath(FolderName= "Processing",
                                    FolderPath=FolderPath)

    SubFolder <- file.path(Folder, "PoolEstimation")
    logger::log_info('Selected output directory: `%s`.', SubFolder)
    if (!dir.exists(SubFolder)) {
      logger::log_trace('Creating directory: `%s`.', SubFolder)
      dir.create(SubFolder)
    }
  }

  ## ------------------ Prepare the data ------------------- ##
  #InputData files:
  if(is.null(SettingsFile_Sample)==TRUE){
    PoolData <- InputData
    PoolData[PoolData == 0] <- NA
  }else{
    PoolData <- InputData[SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["PoolSamples"]],]
    PoolData[PoolData == 0] <- NA
  }


  ###################################################################################################################################
  ## ------------------ Coefficient of Variation ------------------- ##
  logger::log_trace('Calculating coefficient of variation.')
  result_df <- apply(PoolData, 2,  function(x) { (sd(x, na.rm =T)/  mean(x, na.rm =T))*100 }  ) %>% t()%>% as.data.frame()
  rownames(result_df)[1] <- "CV"

  NAvector <- apply(PoolData, 2,  function(x) {(sum(is.na(x))/length(x))*100 })# Calculate the NAs

  # Create Output DF
  result_df_final <- result_df %>%
    t()%>% as.data.frame() %>% dplyr::rowwise() %>%
    dplyr::mutate(HighVar = CV > CutoffCV) %>% as.data.frame()

  result_df_final$MissingValuePercentage <- NAvector

  rownames(result_df_final)<- colnames(InputData)
  result_df_final_out <- tibble::rownames_to_column(result_df_final,"Metabolite" )

  # Remove Metabolites from InputData based on CutoffCV
  logger::log_trace('Applying CV cut-off.')
  if(is.null(SettingsFile_Sample)==FALSE){
      unstable_metabs <- rownames(result_df_final)[result_df_final[["HighVar_Metabs"]]]
      if(length(unstable_metabs)>0){
        filtered_Input_data <- InputData %>% dplyr::select(!unstable_metabs)
      }else{
        filtered_Input_data <- NULL
  }
      }else{
    filtered_Input_data <- NULL
  }

  ## ------------------ QC plots ------------------- ##
  # Start QC plot list
  logger::log_info('Plotting QC plots.')
  PlotList <- list()

  # 1. Pool Sample PCA
  logger::log_trace('Pool sample PCA.')
  dev.new()
  if(is.null(SettingsFile_Sample)==TRUE){
    pca_data <- PoolData
    pca_QC_pool <-invisible(VizPCA(InputData=pca_data,
                                               PlotName = "QC Pool samples",
                                               SaveAs_Plot =  NULL))
  }else{
    pca_data <- merge(SettingsFile_Sample %>% dplyr::select(SettingsInfo[["Conditions"]]), InputData, by=0) %>%
      tibble::column_to_rownames("Row.names") %>%
      dplyr::mutate(Sample_type = dplyr::case_when(.data[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["PoolSamples"]] ~ "Pool",
                                     TRUE ~ "Sample"))

    pca_QC_pool <-invisible(VizPCA(InputData=pca_data %>%dplyr::select(-all_of(SettingsInfo[["Conditions"]]), -Sample_type),
                                               SettingsInfo= c(color="Sample_type"),
                                               SettingsFile_Sample= pca_data,
                                               PlotName = "QC Pool samples",
                                               SaveAs_Plot =  NULL))
  }
  dev.off()
  PlotList [["PCAPlot_PoolSamples"]] <- pca_QC_pool[["Plot_Sized"]][["Plot_Sized"]]


  # 2. Histogram of CVs
  logger::log_trace('CV histogram.')
  HistCV <-suppressWarnings(invisible(ggplot(result_df_final_out, aes(CV)) +
                        geom_histogram(aes(y=after_stat(density)), color="black", fill="white")+
                        geom_vline(aes(xintercept=CutoffCV),
                                   color="darkred", linetype="dashed", size=1)+
                        geom_density(alpha=.2, fill="#FF6666") +
                        labs(title="CV for metabolites of Pool samples",x="Coefficient of variation (CV%)", y = "Frequency")+
                        theme_classic()))

  HistCV_Sized <- plotGrob_Processing(InputPlot =  HistCV, PlotName= "CV for metabolites of Pool samples", PlotType= "Hist")
  PlotList [["Histogram_CV-PoolSamples"]] <- HistCV_Sized

  # 2. ViolinPlot of CVs
  logger::log_trace('CV violin plot.')
  #Make Violin of CVs
  Plot_cv_result_df <- result_df_final_out %>%
    dplyr::mutate(HighVar = ifelse((CV > CutoffCV)==TRUE, paste("> CV", CutoffCV, sep=""), paste("< CV", CutoffCV, sep="")))

  ViolinCV <- invisible(ggplot( Plot_cv_result_df, aes(y=CV, x=HighVar, label=Plot_cv_result_df$Metabolite))+
                          geom_violin(alpha = 0.5 , fill="#FF6666")+
                          geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
                          ggrepel::geom_text_repel(aes(label = ifelse(Plot_cv_result_df$CV > CutoffCV,
                                                                      as.character(Plot_cv_result_df$Metabolite), '')),
                                                   hjust = 0, vjust = 0,
                                                   box.padding = 0.5, # space between text and point
                                                   point.padding = 0.5, # space around points
                                                   max.overlaps = Inf) + # allow for many labels
                          labs(title="CV for metabolites of Pool samples",x="Metabolites", y = "Coefficient of variation (CV%)")+
                          theme_classic())

  ViolinCV_Sized <- plotGrob_Processing(InputPlot = ViolinCV, PlotName= "CV for metabolites of Pool samples", PlotType= "Violin")

  PlotList [["ViolinPlot_CV-PoolSamples"]] <- ViolinCV_Sized

  ###################################################################################################################################
  ## ------------------ Return and Save ------------------- ##
  #Save
  logger::log_info('Preparing saved and returned data.')
  if(is.null(filtered_Input_data)==FALSE){
    DF_list <- list("InputData" = InputData, "Filtered_InputData" = filtered_Input_data, "CV" = result_df_final_out )
  }else{
    DF_list <- list("InputData" = InputData, "CV" = result_df_final_out)
  }
  ResList <- list("DF"= DF_list,"Plot"=PlotList)

  #Save
  DF_list[["InputData"]]<-  DF_list[["InputData"]]%>%tibble::rownames_to_column("Code")

  logger::log_info(
    'Saving results: [SaveAs_Table=%s, SaveAs_Plot=%s, FolderPath=%s].',
    SaveAs_Table,
    SaveAs_Plot,
    SubFolder
  )
  SaveRes(InputList_DF=DF_list,
                      InputList_Plot = PlotList,
                      SaveAs_Table=SaveAs_Table,
                      SaveAs_Plot=SaveAs_Plot,
                      FolderPath= SubFolder,
                      FileName= "PoolEstimation",
                      CoRe=FALSE,
                      PrintPlot=PrintPlot)

  #Return
  logger::log_info('Finished pool estimation.')
  invisible(return(ResList))
}


################################################################################################
### ### ### PreProcessing helper function: FeatureFiltering ### ### ###
################################################################################################

#' FeatureFiltering
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile"). Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "BiologicalReplicates" including numerical values. For CoRe = TRUE add CoRe_media = "Columnname_Input_SettingsFile", which specifies the name of the media controls in the Conditions.
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed.Should not be normalised to media blank. Provide information about control media sample names via SettingsInfo "CoRe_media" samples. \strong{Default = FALSE}
#' @param FeatureFilt \emph{Optional: } If NULL, no feature filtering is performed. If set to "Standard" then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Input_SettingsFile input file including the individual conditions you want to apply the filtering to (Yang, J et al., 2015). \strong{Default = Modified}
#' @param FeatureFilt_Value \emph{Optional: } Percentage of feature filtering. \strong{Default = 0.8}
#'
#' @return List with two elements: filtered matrix  and features filtered
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")
#' Res <- FeatureFiltering(InputData=Intra[-c(49:58), -c(1:3)]%>% dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)),
#'                                      SettingsFile_Sample=Intra[-c(49:58), c(1:3)],
#'                                      SettingsInfo = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates"))
#'
#' @keywords feature filtering or modified feature filtering
#'
#' @importFrom dplyr filter mutate_all
#' @importFrom magrittr %>% %<>%
#' @importFrom logger log_info log_trace
#'
#' @noRd
#'
FeatureFiltering <-function(InputData,
                            SettingsFile_Sample,
                            SettingsInfo,
                            CoRe=FALSE,
                            FeatureFilt="Modified",
                            FeatureFilt_Value=0.8){
  ## ------------ Create log file ----------- ##
  MetaProViz_Init()


  ## ------------------ Prepare the data ------------------- ##
  feat_filt_data <- as.data.frame(replace(InputData, InputData==0, NA))

  if(CoRe== TRUE){ # remove CoRe_media samples for feature filtering
    feat_filt_data <- feat_filt_data %>% dplyr::filter(!SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] ==SettingsInfo[["CoRe_media"]])
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
  }else if(FeatureFilt ==  "Standard"){
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
  filtered_matrix <- as.data.frame(dplyr::mutate_all(as.data.frame(filtered_matrix), function(x) as.numeric(as.character(x))))

  Filtered_results <- list("DF"= filtered_matrix , "RemovedMetabolites" = features_filtered)
  invisible(return(Filtered_results))
}



################################################################################################
### ### ### PreProcessing helper function: Missing Value imputation ### ### ###
################################################################################################

#' Missing Value Imputation using half minimum value
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile", CoRe_media = "Columnname_Input_SettingsFile"). Column "Conditions" with information about the sample conditions, Column "BiologicalReplicates" including numerical values and Column "Columnname_Input_SettingsFile" is used to specify the name of the media controls in the Conditions.
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed. Should not be normalised to media blank. Provide information about control media sample names via SettingsInfo "CoRe_media" samples.\strong{Default = FALSE}
#' @param MVI_Percentage \emph{Optional: } Percentage 0-100 of imputed value based on the minimum value. \strong{Default = 50}
#'
#' @return DF with imputed values
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")
#' Res <- MVImputation(InputData=Intra[-c(49:58), -c(1:3)]%>% dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)),
#'                                  SettingsFile_Sample=Intra[-c(49:58), c(1:3)],
#'                                  SettingsInfo = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates"))
#'
#' @keywords Half minimum missing value imputation
#'
#' @importFrom dplyr select mutate group_by filter
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble column_to_rownames
#' @importFrom logger log_info log_trace
#'
#' @noRd
#'
MVImputation <-function(InputData,
                        SettingsFile_Sample,
                        SettingsInfo,
                        CoRe=FALSE,
                        MVI_Percentage=50){
  ## ------------------ Prepare the data ------------------- ##
  filtered_matrix <- InputData
  filtered_matrix[filtered_matrix == 0] <- NA

  ## ------------------ Perform MVI ------------------ ##
  # Do MVI for the samples
  message("Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0")

  if(CoRe==TRUE){#remove blank samples
    NA_removed_matrix <- filtered_matrix%>% dplyr::filter(!SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]])

  }else{
    NA_removed_matrix <- filtered_matrix %>% as.data.frame()
  }

  for (feature  in colnames(NA_removed_matrix)){
    feature_data <- merge(NA_removed_matrix[feature] , SettingsFile_Sample %>% dplyr::select(Conditions), by= 0)
    feature_data <- tibble::column_to_rownames(feature_data, "Row.names")

    imputed_feature_data <- feature_data %>%
      dplyr::group_by(Conditions) %>%
      dplyr::mutate(across(all_of(feature), ~{
        if(all(is.na(.))) {
          message("For some conditions all measured samples are NA for " , feature, ". Hence we can not perform half-minimum value imputation per condition for this metabolite and will assume it is a true biological 0 in those cases.")
          return(0)  # Return NA if all values are missing
        } else {
          return(replace(., is.na(.), min(., na.rm = TRUE)*(MVI_Percentage/100)))
        }
      }))

    NA_removed_matrix[[feature]] <- imputed_feature_data[[feature]]
  }

  if(CoRe==TRUE){
    replaceNAdf <- filtered_matrix%>% dplyr::filter(SettingsFile_Sample[[SettingsInfo[["Conditions"]]]] == SettingsInfo[["CoRe_media"]])

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
        message("Metabolites with only NAs (=100%) in Control_media samples are: ",paste(names(OnlyNA_metabs), collapse = ", "), ". Those NAs are set zero as we consider them true zeros")
      }
    }

    # if all values are NA set to 0
    replaceNAdf_zero <- as.data.frame(lapply(replaceNAdf, function(x) if(all(is.na(x))) replace(x, is.na(x), 0) else x))
    colnames(replaceNAdf_zero) <-  colnames(replaceNAdf)
    rownames(replaceNAdf_zero) <-  rownames(replaceNAdf)

    # If there is at least 1 value use the half minimum per feature
    replaceNAdf_Zero_MVI <- apply( replaceNAdf_zero, 2,  function(x) {x[is.na(x)] <-  min(x, na.rm = TRUE)/2
    return(x)
    }) %>% as.data.frame()
    rownames(replaceNAdf_Zero_MVI) <-  rownames(replaceNAdf)

    # add the samples in the original dataframe
    filtered_matrix_res <- rbind(NA_removed_matrix, replaceNAdf_Zero_MVI)
  }else{
    filtered_matrix_res <- NA_removed_matrix
  }

  ## ------------------ Return ------------------ ##
  invisible(return(filtered_matrix_res))
}


################################################################################################
### ### ### PreProcessing helper function: Total ion Count Normalization ### ### ###
################################################################################################

#' Total ion count normalisazion
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile").
#' @param TIC \emph{Optional: } If TRUE, Total Ion Count normalization is performed. If FALSE, only RLA QC plots are returned. \strong{Default = TRUE}
#'
#' @return List with two elements: DF (including output table) and Plot (including all plots generated)
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")
#' Res <- TICNorm(InputData=Intra[-c(49:58), -c(1:3)]%>% dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)),
#'                             SettingsFile_Sample=Intra[-c(49:58), c(1:3)],
#'                             SettingsInfo = c(Conditions = "Conditions"))
#'
#' @keywords total ion count normalisation
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyr pivot_longer
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot geom_boxplot geom_hline labs theme_classic theme_minimal theme annotation_custom aes_string
#' @importFrom logger log_info log_trace
#'
#' @noRd
#'
TICNorm <-function(InputData,
                   SettingsFile_Sample,
                   SettingsInfo,
                   TIC=TRUE){
  ## ------------------ Prepare the data ------------------- ##
  NA_removed_matrix <- InputData
  NA_removed_matrix[is.na(NA_removed_matrix)] <- 0#replace NA with 0

  ## ------------------ QC plot ------------------- ##
  ### Before TIC Normalization
  #### Log() transformation:
  log_NA_removed_matrix <- suppressWarnings(log(NA_removed_matrix) %>% t() %>% as.data.frame()) # log tranforms the data
  nan_count <- sum(is.nan(as.matrix(log_NA_removed_matrix)))# Count NaN values (produced by log(0))
  if (nan_count > 0) {# Issue a custom warning if NaNs are present
    warning(paste("For the RLA plot before/after TIC normalisation we have to perform log() transformation. This resulted in", nan_count, "NaN values due to 0s in the data."))
  }

  medians <- apply(log_NA_removed_matrix, 2, median) # get median
  RLA_data_raw <- log_NA_removed_matrix - medians   # Subtract the medians from each column
  RLA_data_long <- tidyr::pivot_longer(RLA_data_raw, cols = everything(), names_to = "Group")
  names(RLA_data_long)<- c("Samples", "Intensity")
  RLA_data_long <- as.data.frame(RLA_data_long)
  for (row in 1:nrow(RLA_data_long)){ # add conditions
    RLA_data_long[row, SettingsInfo[["Conditions"]]] <- SettingsFile_Sample[rownames(SettingsFile_Sample) %in%RLA_data_long[row,1],SettingsInfo[["Conditions"]]]
  }

  # Create the ggplot boxplot
  RLA_data_raw <- ggplot2::ggplot(RLA_data_long, ggplot2::aes_string(x = "Samples", y = "Intensity", color = SettingsInfo[["Conditions"]])) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "solid") +
    ggplot2::labs(title = "Before TIC Normalization")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggplot2::theme(legend.position = "none")

  #RLA_data_raw_Sized <- plotGrob_Processing(InputPlot = RLA_data_raw, PlotName= "Before TIC Normalization", PlotType= "RLA")

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
    log_Data_TIC  <- suppressWarnings(log(Data_TIC) %>% t() %>% as.data.frame()) # log tranforms the data
    medians <- apply(log_Data_TIC, 2, median)
    RLA_data_norm <- log_Data_TIC - medians   # Subtract the medians from each column
    RLA_data_long <- tidyr::pivot_longer(RLA_data_norm, cols = everything(), names_to = "Group")
    names(RLA_data_long)<- c("Samples", "Intensity")
    for (row in 1:nrow(RLA_data_long)){ # add conditions
      RLA_data_long[row, SettingsInfo[["Conditions"]]] <- SettingsFile_Sample[rownames(SettingsFile_Sample) %in%RLA_data_long[row,1],SettingsInfo[["Conditions"]]]
    }

    # Create the ggplot boxplot
    RLA_data_norm <- ggplot2::ggplot(RLA_data_long, ggplot2::aes_string(x = "Samples", y = "Intensity", color = SettingsInfo[["Conditions"]])) +
      ggplot2::geom_boxplot() +
      ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "solid") +
      ggplot2::labs(title = "After TIC Normalization")+
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      ggplot2::theme(legend.position = "none")

    #RLA_data_norm_Sized <- plotGrob_Processing(InputPlot = RLA_data_norm, PlotName= "After TIC Normalization", PlotType= "RLA")

    #Combine Plots
    dev.new()
    norm_plots <- suppressWarnings(gridExtra::grid.arrange(RLA_data_raw+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))+ ggplot2::theme(legend.position = "none"),
                                                           RLA_data_norm+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))+ ggplot2::theme(legend.position = "none"),
                                                           ncol = 2))
    dev.off()
    norm_plots <- ggplot2::ggplot() +ggplot2::theme_minimal()+ ggplot2::annotation_custom(norm_plots)


    ## ------------------ Return ------------------ ##
    Output_list <- list("DF" = list("Data_TIC"=as.data.frame(Data_TIC)),"Plot"=list( "norm_plots"=norm_plots, "RLA_AfterTICNorm"=RLA_data_norm,  "RLA_BeforeTICNorm" = RLA_data_raw ))
    invisible(return(Output_list))
  }else{
    ## ------------------ Return ------------------ ##
    Output_list <- list("Plot"=list("RLA_BeforeTICNorm" = RLA_data_raw))
    invisible(return(Output_list))
  }
}

################################################################################################
### ### ### PreProcessing helper function: CoRe nomalisation ### ### ###
################################################################################################

#' Consumption Release Normalisation
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", CoRe_norm_factor = "Columnname_Input_SettingsFile", CoRe_media = "Columnname_Input_SettingsFile"). Column CoRe_norm_factor is used for normalization and CoRe_media is used to specify the name of the media controls in the Conditions.
#'
#' @return List with two elements: DF (including output table) and Plot (including all plots generated)
#'
#' @examples
#' Media <- ToyData("CultureMedia_Raw")%>% subset(!Conditions=="Pool")%>% dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .))
#' Res <- CoReNorm(InputData= Media[, -c(1:3)],
#'                             SettingsFile_Sample= Media[, c(1:3)],
#'                             SettingsInfo = c(Conditions = "Conditions", CoRe_norm_factor = "GrowthFactor", CoRe_media = "blank"))
#'
#' @keywords Consumption Release Metaqbolomics, Normalisation, Exometabolomics
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot geom_histogram geom_vline geom_density labs theme_classic geom_violin geom_dotplot
#' @importFrom logger log_info log_trace
#' @importFrom dplyr case_when summarise_all mutate rowwise mutate_all select filter pull
#'
#'
#' @noRd
#'
#'
CoReNorm <-function(InputData,
                    SettingsFile_Sample,
                    SettingsInfo){
  ## ------------------ Prepare the data ------------------- ##
  Data_TIC <- InputData
  Data_TIC[is.na(Data_TIC)] <- 0

  ## ------------------ Perform QC ------------------- ##
  Conditions <- SettingsFile_Sample[[SettingsInfo[["Conditions"]]]]
  CoRe_medias <- Data_TIC[grep(SettingsInfo[["CoRe_media"]], Conditions),]

  if(dim(CoRe_medias)[1]==1){
    warning("Only 1 CoRe_media sample was found. Thus, the consistency of the CoRe_media samples cannot be checked. It is assumed that the CoRe_media samples are already summed.")
    CoRe_media_df <- CoRe_medias %>% t() %>% as.data.frame()
    colnames(CoRe_medias) <- "CoRe_mediaMeans"
  }else{
    ######################################################################################
    ## ------------------ QC Plots
    PlotList <- list()
    ##-- 1. PCA Media_control
    media_pca_data <- merge(x= SettingsFile_Sample %>% dplyr::select(SettingsInfo[["Conditions"]]), y= Data_TIC, by=0) %>%
      tibble::column_to_rownames("Row.names") %>%
      dplyr::mutate(Sample_type = dplyr::case_when(Conditions == SettingsInfo[["CoRe_media"]] ~ "CoRe_media",
                                     TRUE ~ "Sample"))

    media_pca_data[is.na( media_pca_data)] <- 0

    dev.new()
    pca_QC_media <-invisible(VizPCA(InputData=media_pca_data %>%dplyr::select(-SettingsInfo[["Conditions"]], -Sample_type),
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

    CutoffCV <- 30
    result_df <- result_df %>% t()%>%as.data.frame() %>% dplyr::rowwise() %>%
      dplyr::mutate(HighVar = CV > CutoffCV) %>% as.data.frame()
    rownames(result_df)<- colnames(CoRe_medias)

    # calculate the NAs
    NAvector <- apply(CoRe_medias, 2,  function(x) { (sum(is.na(x))/length(x))*100 })
    result_df$MissingValuePercentage <- NAvector

    cv_result_df <- result_df

    HighVar_metabs <- sum(result_df$HighVar == TRUE)
    if(HighVar_metabs>0){
      message(paste0(HighVar_metabs, " of variables have high variability (CV > 30) in the CoRe_media control samples. Consider checking the pooled samples to decide whether to remove these metabolites or not."))
    }

    #Make histogram of CVs
    HistCV <- invisible(ggplot2::ggplot(cv_result_df, aes(CV)) +
                          ggplot2::geom_histogram(aes(y=after_stat(density)), color="black", fill="white")+
                          ggplot2::geom_vline(aes(xintercept=CutoffCV),
                                     color="darkred", linetype="dashed", linewidth=1)+
                          ggplot2::geom_density(alpha=.2, fill="#FF6666") +
                          ggplot2::labs(title="CV for metabolites of control media samples (no cells)",x="Coefficient of variation (CV)", y = "Frequency")+
                          ggplot2::theme_classic())

    HistCV_Sized <- plotGrob_Processing(InputPlot = HistCV, PlotName= "CV for metabolites of control media samples (no cells)", PlotType= "Hist")

    PlotList[["Histogram_CoReMediaCV"]] <- HistCV_Sized

    #Make Violin of CVs
    Plot_cv_result_df <- cv_result_df %>%
      dplyr::mutate(HighVar = ifelse(HighVar == TRUE, "> CV 30", "< CV 30"))

    ViolinCV <- invisible(ggplot2::ggplot(Plot_cv_result_df, aes(y=CV, x=HighVar, label=row.names(cv_result_df)))+
                            ggplot2::geom_violin(alpha = 0.5 , fill="#FF6666")+
                            ggplot2::geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
                            ggrepel::geom_text_repel(aes(label = ifelse(Plot_cv_result_df$CV > CutoffCV,
                                                               as.character(row.names(Plot_cv_result_df)), '')),
                                            hjust = 0, vjust = 0,
                                            box.padding = 0.5, # space between text and point
                                            point.padding = 0.5, # space around points
                                            max.overlaps = Inf) + # allow for many labels
                            ggplot2::labs(title="CV for metabolites of control media samples (no cells)",x="Metabolites", y = "Coefficient of variation (CV)")+
                            ggplot2::theme_classic())

    ViolinCV_Sized <- plotGrob_Processing(InputPlot = ViolinCV, PlotName= "CV for metabolites of control media samples (no cells)", PlotType= "Violin")
    PlotList[["CoRe_Media_CV_Violin"]] <- ViolinCV_Sized

    ######################################################################################
    ## ------------------ Outlier testing
    if(dim(CoRe_medias)[1]>=3){
      Outlier_data <- CoRe_medias
      Outlier_data <- Outlier_data %>% dplyr::mutate_all(.funs = ~ FALSE)

      while(HighVar_metabs>0){
        #remove the furthest value from the mean
        if(HighVar_metabs>1){
          max_var_pos <-  CoRe_medias[,result_df$HighVar == TRUE]  %>%
            as.data.frame() %>%
            dplyr::mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
            dplyr::summarise_all(.funs = ~ which.max(abs(.)))
        }else{
          max_var_pos <-  CoRe_medias[,result_df$HighVar == TRUE]  %>%
            as.data.frame() %>%
            dplyr::mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
            dplyr::summarise_all(.funs = ~ which.max(abs(.)))
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

        result_df <- result_df %>% t()%>%as.data.frame() %>% dplyr::rowwise() %>%
          dplyr::mutate(HighVar = CV > CutoffCV) %>% as.data.frame()
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
        contingency_table[2, 2] <- dim(data_cont %>% dplyr::select(!all_of(sample)))[1]*dim(data_cont %>% dplyr::select(!all_of(sample)))[2] -sum( rowSums(data_cont) - current_sample)

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

      contingency_data_contframe <- contingency_data_contframe %>% dplyr::mutate(Total = rowSums(contingency_data_contframe))
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
      CoRe_medias <- CoRe_medias %>% dplyr::filter(!rownames(CoRe_medias) %in% different_samples)
    }
    CoRe_media_df <- as.data.frame(data.frame("CoRe_mediaMeans"=  colMeans(CoRe_medias, na.rm = TRUE)))
  }

  cv_result_df <- tibble::rownames_to_column(cv_result_df, "Metabolite")

  ######################################################################################
  ##------------------------ Substract mean (media control) from samples
  message("CoRe data are normalised by substracting mean (blank) from each sample and multiplying with the CoRe_norm_factor")
  ##-- Check CoRe_norm_factor
  if(("CoRe_norm_factor" %in% names(SettingsInfo))){
    CoRe_norm_factor <-   SettingsFile_Sample %>% dplyr::filter(!!as.name(SettingsInfo[["Conditions"]])!=SettingsInfo[["CoRe_media"]]) %>% dplyr::select(SettingsInfo[["CoRe_norm_factor"]]) %>%dplyr::pull()
    if(var(CoRe_norm_factor) ==  0){
      warning("The growth rate or growth factor for normalising the CoRe result, is the same for all samples")
    }
  }else{
    CoRe_norm_factor <- as.numeric(rep(1,dim(SettingsFile_Sample %>% dplyr::filter(!!as.name(SettingsInfo[["Conditions"]])!=SettingsInfo[["CoRe_media"]]))[1]))
  }

  # Remove CoRe_media samples from the data
  Data_TIC <- merge(SettingsFile_Sample, Data_TIC, by="row.names")%>%
    dplyr::filter(!!as.name(SettingsInfo[["Conditions"]])!=SettingsInfo[["CoRe_media"]])%>%
    tibble::column_to_rownames("Row.names")%>%
    dplyr::select(-1:-ncol(SettingsFile_Sample))

  Data_TIC_CoReNorm_Media <- as.data.frame(t( apply(t(Data_TIC),2, function(i) i-CoRe_media_df$CoRe_mediaMeans)))  #Subtract from each sample the CoRe_media mean
  Data_TIC_CoReNorm <- as.data.frame(apply(Data_TIC_CoReNorm_Media, 2, function(i) i*CoRe_norm_factor))

  #Remove CoRe_media samples from the data
  #Input_SettingsFile <- Input_SettingsFile[Input_SettingsFile$Conditions!="CoRe_media",]
  #Conditions <- Conditions[!Conditions=="CoRe_media"]

  ######################################################################################
  ##------------------------ Return Plots and Data
  if(dim(CoRe_medias)[1]>=3){
  DF_list <- list("CV_CoRe_blank" = cv_result_df, "Contigency_table_CoRe_blank" = contingency_data_contframe, "Core_Norm" = Data_TIC_CoReNorm)
  } else{
    DF_list <- list("CV_CoRe_blank" = cv_result_df, "Core_Norm" = Data_TIC_CoReNorm)
  }

  #Return
  Output_list <- list("DF"= DF_list,"Plot"=PlotList)
  invisible(return(Output_list))
}


################################################################################################
### ### ### PreProcessing helper function: Outlier detection ### ### ###
################################################################################################

#' OutlierDetection
#'
#' @param InputData DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected and consider converting any zeros to NA unless they are true zeros.
#' @param SettingsFile_Sample DF which contains information about the samples, which will be combined with the input data based on the unique sample identifiers used as rownames.
#' @param SettingsInfo  Named vector containing the information about the names of the experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile", Biological_Replicates="ColumnName_Plot_SettingsFile", CoRe_media = "Columnname_Input_SettingsFile"). Column "Conditions" with information about the sample conditions, Column "BiologicalReplicates" including numerical values and Column "Columnname_Input_SettingsFile" is used to specify the name of the media controls in the Conditions.
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed. If not normalised yet, provide information about control media sample names via SettingsInfo "CoRe_media" samples. \strong{Default = FALSE}
#' @param HotellinsConfidence \emph{Optional: } Confidence level for Hotellin's T2 test. \strong{Default = 0.99}
#'
#' @return List with two elements: : DF (including output tables) and Plot (including all plots generated)
#'
#' @examples
#' Intra <- ToyData("IntraCells_Raw")
#' Res <- OutlierDetection(InputData=Intra[-c(49:58), -c(1:3)]%>% dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)),
#'                                      SettingsFile_Sample=Intra[-c(49:58), c(1:3)],
#'                                      SettingsInfo = c(Conditions = "Conditions", Biological_Replicates = "Biological_Replicates"))
#'
#' @keywords Hotellins T2 outlier detection
#'
#' @importFrom dplyr relocate case_when mutate mutate_all select
#' @importFrom inflection uik
#' @importFrom factoextra fviz_screeplot
#' @importFrom ggplot2 ggplot theme_classic theme geom_vline annotate geom_line geom_point geom_hline scale_y_continuous ggtitle scale_linetype_discrete
#' @importFrom qcc mqcc
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom hash values keys
#' @importFrom logger log_info log_trace
#'
#' @noRd
#'
OutlierDetection <-function(InputData,
                            SettingsFile_Sample,
                            SettingsInfo,
                            CoRe=FALSE,
                            HotellinsConfidence=0.99){
  # Message:
  message("Identification of outlier samples is performed using Hotellin's T2 test to define sample outliers in a mathematical way (Confidence = 0.99 ~ p.val < 0.01) REF: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3), 360–378, doi:https://doi.org/10.1214/aoms/1177732979.")
  message(paste("HotellinsConfidence value selected:", HotellinsConfidence))

  # Load the data:
  data_norm <- InputData%>%
    dplyr::mutate_all(~ replace(., is.nan(.), 0))
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
      data_norm <- data_norm %>% dplyr::select(-all_of(metab))
    }

    ##---  PCA
    PCA.res <- prcomp(data_norm, center =  TRUE, scale. =  TRUE)
    outlier_PCA_data <- data_norm
    outlier_PCA_data$Conditions <- Conditions

    dev.new()
    pca_outlier <-invisible(VizPCA(InputData=data_norm,
                                               SettingsInfo= c(color=SettingsInfo[["Conditions"]]),
                                               SettingsFile_Sample= outlier_PCA_data,
                                               PlotName = paste("PCA outlier test filtering round ",loop),
                                               SaveAs_Plot =  NULL))

    if(loop==1){
      pca_outlierloop1 <- pca_outlier[["Plot_Sized"]][[1]]
    }
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
                                            linecolor = "black",linetype = 1) +
      ggplot2::theme_classic()+
      ggplot2::geom_vline(xintercept = npcs+0.5, linetype = 2, color = "red") +
      ggplot2::annotate("text", x = c(1:20),y = -0.8,label = screeplot_cumul,col = "black", size = 1.75)

    #screeplot_Sized <- plotGrob_Processing(InputPlot = screeplot, PlotName= paste("PCA Explained variance plot filtering round ",loop, sep = ""), PlotType= "Scree")

    if(loop==1){
      scree_outlierloop1 <-screeplot
    }
    dev.new()

    outlier_plot_list[[paste("ScreePlot_round",loop,sep="")]] <- screeplot # save plot
    dev.off()

    ##--- HotellingT2 test for outliers
    data_hot <- as.matrix(PCA.res$x[,1:npcs])
    hotelling_qcc <- qcc::mqcc(data_hot, type = "T2.single",labels = rownames(data_hot),confidence.level = HotellinsConfidence, title = paste("Outlier filtering via HotellingT2 test filtering round ",loop,", with ",HotellinsConfidence, "% Confidence",  sep = ""), plot = FALSE)
    HotellingT2plot_data <- as.data.frame(hotelling_qcc$statistics)
    HotellingT2plot_data <- tibble::rownames_to_column(HotellingT2plot_data, "Samples")
    colnames(HotellingT2plot_data) <- c("Samples", "Group summary statisctics")
    outlier <- HotellingT2plot_data %>% dplyr::filter(HotellingT2plot_data$`Group summary statisctics`>hotelling_qcc$limits[2])
    limits <- as.data.frame(hotelling_qcc$limits)
    legend <- colnames(HotellingT2plot_data[2])
    LegendTitle = "Limits"

    HotellingT2plot <- ggplot2::ggplot(HotellingT2plot_data, aes(x = Samples, y = `Group summary statisctics`, group = 1, fill = ))
    HotellingT2plot <- HotellingT2plot +
      ggplot2::geom_point(aes(x = Samples,y = `Group summary statisctics`), color = 'blue', size = 2) +
      ggplot2::geom_point(data = outlier, aes(x = Samples,y = `Group summary statisctics`), color = 'red',size = 3) +
      ggplot2::geom_line(linetype = 2)

    #draw the horizontal lines corresponding to the LCL,UCL
    HotellingT2plot <- HotellingT2plot +
      ggplot2::geom_hline(aes(yintercept = limits[,1]), color = "black", data = limits,  show.legend = F) +
      ggplot2::geom_hline(aes(yintercept = limits[,2], linetype = "UCL"), color = "red", data = limits, show.legend = T) +
      ggplot2::scale_y_continuous(breaks = sort(c(ggplot_build(HotellingT2plot)$layout$panel_ranges[[1]]$y.major_source, c(limits[,1],limits[,2]))))

    HotellingT2plot <- HotellingT2plot +
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
      ggplot2::ggtitle(paste("Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"))+
      ggplot2::scale_linetype_discrete(name = LegendTitle,)+
      ggplot2::theme(plot.title = element_text(size = 13))+#, face = "bold")) +
      ggplot2::theme(axis.text = element_text(size = 7))
    #HotellingT2plot_Sized <- plotGrob_Processing(InputPlot = HotellingT2plot, PlotName= paste("Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"), PlotType= "Hotellings")

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

  data_norm_filtered_full <- data_norm_filtered_full %>% dplyr::relocate(Outliers) #Put Outlier columns in the front
  data_norm_filtered_full <- merge(SettingsFile_Sample, data_norm_filtered_full,  by = 0) # add the design in the output df (merge by rownames/sample names)
  rownames(data_norm_filtered_full) <- data_norm_filtered_full$Row.names
  data_norm_filtered_full$Row.names <- c()

  ##-- 2.  Quality Control (QC) PCA
  MetaData_Sample <- data_norm_filtered_full %>%
    dplyr::mutate(Outliers = dplyr::case_when(Outliers == "no" ~ 'no',
                                Outliers == "Outlier_filtering_round_1" ~ ' Outlier_filtering_round = 1',
                                Outliers == "Outlier_filtering_round_2" ~ ' Outlier_filtering_round = 2',
                                Outliers == "Outlier_filtering_round_3" ~ ' Outlier_filtering_round = 3',
                                Outliers == "Outlier_filtering_round_4" ~ ' Outlier_filtering_round = 4',
                                TRUE ~ 'Outlier_filtering_round = or > 5'))
  MetaData_Sample$Outliers <- relevel(as.factor(MetaData_Sample$Outliers), ref="no")

  # 1. Shape Outliers
  if(length(sample_outliers)>0){
    dev.new()
    pca_QC <-invisible(VizPCA(InputData=as.data.frame(InputData)%>%dplyr::select(-zero_var_metab_export_df$Metabolite),
                                          SettingsInfo= c(color=SettingsInfo[["Conditions"]], shape = "Outliers"),
                                          SettingsFile_Sample= MetaData_Sample ,
                                          PlotName = "Quality Control PCA Condition clustering and outlier check",
                                          SaveAs_Plot =  NULL))
    dev.off()
    outlier_plot_list[["QC_PCA_and_Outliers"]] <- pca_QC[["Plot_Sized"]][[1]]
  }

  # 2. Shape Biological replicates
  if("Biological_Replicates" %in% names(SettingsInfo)){
    dev.new()
    pca_QC_repl <-invisible(VizPCA(InputData=as.data.frame(InputData)%>%dplyr::select(-zero_var_metab_export_df$Metabolite),
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
