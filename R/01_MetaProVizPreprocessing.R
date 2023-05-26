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
#' @param Input DF which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Experimental_design DF which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates"
#' @param Feature_Filtering \emph{Optional: }If set to "None" then no feature filtering is performed. If set to Standard then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Experimental_design including the individual conditions you want to apply the filtering to (Yang, J et al., 2015).\strong{Default = TRUE} \strong{Default = Modified}
#' @param Feature_Filt_Value \emph{Optional: } Percentage of feature filtering (Bijlsma S. et al., 2006).\strong{Default = 0.8}
#' @param TIC_Normalization \emph{Optional: } If TRUE, Total Ion Count normalization is performed. \strong{Default = TRUE}
#' @param HotellinsConfidence \emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. \strong{Default = 0.99}
#' @param ExportQCPlots \emph{Optional: } Select whether the quality control (QC) plots will be exported. \strong{Default = TRUE}
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed and the CoRe value will be calculated. Please consider providing a Normalisation factor column called "CoRe_norm_factor" in your "Experimental_design" DF, where the column "Conditions" matches. Th normalisation factor must be a numerical value obtained from growth rate that has been obtained from a growth curve or growth factor that was obtained by the ratio of cell count/protein quantification at the start point to cell count/protein quantification at the end point.. Additionally control media samples have to be available in the "Input" DF and defined as "blank" samples in the "Conditions" column in the "Experimental_design" DF, e.g. "blank_1", "blank_2". \strong{Default = FALSE}
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf, jpeg, tiff, bmp. \strong{Default = svg}
#'
#' @keywords 80% filtering rule, Missing Value Imputation, Total Ion Count normalization, PCA, HotellingT2, multivariate quality control charts,
#' @export


###################################################
### ### ### Metabolomics pre-processing ### ### ###
###################################################


Preprocessing <- function(Input_data,
                          Experimental_design,
                          Feature_Filtering = "Modified",
                          Feature_Filt_Value = 0.8,
                          TIC_Normalization = TRUE,
                          HotellinsConfidence = 0.99,
                          ExportQCPlots = TRUE,
                          CoRe = FALSE,
                          Save_as = "svg"
                          ){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "factoextra", "qcc", "ggplot2","hash", "inflection")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  suppressMessages(library(tidyverse))
  
  # library(tidyverse) # general scripting
  # library(factoextra) # visualize PCA
  # library(qcc) # for hotelling plots
  # library(ggplot2) # For visualization PCA
  # library(hash) # Dictionary in R for making column of outliers
  # library(inflection) # For finding inflection point/ Elbow knee /PCA component selection # https://cran.r-project.org/web/packages/inflection/inflection.pdf # https://deliverypdf.ssrn.com/delivery.php?ID = 454026098004123081018105104090015093000085002012023032095093077109069092095000114006057018122039107109012089110120018031068078025094036037013095100070100076109026029024044005068010070117123085122016083112098002109001027028000024115096122101001083084026&EXT = pdf&INDEX = TRUE # https://arxiv.org/abs/1206.5478

  ## ------------------ Run ------------------- ##
  
  #######################################################################################
  ### ### ### Check Input Information and add Experimental_design information ### ### ###

  #1.  Inpur data
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        Input_data <- Input_data
      )
    }
  }

  #2. Conditions
  if ( "Conditions" %in% colnames(Experimental_design)){   # Parse Condition and Replicate information
    Conditions <- Experimental_design$Conditions
  }else{
    Conditions <- NULL
  }

  #3. Core parameters
  if (CoRe ==  TRUE){   # parse CoRe normalisation factor
    message("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
    if ("CoRe_norm_factor" %in% colnames(Experimental_design)){
      CoRe_norm_factor <- Experimental_design$CoRe_norm_factor
    }else{
      warning("No growth rate or growth factor provided for normalising the CoRe result, hence CoRe_norm_factor set to 1 for each sample")
      CoRe_norm_factor <- as.numeric(rep(1,length(Experimental_design$Conditions)))
    }
    if (length(CoRe_norm_factor) !=  length(Experimental_design$Conditions)){ # Check is the length of normalization factor and conditions in 1 to 1
      stop("The CoRe_norm_factor length is different from the amount of samples. Please input a vector with a value for each sample. Blanks should take a value of 1.")
    }
    if(length(grep("blank", Conditions)) < 1){     # Check for blank samples
      stop("No blank samples were provided in the 'Conditions' in the Experimental design'. For a CoRe experiment control media samples have to also be measured and be added in the 'Conditions'
           column labeled as 'blank' (see @param section). Please make sure that you used the correct labelling or whether you need CoRE = TRUE for your analysis")
    }
  }
  
  #4. General parameters
  `%notin%` <- Negate(`%in%`) # Create a not in function
  Feature_Filtering_options <- c("Standard","Modified", "none")
  if(Feature_Filtering %notin% Feature_Filtering_options ){
    stop("Check input. The selected Feature_Filtering option is not valid. Please select one of the folowwing: ",paste(Feature_Filtering_options,collapse = ", "),"." )
  }
  if( is.numeric(Feature_Filt_Value) == FALSE |Feature_Filt_Value > 1 | Feature_Filt_Value < 0){
    stop("Check input. The selected Filtering value should be numeric and between 0 and 1.")
  }
  if(is.logical(TIC_Normalization) == FALSE){
    stop("Check input. The TIC_Normalization value should be either =TRUE if TIC normalization is to be performed or =FALSE if no data normalization is to be applied.")
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
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %notin% Save_as_options){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
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
  
  message("Feature filtering is performed to reduce missing values that can bias the analysis and cause methods to underperform, which leads to low precision in the statistical analysis. REF: Steuer et. al. (2007), Methods Mol Biol. 358:105-26., doi:10.1007/978-1-59745-244-1_7.")

  if (Feature_Filtering ==  "Modified"){
    message("Here we apply the modified 80%-filtering rule that takes the class information (Column `Conditions`) into account, which additionally reduces the effect of missing values. REF: Yang et. al., (2015), doi: 10.3389/fmolb.2015.00004)")
    message(paste("filtering value selected:", Feature_Filt_Value))

    Input_data <- as.data.frame(Input_data)
    unique_conditions <- unique(Conditions) # saves the different conditions

    if(is.null(unique(Conditions)) ==  TRUE){
      stop("Condition information is missing from the Experimental design")
    }
    if(length(unique(Conditions)) ==  1){
      stop("To perform the Modified feature filtering there have to be at least 2 different Conditions in the COndition column in the Experimental design. Consider using the Standard feature filtering option.")
    }

    miss <- c()
    message("***Performing modified feature filtering***")
    for (i in unique_conditions){
      split_Input <- split(Input_data, Conditions) # splits data frame into a list of dataframes by condition

      for (m in split_Input){ # Select metabolites to be filtered for different conditions
        for(i in 1:ncol(m)) {
          if(length(which(is.na(m[,i]))) > Feature_Filt_Value*nrow(m))  ## Check complete.case instead of is.na. it is faster and you dont have to use which
            miss <- append(miss,i)
        }
      }
    }
  
    if(length(miss) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      filtered_matrix <- Input_data
      feat_file_res <- "There where no metabolites exluded"
      write.table(feat_file_res,row.names =  FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    } else {
      message(paste( length(unique(miss)) ,"metabolites where removed:"))
      message(unique(colnames(Input_data)[miss]))
      filtered_matrix <- Input_data[,-miss]
      write.table(unique(colnames(Input_data)[miss]),row.names = FALSE, file =  paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    }
  }
  if (Feature_Filtering ==  "Standard"){
    message("Here we apply the so-called 80%-filtering rule is used, which removes metabolites with missing values in more than 80% of samples. REF: Smilde et. al. (2005), Anal. Chem. 77, 6729–6736., doi:10.1021/ac051080y")
    message(paste("filtering value selected:", Feature_Filt_Value))
    split_Input <- Input_data

    miss <- c()
    message("***Performing standard feature filtering***")
    for(i in 1:ncol(split_Input)) { # Select metabolites to be filtered for one condition
      if(length(which(is.na(split_Input[,i]))) > (Feature_Filt_Value)*nrow(split_Input))
        miss <- append(miss,i)
    }
    
    if(length(miss) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      filtered_matrix <- Input_data
      feat_file_res <- "There where no metabolites exluded"
      write.table(feat_file_res,row.names =  FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    } else {
      message(paste( length(unique(miss)) ,"metabolites where removed:"))
      message(unique(colnames(Input_data)[miss]))
      filtered_matrix <- Input_data[,-miss]
      write.table(unique(colnames(Input_data)[miss]),row.names =  FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
    }
  }
  if (Feature_Filtering ==  "None"){
    warning("No feature filtering is selected")
    filtered_matrix <- as.data.frame(Input_data)
  }
  
  filtered_matrix <- as.matrix(mutate_all(as.data.frame(filtered_matrix), function(x) as.numeric(as.character(x))))
  
  
  ################################################
  ### ### ###Zero variance metabolites ### ### ###
  
  filtered_matrix <- filtered_matrix %>% as.data.frame()
  zero_var_metabs <- filtered_matrix[, sapply(filtered_matrix, var,na.rm =  TRUE) ==  0] %>% colnames()
  zero_var_metabs <- paste(zero_var_metabs, collapse = ', ')
  if (length(zero_var_metabs) > 1){
    message(paste("There are metabolits with Zero variance. These are: ",zero_var_metabs ))
  }
  zero_var_removed <- filtered_matrix[, sapply(filtered_matrix, var,na.rm  =  TRUE) !=  0] # Remove zero variance metabolites
  filtered_matrix <- zero_var_removed
  
  
  ################################################
  ### ### ### Missing value Imputation ### ### ###

  message("Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0")
  NA_removed_matrix <- replace(filtered_matrix, filtered_matrix %in% NA, ((min(filtered_matrix, na.rm =  TRUE))/2))

  
  #######################################################
  ### ### ### Total Ion Current Normalization ### ### ###

  if (TIC_Normalization ==  TRUE){ # Ask if we want to normalize for total ion counts
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

  if (CoRe ==  TRUE){
    blankMeans <- colMeans( Data_TIC[grep("blank", Conditions),])
    blankSd <- as.data.frame( apply(Data_TIC[grep("blank", Conditions),], 2, sd) )

    blank_df <- as.data.frame(data.frame(blankMeans, blankSd))
    names(blank_df) <- c("blankMeans", "blankSd")
    blank_df$CV <- blank_df$blankSd / blank_df$blankMeans

    CV <- (sum(blank_df[blank_df$blankMeans !=  0,]$CV > 1)/length(blank_df[blank_df$blankMeans !=  0,]$CV))*100     # CV is the sd/mean and it sa measure of variability. above 1 is high and below is ok.
    message(paste0(CV, " of variables have very high variability in the blank samples"))
    Data_TIC_CoRe <- as.data.frame(t( apply(t(Data_TIC),2, function(i) i-blank_df$blankMeans)))  #Subtract from each sample the blank mean
    message("CoRe data are normalised using CoRe_norm_factor")
    Data_TIC <- apply(Data_TIC_CoRe, 2, function(i) i*CoRe_norm_factor)
    
    if (var(CoRe_norm_factor) ==  0){
      warning("The growth rate or growth factor for normalising the CoRe result, is the same for all samples") 
    }
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
  k =  1
  a =  1
  for (loop in 1:Outlier_filtering_loop){   # here we do 10 rounds of hotelling filtering

    metabolite_var <- as.data.frame( apply(data_norm, 2, var) %>% t()) # calculate each metabolites variance
    metabolite_zero_var_list <- list( colnames(metabolite_var)[which(metabolite_var[1,]==0)]) # takes the names of metablolites with zero variance and puts them in list
    
    if(sum(metabolite_var[1,]==0)==0){
      metabolite_zero_var_total_list[loop] <- 0
    } else if(sum(metabolite_var[1,]==0)>0){
    metabolite_zero_var_total_list[loop] <-metabolite_zero_var_list
    }
    
    for (metab in metabolite_zero_var_list){  # Remove the metabolites with zero variance from the data to do PCA
      data_norm <- data_norm %>% select(-metab)
    }
   
    ### ### PCA  ### ###
    PCA.res <- prcomp(data_norm, center =  TRUE, scale. =  TRUE) 
    outlier_PCA_data <- data_norm 
    outlier_PCA_data$Conditions <- Conditions
    pca.obj <- prcomp(data_norm, center =  TRUE, scale. =  TRUE)
    dtp <- data.frame('Conditions' = outlier_PCA_data$Conditions, pca.obj$x[,c(1,2)])
    pca_outlier <- ggplot(data = dtp) +
      geom_point(aes(x = PC1, y = PC2, colour = Conditions), size = 4, alpha = 0.8, show.legend = TRUE) +  
      ggtitle(paste("PCA outlier test filtering round ",loop))+
      theme_classic()+
      geom_hline(yintercept =  0,  color = "black", linewidth =  0.1)+
      geom_vline(xintercept =  0,  color = "black", linewidth = 0.1)+
      geom_text(aes(x = PC1, y = PC2, label = rownames(outlier_PCA_data)),hjust = 0.3, vjust = -0.5,size = 3,alpha = 0.6 )+
      scale_x_continuous(paste("PC1 ",summary(PCA.res)$importance[2,][[1]]*100,"%")) +
      scale_y_continuous(paste("PC2 ",summary(PCA.res)$importance[2,][[2]]*100,"%"))
    
    
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

    plot(screeplot)
    outlier_plot_list[[k]] <- recordPlot() # save plot
    dev.off()
    k = k+1

    ### ### HotellingT2 test for outliers ### ###
    data_hot <- as.matrix(PCA.res$x[,1:npcs])
    message("***Checking for outliers***")
    hotelling_qcc <- qcc::mqcc(data_hot, type = "T2.single",labels = rownames(data_hot),confidence.level = HotellinsConfidence, title = paste("Outlier filtering via HotellingT2 test filtering round ",loop,", with 99% Confidence",  sep = ""), plot = FALSE)
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

    plot(HotellingT2plot)
    outlier_plot_list[[k]] <- recordPlot()
    dev.off()
    k = k+1

    ### Save the outlier detection plots in the outlier detection folder
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "/PCA_OD_round_" ,a ,".", Save_as, sep = ""),
           plot = pca_outlier, width = 10,height = 8)
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "//Scree_plot_OD_round_" ,a ,".",Save_as, sep = ""),
           plot = screeplot, width = 10,height = 8)
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "/Hotelling_OD_round_" ,a ,".", Save_as, sep = ""),
           plot = HotellingT2plot, width = 10,height = 8)
    a = a+1

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
    message("There are sample outliers that were removed from the data") #This was a warning
    for (i in 1:length(sample_outliers)  ){
      message("Filtering round ",i ," Outlier Samples: ", paste( head(sample_outliers[[i]]) ," "))
    }
  }else{message("No sample outliers were found")}
  
  zero_var_metab_export_df <- data.frame(1,2)
  names(zero_var_metab_export_df) <- c("Filtering round","Metabolite")

  # Print zero variance metabolites
  count = 1
  for (i in 1:length(metabolite_zero_var_total_list)){  
    if (metabolite_zero_var_total_list[[i]] != 0){
      message("Filtering round ",i ,". Zero variance metabolites identified: ", paste( metabolite_zero_var_total_list[[i]] ," "))
      zero_var_metab_export_df[count,"Filtering round"] <- paste(i)
      zero_var_metab_export_df[count,"Metabolite"] <- paste(metabolite_zero_var_total_list[[i]])
      count = count +1
    }
  }
  write.table(zero_var_metab_export_df,row.names = FALSE, file =  paste(Results_folder_Preprocessing_folder,"/Zero_variance_metabolites",".csv",sep =  "")) # save zero var metabolite list
  
  
  #############################################
  ### ### ### Make Output Dataframe ### ### ###
  
  total_outliers <- hash::hash() # make a dictionary
  if(length(sample_outliers) > 0){ # Create columns with outliers to merge to output dataframe
    for (i in 1:length(sample_outliers)  ){
      total_outliers[[paste("Outlier_filtering_round_",i, sep = "")]] <- sample_outliers[i]
    }
  }
  
  data_norm_filtered_full <- as.data.frame(Data_TIC)
  
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
  data_norm_filtered_full <- merge(Experimental_design, data_norm_filtered_full,  by = 0) # add the design in the output df (merge by rownames/sample names)
  rownames(data_norm_filtered_full) <- data_norm_filtered_full$Row.names
  data_norm_filtered_full$Row.names <- c()

  
  ################################################
  ### ### ### Quality Control (QC) PCA ### ### ###

  QC_PCA_data <- Data_TIC %>% as.data.frame()
  QC_PCA_data$Conditions <- Experimental_design$Conditions
  if (is.null(Experimental_design$Biological_Replicates) != TRUE){
    QC_PCA_data$Biological_Replicates <-  as.character(Experimental_design$Biological_Replicates)
  }
  pca.obj <- prcomp(Data_TIC, center = TRUE, scale. = TRUE)
  
  dtp <- merge(data_norm_filtered_full %>% select(Conditions, Outliers), pca.obj$x[,1:2],by = 0)
  dtp <- dtp%>%  mutate(Outliers = case_when(Outliers == "no" ~ 'no',
                                  Outliers == "Outlier_filtering_round_1" ~ ' Outlier_filtering_round = 1',
                                  Outliers == "Outlier_filtering_round_2" ~ ' Outlier_filtering_round = 2',
                                  Outliers == "Outlier_filtering_round_3" ~ ' Outlier_filtering_round = 3',
                                  Outliers == "Outlier_filtering_round_4" ~ ' Outlier_filtering_round = 4',
                                  TRUE ~ 'Outlier_filtering_round = or > 5'))
  
  dtp$Outliers <- relevel( as.factor(dtp$Outliers), ref="no")

  ## PCA conditions and Outlier
  pca_QC <- ggplot(data = dtp) +
    geom_point(aes(x = PC1, y = PC2, color = Conditions, shape = Outliers), size = 4, alpha = 0.8) +  
    ggtitle("Quality COntrol PCA Condition clustering and Outlier check")+
    theme_classic()+
    geom_hline(yintercept = 0,  color = "black", linewidth = 0.1)+
    geom_vline(xintercept = 0,  color = "black", linewidth = 0.1)+
    geom_text(aes(x = PC1, y = PC2, label = rownames(QC_PCA_data)),hjust = 0.3, vjust = -0.5,size = 3,alpha = 0.6 )+
    scale_x_continuous(paste("PC1 ",summary(pca.obj)$importance[2,][[1]]*100,"%")) +
    scale_y_continuous(paste("PC2 ",summary(pca.obj)$importance[2,][[2]]*100,"%"))

  if (ExportQCPlots == TRUE){
   ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/PCA_Condition_Clustering.",Save_as), plot = pca_QC, width = 10,  height = 8)
  }
  
  if(is.null(Experimental_design$Biological_Replicates)!= TRUE){
    ### ### QC PCA color for replicates
    pca_QC_repl <- ggplot(data = dtp) +
      geom_point(aes(x = PC1, y = PC2, colour = Experimental_design$Conditions, shape = as.factor(Experimental_design$Biological_Replicates)), size = 4, alpha = 0.8) +  
      ggtitle("Quality Control PCA replicate spread check")+
      theme_classic()+
      geom_hline(yintercept = 0,  color = "black", linewidth = 0.1)+
      geom_vline(xintercept = 0,  color = "black", linewidth = 0.1)+
      geom_text(aes(x = PC1, y = PC2, label = rownames(QC_PCA_data)),hjust = 0.3, vjust = -0.5,size = 3,alpha = 0.6 )+
      scale_x_continuous(paste("PC1 ",summary(pca.obj)$importance[2,][[1]]*100,"%")) +
      scale_y_continuous(paste("PC2 ",summary(pca.obj)$importance[2,][[2]]*100,"%"))
    
    if (ExportQCPlots == TRUE){
      ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/PCA_replicate_distribution.",Save_as), plot = pca_QC_repl, width = 10,  height = 8)
    }
  }

  
  #########################################################
  ### ### ###  Make list with output dataframes ### ### ###

  output_list <- list()  #Here we make a list in which we will save the output
  preprocessing_output_list <- list(Experimental_design = Experimental_design, Raw_data = as.data.frame(Input_data), Processed_data = data_norm_filtered_full)

  ##Write to file
  preprocessing_output_list_out <- lapply(preprocessing_output_list, function(x) rownames_to_column(x, "Sample_ID")) #  # use this line to make a sample_ID column in each dataframe
  writexl::write_xlsx(preprocessing_output_list_out, paste(Results_folder_Preprocessing_folder, "/Preprocessing_output.xlsx", sep = ""))#,showNA = TRUE)

  # Return the result
  message("Done")
  return(preprocessing_output_list)
}






############################################################
### ### ### Merge analytical replicates function ### ### ###
############################################################


#' @param Input Dataframe which contains unique sample identifiers as row names the Experimental design and the metabolite numerical values in columns with metabolite identifiers as column names. Needs to have Conditions, Biological_Replicates and Analytical_Replicate columns
#' @keywords Analyrical Replicate Merge
#' @export Dataframe with merged Analytical Replicates based on Biological Replicates and Conditions


ReplicateSum <- function(Input_data){
  
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  ###############################################
  ### ### ### Check Input Information ### ### ###
  
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } 

  # Parse Condition and Replicate information
  if ( "Conditions" %in% colnames(Input_data)){
    Conditions <- Input_data$Conditions
  }else{
    Conditions <- NULL
  }
  if ( "Biological_Replicates" %in% colnames(Input_data)){
    Biological_Replicates <- Input_data$Biological_Replicates
  }else{
    Biological_Replicates <- NULL
  }
  if ( "Analytical_Replicates" %in% colnames(Input_data)){
    Analytical_Replicates <- Input_data$Analytical_Replicates
  }else{
    Analytical_Replicates <- NULL
  }
  
  ############################################
  ### ### ### Create Output folder ### ### ###
  
  # This searches for a folder called "Results" within the current working directory and if its not found it creates one
  Results_folder = paste(getwd(), "/MetaProViz_Results_",Sys.Date(),  sep =  "")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
  Results_folder_Preprocessing_folder = paste(Results_folder, "/MetaProViz_Preprocessing", sep = "")
  if (!dir.exists(Results_folder_Preprocessing_folder)) {dir.create(Results_folder_Preprocessing_folder)}  # check and create folder
  
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
  

  
  Input_data_numeric_summed <- merge(nReplicates,Input_data_numeric_summed, by = c("Conditions","Biological_Replicates"))
  
  # Export result
  writexl::write_xlsx(Input_data_numeric_summed, paste(Results_folder_Preprocessing_folder, "/ReplicateSum_output.xlsx", sep = ""))#,showNA = TRUE)
  
  # Return the result
  message("Done")
  return(Input_data_numeric_summed)
}



