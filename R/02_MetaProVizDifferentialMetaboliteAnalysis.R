#' Metabolomics pre-processing pipeline
#' @author Prymidis Dimitrios, Schmidt Christina
#' Date: "2022-10-28"
#'
#' This script allows you to perform differential metabolite analysis
#'
#' @param Input Data matrix which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Experimental_design Data matrix which contains information about the samples. Conditions, replicates etc. (i.e. "1 and 2" or "N and T" or
#' "Normal and Tumor").
#' @param padj_test String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing
#' @param Output_Name String which contains the name of the output file of the DMA
#' @param Input_Pathways A dataframe which contains the pathway information for each metabolite
#' #'
#' @keywords Differential Metabolite Analysis,
#' @export
#'
#'
#'
# Load libraries
#library(tidyverse) # general scripting
#library(gtools) # for calculating Log2FC
#library(EnhancedVolcano) # for Volcano plots


########################################################
### ### ### Differential Metabolite Analysis ### ### ###
########################################################

MetaProVizDMA <-function(Input_data =Data,
                         Experimental_design= Design,
                         padj_test="fdr", # search:?p.adjust for more # Methods:"BH" for Benjamini & Hochberg, "fdr", "bonferroni", "holm",
                         #Output_Name="DMA_Output_Condition1-versus-Condition2", not used I generate this name in the functions from the Conditions
                         Condition1 = Condition1,
                         Condition2 = Condition2,
                         STAT_pval ="t-test", # or wilcox-test
                         padjVALs = "0.05", #    padjVALs =   c("0.01","0.05","0.1")
                         log2FCs = "0.5", #    log2FCs = c("0.5","1","2")
                         Input_Pathways = NULL,
                         OutputName= '',
                         CoRe=FALSE,
                         Save_as = svg  #   Save_as = "svg"
){
  

  
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "gtools", "EnhancedVolcano")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  
  ####################################################
  ### ### ### Create Results output folder ### ### ###
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  MetaProViz_results_folder <- file.path(WorkD, name) # Make Results folder
  if (!dir.exists(MetaProViz_results_folder)) {dir.create(MetaProViz_results_folder)}
  MetaProViz_results_folder_DMA_folder <- file.path(MetaProViz_results_folder,"MetaProVizDMA") # Make DMA results folder
  if (!dir.exists(MetaProViz_results_folder_DMA_folder)) {dir.create(MetaProViz_results_folder_DMA_folder)} 
  MetaProViz_results_folder_Conditions <- file.path(MetaProViz_results_folder_DMA_folder,paste0(Condition1,"_vs_",Condition2)) # Make comparison folder
  if (!dir.exists(MetaProViz_results_folder_Conditions)) {dir.create(MetaProViz_results_folder_Conditions)} 

  
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  Save_as= deparse(substitute(Save_as))
  
  
  if(length(Condition1)>1 |length(Condition2)>1 ){
    Pooled=TRUE
  }else{
    Pooled=FALSE
  }
  
  # This is when we use the data_processed_summed data
  if ('Conditions' %in% colnames(Input_data)){
    # Cond <- Input_data$Conditions
    Input_data<- select_if(Input_data, is.numeric)%>%
      select(-Analytical_Replicates, -Biological_Replicates)
    
  }else{ # This is when we use the data_processed
    # Take only the colums that are numeric
    Input_data<- select_if(Input_data, is.numeric)
    # Parse from the Experimental_design only the sample kept
    Experimental_design <- Experimental_design %>%
      filter(rownames(Experimental_design) %in% rownames(Input_data))
  }
  
  # select samples according to the input conditions
  C1 <- Experimental_design %>% filter(Experimental_design$Conditions %in% Condition1) %>% rownames()
  C2 <- Experimental_design %>% filter(Experimental_design$Conditions %in% Condition2) %>% rownames()
  Data_DMA_Cond1 <- Input_data %>% filter(rownames(Input_data) %in% C1)
  Data_DMA_Cond2 <- Input_data %>% filter( rownames(Input_data) %in% C2)
  
  
  Log2FC_Condition1 <- Input_data %>% filter(rownames(Input_data) %in% C1) %>%  summarise_all("mean")
  Log2FC_Condition2 <- Input_data %>% filter(rownames(Input_data) %in% C2) %>%  summarise_all("mean")
  #add +1 count to the metabolites with 0 because otherwise the Log2FC wives NA
  Log2FC_Condition1[,which(Log2FC_Condition1[1,]==0)] <- Log2FC_Condition1[,which(Log2FC_Condition1[1,]==0)]+1
  Log2FC_Condition2[,which(Log2FC_Condition2[1,]==0)] <- Log2FC_Condition2[,which(Log2FC_Condition2[1,]==0)]+1
  
  message("*** Performing Differential Metabolite Analysis ***")
  
  ########################################
  ### ### ### Calculate Log2FC ### ### ###
  
  if(CoRe==TRUE){
    #caluclate the log2 by using the distance between the means
    #test code delete it
    temp1 <- Log2FC_Condition1
    temp2 <- Log2FC_Condition2
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
  }
    
  FC_C1vC2 <- mapply(gtools::foldchange,Log2FC_Condition1,Log2FC_Condition2)
  Log2FC_C1vC2 <- as.data.frame(gtools::foldchange2logratio(FC_C1vC2, base=2))
  Log2FC_C1vC2 <- cbind(rownames(Log2FC_C1vC2), data.frame(Log2FC_C1vC2, row.names=NULL))
  names(Log2FC_C1vC2)[names(Log2FC_C1vC2) == "rownames(Log2FC_C1vC2)"] <- "Metabolite"
  names(Log2FC_C1vC2)[names(Log2FC_C1vC2) == "gtools..foldchange2logratio.FC_C1vC2..base...2."] <- "Log2FC"

  
  # remove Nas and Infs
  # Have this or the +1 in the counts some lines above.
  #Log2FC_C1vC2 <-  Log2FC_C1vC2[complete.cases(Log2FC_C1vC2),]
  #Log2FC_C1vC2 <- Log2FC_C1vC2[!is.infinite(Log2FC_C1vC2$Log2FC),]
  
  ############################################
  ### ### ### Check data normality ### ### ###
  
  # Before Hypothesis testing, we have to decide whether to use a parametric or a non parametric test.
  # To test the data normality we use the Shapiro test. If most of the metabolites follow a normal distribution we will use parametric test. IF not then we
  # will use a non parametric test
  
  # Select data
  shaptest <-  Input_data %>%
    filter(Experimental_design$Conditions %in% Condition1 | Experimental_design$Conditions %in% Condition2)
  
  
  #  we have to remove features with zero variance if there are any.
  temp<- as.vector(sapply(shaptest, function(x) var(x)) == 0)
  shaptest <- shaptest[,!temp]
  
  shapirotres <- as.data.frame(sapply(shaptest, function(x) shapiro.test(x))) # do the test for each metabolite
  shapirotres <- as.data.frame(t(shapirotres))
  
  #  The data are normal if the p-value is above 0.05.
  if (sum(shapirotres$p.value > 0.05) > sum(shapirotres$p.value < 0.05)){
    message(paste0( "Most metabolites follow normal distribution" ," (" ,format(round(sum(shapirotres$p.value > 0.05)/dim(shaptest)[2],2), nsmall = 2), "%).", "Thus, a t-test (parametric Hypothesis testing) should be used."  ))
    #STAT_pval =t.test
    if(STAT_pval =="wilcox-test"){
      warning("The data are normally distributed. However, 'wilcox.test' is selected for Hypothesis testing. Please change it to 't-test'.")
    }
    
  }else{
    message(paste0( "Most metabolites are not normally distributed" ," (" ,format(round(sum(shapirotres$p.value < 0.05)/dim(shaptest)[2],2), nsmall = 2), "%).", "Thus, wilcox.test (non parametric Hypothesis testing) should be used."))
    #STAT_pval =wilcox.test
    
    if(STAT_pval =="t-test"){
      warning("The data are not normally distributed. However, 't.test' is selected for Hypothesis testing. Please change it to 'wilcox-test'.")
    }
  }
  
  # check number of samples for wilcoxons test
  if(STAT_pval =="wilcox-test"){
    if(dim(Data_DMA_Cond1)[1]<5 | dim(Data_DMA_Cond2)[1]<5){
      stop(paste0( "Very small number of samples for using wilcox-test. Consider using, we should use a t-test."  ))
    }
  }
    
  
  # select the test for hypothesis testing
  if(STAT_pval=="t-test"){
    STAT_pval= t.test
  }else if(STAT_pval=="wilcox-test") {
    STAT_pval=wilcox.test
  }else{
    stop("Please select an apropriate hypothesis testing option.")
  }
  
  
  ##################################################
  ### ### ### Perform Hypothesis testing ### ### ###
  
  T_C1vC2 <-mapply(STAT_pval, x= as.data.frame(Data_DMA_Cond2), y = as.data.frame(Data_DMA_Cond1), SIMPLIFY = F)
  
  VecPVAL_C1vC2 <- c()
  for(i in 1:length(T_C1vC2)){
    p_value <- unlist(T_C1vC2[[i]][3])
    VecPVAL_C1vC2[i] <- p_value
  }
  Metabolite <- colnames(Data_DMA_Cond2)
  PVal_C1vC2 <- data.frame(Metabolite, VecPVAL_C1vC2)
  #p-adjusted
  #library(stats)
  VecPADJ_C1vC2 <- p.adjust((PVal_C1vC2[,2]),method = padj_test, n = length((PVal_C1vC2[,2])))
  PADJ_C1vC2 <- data.frame(Metabolite, VecPADJ_C1vC2)
  STAT_C1vC2 <- merge(PVal_C1vC2,PADJ_C1vC2, by="Metabolite")
  STAT_C1vC2 <- merge(Log2FC_C1vC2,STAT_C1vC2, by="Metabolite")
  names(STAT_C1vC2)[names(STAT_C1vC2) == "VecPVAL_C1vC2"] <- "p.val"
  names(STAT_C1vC2)[names(STAT_C1vC2) == "VecPADJ_C1vC2"] <- "p.adj"
  
  # add consumption release info to the result
  if(CoRe==TRUE){
    CoRe_info <- t(CoRe_info) %>% as.data.frame()
    CoRe_info <- rownames_to_column(CoRe_info, "Metabolite")
    #CoRe_info <- CoRe_info[,c(1,4)]
    names(CoRe_info)[2] <- paste("Mean", "CoRe", Condition1, sep="_")
    names(CoRe_info)[3] <- paste("Mean", "CoRe", Condition2, sep="_")
    names(CoRe_info)[4] <- "CoRe"
    STAT_C1vC2 <- merge(STAT_C1vC2,CoRe_info,by= "Metabolite")
  }
  
  #########################################################################
  ### ### ###  Check and Add pathway information to DMA results ### ### ###
  
  if(is.null(Input_Pathways)!=TRUE){
    if('Metabolite' %in% colnames(Input_Pathways) & 'Pathway' %in% colnames(Input_Pathways)){
      STAT_C1vC2$Metabolite %in% Input_Pathways$Metabolite
      STAT_C1vC2<- merge(STAT_C1vC2,Input_Pathways,by="Metabolite", all.x=T)
      STAT_C1vC2$Pathway[  is.na(STAT_C1vC2$Pathway)] <- "unknown"    # Merge the pathways to DMA result. All non matching metabolites get NA that are changed into "unknown".
    }else{
      stop("The pathway data must have 2 columns named: Metabolite and Pathway")
    }
  }else{
    message("No pathway information was added")
  }

  DMA_Output <- STAT_C1vC2
  
  padjVALs = as.numeric(unlist(strsplit(padjVALs, ",")))
  log2FCs = as.numeric(unlist(strsplit(log2FCs, ",")))
  DMA_Output_out_list <- list(all = DMA_Output)
  for (log2.FC in log2FCs){
    for (padjVAL in padjVALs){
      # Subset results based on different significance thresholds adn save them into different sheets in output excel file
      DMA_Output_padj_Log2FC <-  DMA_Output %>% filter(p.adj < padjVAL & abs(Log2FC) > log2.FC)# %>% order(DMA_Heart$Log2FC)
      #DMA_Output_padj_Log2FC <- DMA_Output[which(DMA_Output$p.adj < padjVAL & abs(DMA_Output$Log2FC) > log2.FC),]
      DMA_Output_padj_Log2FC <- DMA_Output_padj_Log2FC[order(DMA_Output_padj_Log2FC$Log2FC), ]
      
      # Save data
      if (nrow(DMA_Output_padj_Log2FC) > 0) {
        
        sheet = paste0("abslog2fold>", log2.FC, " padj<", padjVAL)
        DMA_Output_out_list[[sheet]] <- DMA_Output_padj_Log2FC
        
      }
    }
  }
  
  # Save the DMA results table
  xlsDMA <- file.path(MetaProViz_results_folder_Conditions,paste0("DMA_Output_",Condition1,"_vs_",Condition2,"_", OutputName, ".xlsx"))
  
  if (Pooled == TRUE){
    writexl::write_xlsx(DMA_Output_out_list,xlsDMA, col_names = TRUE)
    
  } else{
    writexl::write_xlsx(DMA_Output_out_list,xlsDMA, col_names = TRUE)
  }

# Make a simple Volcano plot
  for (padjVAL in padjVALs){
    for (log2.FC in log2FCs){
  VolcanoPlot<- EnhancedVolcano::EnhancedVolcano(DMA_Output,
                                 lab = DMA_Output$Metabolite,#Metabolite name
                                 x = "Log2FC",#Log2FC
                                 y = "p.adj",#p-value or q-value
                                 xlab = bquote(~Log[2]~ FC),
                                 ylab = bquote(~-Log[10]~p.adj),#(~-Log[10]~adjusted~italic(P))
                                 pCutoff = padjVAL,
                                 FCcutoff = log2.FC,#Cut off Log2FC, automatically 2
                                 pointSize = 3,
                                 labSize = 2,
                                 titleLabSize = 16,
                                 # colCustom = c("black", "grey", "grey", "red"),
                                 colAlpha = 0.7,
                                 title= paste0(Condition1,"-vs-",Condition2),
                                 subtitle = bquote(italic("Differential metabolite analysis")),
                                 caption = paste0("total = ", nrow(DMA_Output), " Metabolites"),
                                 xlim =  c(min(DMA_Output$Log2FC[is.finite(DMA_Output$Log2FC )])-0.2,max(DMA_Output$Log2FC[is.finite(DMA_Output$Log2FC )])+0.2  ),
                                 ylim = c(0,(ceiling(-log10(Reduce(min,DMA_Output$p.adj))))),
                                 cutoffLineType = "dashed",
                                 cutoffLineCol = "black",
                                 cutoffLineWidth = 1,
                                 legendLabels=c('No changes',paste(log2.FC,"< |Log2FC|"),paste("p.adj <",padjVAL) , paste('p.adj<',padjVAL,' &',log2.FC,"< |Log2FC|")),
                                 legendPosition = 'right',
                                 legendLabSize = 8,
                                 legendIconSize =4,
                                 gridlines.major = FALSE,
                                 gridlines.minor = FALSE,
                                 drawConnectors = FALSE)
  OutputPlotName = paste0(OutputName,"_padj_",padjVAL,"log2FC_",log2.FC)
  
  
  volcanoDMA <- file.path(MetaProViz_results_folder_Conditions,paste0( "Volcano_Plot_",Condition1,"-versus-",Condition2,"_", OutputPlotName,".",Save_as))
  if(Pooled==TRUE){
    ggsave(volcanoDMA,plot=VolcanoPlot, width=12, height=9)
  }else{
    ggsave(volcanoDMA,plot=VolcanoPlot, width=12, height=9)
  }
  
  
    }
  }
  return(DMA_Output)
}



########################################
### ### ### Run the Analysis ### ### ###

#DMA_output <- DMA(Input_data =preprocessing_output[["data_processed"]] ,
# Experimental_design= preprocessing_output[["Experimental_design"]],
#padj_test="BH",
# Condition1 = "Control",
# Condition2 = "Rot",
# Input_Pathways = pathways)





#' Notes/ Things to add
#' Check for heteroscedasticity alongside normality?? in order to apply a transformation, log if n is small and vst if n is large. But how smass is small? For this
#' the problem will be solved as we will do thr apgelm normalization is the preprocessing step. So it should be solved ina previous step.

# Let the user decide in which file format the df is safed (xlsx/csv/tsv/...)
# Let the user decide which test to use when normally distributed (default=t.test) or when not nromally disctributed (default=wilcox.test)
# Question to check: If the amount of normally distributed to not-normally distributed metabolites is 51%:49% is it ok to perfom a t.test? (is there a threshold for those decisions?)



