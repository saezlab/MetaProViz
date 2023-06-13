#' Metabolomics pre-processing pipeline
#' @author Prymidis Dimitrios, Schmidt Christina
#' Date: "2022-10-28"
#'
#' This script allows you to perform metabolite clustering analysis
#'
#' @param Input_data Dataframe which contains metabolites in rows and Log fold changes, pvalues and padjusted values in columns.
#' @param Input_data2 Dataframe same as Input_data for another comparison.
#' @param Output_Name String which is added to the output files of the Metabolic Clusters.
#' @param Condition1 String which contains the name of the first condition.
#' @param Condition2 String which contains the name of the second condition.
#' @param pCutoff Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param FCcutoff Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param test String which selects pvalue or padj for significance.  \strong{Default = padj}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used. \strong{Default = FALSE}
#' @param plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg or pdf. \strong{Default = svg}
#'
#'
#' @keywords Metabolic Clusters,
#' @export


##################################################
### ### ### Metabolite Clustering Analysis ### ### ###
##################################################

MCA <- function(Input_data,
                Input_data2,
                Output_Name = "",
                Condition1,
                Condition2,
                pCutoff= 0.05 ,
                FCcutoff=0.5,
                test = "p.adj",
                plot = TRUE,
                CoRe=FALSE,
                Save_as = "svg"
               ){
  
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "alluvial")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))
  
  
  ## ------------ Check Input files ----------- ##
  for(data in list(Input_data, Input_data2)){
    if(any(duplicated(row.names(data)))==TRUE){
      stop("Duplicated row.names of Input_data, whilst row.names must be unique")
    } 
  }
  if( is.numeric(pCutoff)== FALSE |pCutoff > 1 | pCutoff < 0){
    stop("Check input. The selected pCutoff value should be numeric and between 0 and 1.")
  }
  if( is.numeric(FCcutoff)== FALSE  | FCcutoff < 0){
    stop("Check input. The selected pCutoff value should be numeric and between 0 and +oo.")
  }
  if(test != "p.val" & test != "p.adj" ){
    stop("Check input. The selected test option for assessing significance is not valid. Please select one of the following: p.adj, p.val.")
  }
  for(Input_data in list(Input_data, Input_data2)){
    if(test %in% colnames(Input_data) == FALSE){
      stop("Check Input data. There is no column ", test, " for assessing significance.")
    } 
  }
  if(is.logical(plot) == FALSE){
    stop("Check input. The plot value should be either =TRUE if an Alluvial plot presenting the Metabolic Cluster Analysisresults is to be exported or =FALSE if not.")
  }
  if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for analysis of Consuption/Release experiment or =FALSE if not.")
  }
  Save_as_options <- c("svg","pdf")
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the following: ",paste(Save_as_options,collapse = ", "),"." )
  }
  
  
  ####################################################
  ### ### ### Create Results output folder ### ### ###
  
  name <- paste0("MetaProViz_Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name) # Make Results folder
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  Results_folder_MCA_folder = file.path(Results_folder, "MCA") # select name of result directory
  if (!dir.exists(Results_folder_MCA_folder)) {dir.create(Results_folder_MCA_folder)}  # check and create folder
  
  
  if (CoRe == TRUE){
    
    Name = paste0("MCA_Output_",gsub(" ", "_",Condition1),"_with_",gsub(" ", "_",Condition2), sep = "")
    
    C1 <- Input_data
    C1 <- na.omit(C1)
    C2 <- Input_data2
    C2 <- na.omit(C2)
    
    #C1 <- C1 %>% drop_na()
    #C2 <- C2 %>% drop_na()
    
    # Make intracellular data C1 and CoRe C2
    if ("CoRe" %in% colnames(C1)){
      DMA_Intra <- C2
      DMA_CoRe <- C1
    }else if ("CoRe" %in% colnames(C2)){
      DMA_Intra <- C1
      DMA_CoRe <- C2
    }else{
      stop("No CoRe column was found")
    }
    
    # Remove the metabolite Means column from CoRe
    # grep("Mean",colnames(DMA_CoRe))
    # DMA_CoRe <- DMA_CoRe[,-grep("Mean",colnames(DMA_CoRe))]
    
    # 0. Overall Regulation: label the selection of metabolites that change or where at least we have a change in one of the two conditions
    
    # Merge the Intracellular and CoRe results
    DMA_Intra <- DMA_Intra%>%
      rename("Log2FC_Intra"="Log2FC",
             "p.val_Intra"="p.val",
             "p.adj_Intra"="p.adj")
    
    DMA_CoRe <- DMA_CoRe%>%
      rename("Log2FC_CoRe"="Log2FC",
             "p.val_CoRe"="p.val",
             "p.adj_CoRe"="p.adj")
    
    
    DMA <- merge( DMA_Intra, DMA_CoRe, by="Metabolite")
    
    
    # Define the clusters based on p.val or p.adj
    if(test== "p.val"){
      DMA <- DMA%>%
        mutate(`Intracellular Change` = case_when(Log2FC_Intra >= 0.5 & p.val_Intra< 0.05 ~ 'UP',
                                                  Log2FC_Intra <= -0.5 & p.val_Intra< 0.05 ~ 'DOWN',
                                                  TRUE ~ 'No Change'))%>%
        mutate(`CoRe Change` = case_when(Log2FC_CoRe >= 0.5 & p.val_CoRe < 0.05 & CoRe =="Released" ~ 'Release UP',
                                         Log2FC_CoRe <= -0.5 & p.val_CoRe < 0.05 & CoRe =="Released" ~ 'Release DOWN',
                                         Log2FC_CoRe >= 0.5 & p.val_CoRe < 0.05 & CoRe =="Consumed"~ 'Consume UP',
                                         Log2FC_CoRe <= -0.5 & p.val_CoRe < 0.05 & CoRe =="Consumed"~ 'Consume DOWN',
                                         Log2FC_CoRe >= 0.5 & p.val_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed UP',
                                         Log2FC_CoRe >= 0.5 & p.val_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed UP',
                                         Log2FC_CoRe <= -0.5 & p.val_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed DOWN',
                                         Log2FC_CoRe <= -0.5 & p.val_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed DOWN',
                                         TRUE ~ 'No Change'))
      
      # match("Released" ,str_split(DMA$CoRe[a], " ", simplify = TRUE),nomatch= FALSE) < 
      #   match("Consumed" ,str_split(DMA$CoRe[a], " ", simplify = TRUE),nomatch= FALSE) &
      #   match("Released" ,str_split(DMA$CoRe[a], " ", simplify = TRUE),nomatch= FALSE) !=0 &
      #   match("Consumed" ,str_split(DMA$CoRe[a], " ", simplify = TRUE),nomatch= FALSE) !=0
      
    }else{ # else if(test== "p.adj")
      DMA <- DMA%>%
        mutate(`Intracellular Change` = case_when(Log2FC_Intra >= 0.5 & p.adj_Intra< 0.05 ~ 'UP',
                                                  Log2FC_Intra <= -0.5 & p.adj_Intra< 0.05 ~ 'DOWN',
                                                  TRUE ~ 'No Change'))%>%
        mutate(`CoRe Change` = case_when(Log2FC_CoRe >= 0.5 & p.adj_CoRe < 0.05 & CoRe =="Released" ~ 'Release UP',
                                         Log2FC_CoRe <= -0.5 & p.adj_CoRe < 0.05 & CoRe =="Released" ~ 'Release DOWN',
                                         Log2FC_CoRe >= 0.5 & p.adj_CoRe < 0.05 & CoRe =="Consumed"~ 'Consume UP',
                                         Log2FC_CoRe <= -0.5 & p.adj_CoRe < 0.05 & CoRe =="Consumed"~ 'Consume DOWN',
                                         Log2FC_CoRe >= 0.5 & p.adj_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed UP',
                                         Log2FC_CoRe >= 0.5 & p.adj_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed UP',
                                         Log2FC_CoRe <= -0.5 & p.adj_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed DOWN',
                                         Log2FC_CoRe <= -0.5 & p.adj_CoRe < 0.05 & all(c("Released", "Consumed") %in% str_split(DMA$CoRe," ",simplify = TRUE)) ~ 'Released/Consumed DOWN',
                                         TRUE ~ 'No Change'))
    }
    
    
    DMA <- DMA%>%
      mutate(`Metabolite Cluster` = case_when(`Intracellular Change` == "UP" & `CoRe Change`== "Release UP" ~ '1',
                                              `Intracellular Change` == "UP" & `CoRe Change`== "Release DOWN" ~ '2',
                                              `Intracellular Change` == "UP" & `CoRe Change`== "Consume UP" ~ '3',
                                              `Intracellular Change` == "UP" & `CoRe Change`== "Consume DOWN" ~ '4',
                                              `Intracellular Change` == "UP" & `CoRe Change`== "Released/Consumed UP" ~ '5',
                                              `Intracellular Change` == "UP" & `CoRe Change`== "Released/Consumed DOWN" ~ '6',
                                              `Intracellular Change` == "UP" & `CoRe Change`== "No Change" ~ '7',
                                              `Intracellular Change` == "No Change" & `CoRe Change`== "Release UP" ~ '8',
                                              `Intracellular Change` == "No Change" & `CoRe Change`== "Release DOWN" ~ '9',
                                              `Intracellular Change` == "No Change" & `CoRe Change`== "Consume UP" ~ '10',
                                              `Intracellular Change` == "No Change" & `CoRe Change`== "Consume DOWN" ~ '11',
                                              `Intracellular Change` == "No Change" & `CoRe Change`== "Released/Consumed UP" ~ '12',
                                              `Intracellular Change` == "No Change" & `CoRe Change`== "Released/Consumed DOWN" ~ '13',
                                              `Intracellular Change` == "No Change" & `CoRe Change`== "No Change" ~ '14',
                                              `Intracellular Change` == "DOWN" & `CoRe Change`== "Release UP" ~ '15',
                                              `Intracellular Change` == "DOWN" & `CoRe Change`== "Release DOWN" ~ '16',
                                              `Intracellular Change` == "DOWN" & `CoRe Change`== "Consume UP" ~ '17',
                                              `Intracellular Change` == "DOWN" & `CoRe Change`== "Consume DOWN" ~ '18',
                                              `Intracellular Change` == "DOWN" & `CoRe Change`== "Released/Consumed UP" ~ '19',
                                              `Intracellular Change` == "DOWN" & `CoRe Change`== "Released/Consumed DOWN" ~ '20',
                                              `Intracellular Change` == "DOWN" & `CoRe Change`== "No Change" ~ '21',
                                              Log2FC_Intra <= -0.5 & p.val_Intra< 0.05 ~ 'DOWN',
                                              TRUE ~ 'X'))
    
    # Make a df with the clusters and how metabolites are sssigned to them
    Cluster_assignment_df <- data.frame("Cluster" = seq(1,21),
                                        "Intracellular Change"= c(rep("UP",7),rep("No Change",7),rep("DOWN",7)), 
                                        "CoRe Change"= rep(c("Release UP","Release DOWN","Consume UP","Consume DOWN","Released/Consumed UP","Released/Consumed DOWN","No Change"),3))  
    
    Alluvial_DF.final <- DMA
    
    writexl::write_xlsx(Alluvial_DF.final, paste0(Results_folder_MCA_folder,"/",Name,Output_Name,".xlsx", sep = ""))
    write.table(Cluster_assignment_df,row.names =  FALSE, file = paste0(Results_folder_MCA_folder,"/Cluster_assignment_table.csv", sep = "")     )
    
    Alluvial_DF <- Alluvial_DF.final
    
    
  }else{ #  else if(CoRe == FALSE){
    Name = paste0("MCA_Output_",gsub(" ", "_",Condition1),"_with_",gsub(" ", "_",Condition2), sep = "")
    
    C1 <- Input_data
    C1 <- na.omit(C1)
    C1$class <- paste (Condition1)
    C2 <- Input_data2
    C2 <- na.omit(C2)
    C2$class <- paste (Condition2)
    
    #C1 <- C1 %>% drop_na()
    #C2 <- C2 %>% drop_na()
    
    # 0. Overall Regulation: label the selection of metabolites that change or where at least we have a change in one of the two conditions
    if(test== "p.val"){
      C1 <- C1 %>%
        mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                        Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                        TRUE ~ 'No_Change'))
      C2 <- C2%>%
        mutate(MetaboliteChange_Significant1 = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                         Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                         TRUE ~ 'No_Change'))%>%
        rename("class1"="class",
               "Log2FC1"="Log2FC",
               "p.val1"="p.val",
               "p.adj1"="p.adj")
      
    }else{
      C1 <- C1 %>%
        mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                        Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                        TRUE ~ 'No_Change'))
      C2 <- C2%>%
        mutate(MetaboliteChange_Significant1 = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                         Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                         TRUE ~ 'No_Change'))%>%
        rename("class1"="class",
               "Log2FC1"="Log2FC",
               "p.val1"="p.val",
               "p.adj1"="p.adj")
    }
    
    
    
    MergeDF<- merge(C1, C2[,c("Metabolite","MetaboliteChange_Significant1", "class1","Log2FC1","p.val1","p.adj1")], by="Metabolite")%>%
      mutate(Change_Specific = case_when((class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ 'OppositeChange',
                                         (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ 'OppositeChange',
                                         (class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ 'SameDirection_UP',
                                         (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ 'SameDirection_DOWN',
                                         (class== paste(Condition1)& MetaboliteChange_Significant == "No_Change")& (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ paste("ChangeOnly", Condition2, "DOWN", sep="_"),
                                         (class== paste(Condition1)& MetaboliteChange_Significant == "No_Change")& (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ paste("ChangeOnly", Condition2, "UP", sep="_"),
                                         (class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "No_Change") ~ paste("ChangeOnly", Condition1, "UP", sep="_"),
                                         (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "No_Change") ~ paste("ChangeOnly", Condition1, "DOWN", sep="_"),
                                         TRUE ~ 'SameDirection_NoChange'))
    
    
    MergeDF_C1 <-MergeDF %>% select(-c( "MetaboliteChange_Significant1", "class1", "Log2FC1", "p.val1", "p.adj1")) %>%
      unite(col=UniqueID, c(Metabolite, MetaboliteChange_Significant, class), sep = "_", remove = FALSE, na.rm = FALSE)%>%
      mutate(Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                         Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP',
                                         Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_DOWN',
                                         Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                         Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~paste("ChangeOnly", Condition2, "DOWN", sep="_"),
                                         Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") ~paste("ChangeOnly", Condition2, "UP", sep="_"),
                                         Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~paste("Unique", Condition1,"DOWN", sep="_"),
                                         Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~paste("Unique", Condition1, "UP", sep="_"),
                                         TRUE ~ paste("Unique", Condition1, sep="_")))
    
    MergeDF_C2 <- MergeDF %>% select(-c( "MetaboliteChange_Significant", "class", "Log2FC", "p.val", "p.adj")) %>%
      unite(col=UniqueID, c(Metabolite, MetaboliteChange_Significant1, class1), sep = "_", remove = FALSE, na.rm = FALSE)%>%
      rename("class"="class1",
             "MetaboliteChange_Significant"="MetaboliteChange_Significant1",
             "Log2FC"="Log2FC1",
             "p.val"="p.val1",
             "p.adj"="p.adj1")%>%
      mutate(Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                         Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP',
                                         Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_DOWN',
                                         Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                         Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~paste("ChangeOnly", Condition1, "DOWN", sep="_"),
                                         Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~paste("ChangeOnly", Condition1, "UP", sep="_"),
                                         Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~paste("Unique", Condition2,"DOWN", sep="_"),
                                         Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") ~paste("Unique", Condition2, "UP", sep="_"),
                                         TRUE ~ paste("Unique", Condition1, sep="_")))
    
    
    Alluvial_DF <- rbind(MergeDF_C1, MergeDF_C2)
    Alluvial_DF<- Alluvial_DF%>%
      mutate(Amount_Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ paste((sum(Alluvial_DF$Change_Specific=="OppositeChange", na.rm=T))/2),
                                                Change_Specific== 'SameDirection_UP' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_UP", na.rm=T))/2),
                                                Change_Specific== 'SameDirection_DOWN' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_DOWN", na.rm=T))/2),
                                                Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition1, "UP", sep="_"), na.rm=T)),
                                                Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition1, "DOWN", sep="_"), na.rm=T)),
                                                Change_Specific== paste("ChangeOnly", Condition2, "UP" ,sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition2, "UP", sep="_"), na.rm=T)),
                                                Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition2, "DOWN", sep="_"), na.rm=T)),
                                                Change_Specific== 'SameDirection_NoChange' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_NoChange", na.rm=T))/2),
                                                Change_Specific== paste("Unique", Condition1,"DOWN", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition1,"DOWN", sep="_"), na.rm=T)),
                                                Change_Specific== paste("Unique", Condition1,"UP", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition1,"UP", sep="_"), na.rm=T)),
                                                Change_Specific== paste("Unique", Condition2,"DOWN", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition2,"DOWN", sep="_"), na.rm=T)),
                                                Change_Specific== paste("Unique", Condition2,"UP", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition2,"UP", sep="_"), na.rm=T)),
                                                TRUE ~ 'FALSE'))
    Alluvial_DF<- Alluvial_DF%>%
      mutate(Overall_Change = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                        Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP_or_DOWN',
                                        Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_UP_or_DOWN',
                                        Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                        Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition2, sep="_"),
                                        Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition2, sep="_"),
                                        Change_Specific== paste("Unique", Condition1, "DOWN",sep="_") ~paste("Unique", Condition1, sep="_"),
                                        Change_Specific== paste("Unique", Condition1, "UP",sep="_") ~paste("Unique", Condition1, sep="_"),
                                        Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition1, sep="_"),
                                        Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition1, sep="_"),
                                        Change_Specific== paste("Unique", Condition2,"DOWN", sep="_") ~paste("Unique", Condition2, sep="_"),
                                        Change_Specific== paste("Unique", Condition2,"UP", sep="_") ~paste("Unique", Condition2, sep="_"),
                                        TRUE ~ "FALSE"))
    Alluvial_DF <- Alluvial_DF %>%
      mutate(Amount_Overall_Change = case_when(Overall_Change== 'OppositeChange' ~ paste((sum(Alluvial_DF$Overall_Change=="OppositeChange", na.rm=T))/2),
                                               Overall_Change== 'SameDirection_UP_or_DOWN' ~ paste((sum(Alluvial_DF$Overall_Change=="SameDirection_UP_or_DOWN", na.rm=T))/2),
                                               Overall_Change== paste("Unique", Condition1, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("Unique", Condition1, sep="_"), na.rm=T)),
                                               Overall_Change== paste("Unique", Condition2, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("Unique", Condition2, sep="_"), na.rm=T)),
                                               Overall_Change== paste("ChangeOnly", Condition1, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("ChangeOnly", Condition1, sep="_"), na.rm=T)),
                                               Overall_Change== paste("ChangeOnly", Condition2, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("ChangeOnly", Condition2, sep="_"), na.rm=T)),
                                               Overall_Change== 'SameDirection_NoChange' ~ paste((sum(Alluvial_DF$Overall_Change=="SameDirection_NoChange", na.rm=T))/2),
                                               TRUE ~ 'FALSE'))
    
    # Create Alluvial final output df
    C1.final<- Alluvial_DF %>% filter(class == Condition1)
    C1.final <- C1.final %>% select(-c("UniqueID", "class", "Change_Specific","Amount_Change_Specific", "Overall_Change", "Amount_Overall_Change" ))
    C2.final<- Alluvial_DF %>% filter(class == Condition2)
    C2.final <- C2.final %>% select(-c("UniqueID", "class" ))
    Alluvial_DF.final <- merge(C1.final, C2.final, by = "Metabolite")
    names(Alluvial_DF.final) <- gsub(".x",paste(".",Condition1, sep = ""),names(Alluvial_DF.final))
    names(Alluvial_DF.final) <- gsub(".y",paste(".",Condition2, sep = ""),names(Alluvial_DF.final))
    names(Alluvial_DF.final) <- gsub(x = names(Alluvial_DF.final), pattern = "MetaboliteChange_Significant", replacement =  paste("MetaboliteChange_Significant_",test,pCutoff,"logFC",FCcutoff, sep = ""))
    
    ##Write to file
    writexl::write_xlsx(Alluvial_DF.final, paste0(Results_folder_MCA_folder,"/MCA_Output_",Name,Output_Name,".xlsx", sep = ""))
    # write.csv(Alluvial_DF2, paste("AlluvianDF", Output, ".csv", sep="_"), row.names= TRUE)
  }
  
  if(plot == TRUE){
    
    ###########################
    ### Make Alluvial plot ###
    
    #Add Frequency:
    Alluvial_DF[,"Frequency"]  <- as.numeric("1")
    Alluvial_Plot <- Alluvial_DF
    
    ### select Plot parameters
    if (CoRe == TRUE){
      plot_column_names= c("Intracellular Change", "CoRe Change", "Metabolite Cluster", "Metabolite")
      plot_color_variable = "Metabolite Cluster"
      # Remove the No Change No Change from the plot
      Alluvial_Plot <- Alluvial_Plot[!c(Alluvial_Plot$`Intracellular Change`=="No Change" & Alluvial_Plot$`CoRe Change`== "No Change"),]
    }else{
      plot_column_names= c("class", "MetaboliteChange_Significant", "Overall_Change", "Metabolite")
      plot_color_variable = "Overall_Change"
      Alluvial_Plot<- Alluvial_Plot[-which(Alluvial_Plot[,plot_color_variable]=="SameDirection_NoChange"),]
    }
    safe_colorblind_palette = c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882255",  "#6699CC", "#117733", "#888888","red", "white", "#000") # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible,
    
    if(Save_as=="pdf"){
    pdf(paste(Results_folder_MCA_folder,"/AlluvianPlot", Name,Output_Name, ".",Save_as, sep=""), width=12, height=9) # Save_as_var
    }else{
      svg(paste(Results_folder_MCA_folder,"/AlluvianPlot", Name,Output_Name, ".",Save_as, sep=""), width=12, height=9) # Save_as_var
    }
    
    
    par(oma=c(2,2,8,2), mar = c(2, 2, 0.1, 2)+0.1)#https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
    alluvial::alluvial( Alluvial_Plot %>% select(all_of(plot_column_names)), freq=Alluvial_Plot$Frequency,
                        col = case_when(Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[1] ~ safe_colorblind_palette[1],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[2] ~ safe_colorblind_palette[2],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[3] ~ safe_colorblind_palette[3],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[4] ~ safe_colorblind_palette[4],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[5] ~ safe_colorblind_palette[5],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[6] ~ safe_colorblind_palette[6],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[7] ~ safe_colorblind_palette[7],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[8] ~ safe_colorblind_palette[8],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[9] ~ safe_colorblind_palette[9],
                                        Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[10] ~ safe_colorblind_palette[10],
                                        TRUE ~ 'black'),
                        border = case_when(Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[1] ~ safe_colorblind_palette[1],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[2] ~ safe_colorblind_palette[2],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[3] ~ safe_colorblind_palette[3],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[4] ~ safe_colorblind_palette[4],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[5] ~ safe_colorblind_palette[5],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[6] ~ safe_colorblind_palette[6],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[7] ~ safe_colorblind_palette[7],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[8] ~ safe_colorblind_palette[8],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[9] ~ safe_colorblind_palette[9],
                                           Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[10] ~ safe_colorblind_palette[10],
                                           
                                           TRUE ~ 'black'),
                        hide = Alluvial_Plot$Frequency == 0,
                        cex = 0.3,
                        cex.axis=0.5)
    mtext("Selection of metabolites that change in at least one of the two conditions", side=3, line=6, cex=1.2, col="black", outer=TRUE) #https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
    mtext(paste("",Name), side=3, line=5, cex=0.8, col="black", outer=TRUE)
    mtext(paste("Underlying comparison: ",Condition1,"-versus-",Condition2), side=2, line=0, cex=0.8, col="black", outer=TRUE)
    
    if(CoRe==FALSE){
      #mtext("Legend", side=3, line=5, adj=1.0, cex=1, col="black", outer=TRUE)
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[1])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[1]), side=3, line=6, adj=0, cex=0.6, col=safe_colorblind_palette[1], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[2])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[2]), side=3, line=5, adj=0, cex=0.6, col=safe_colorblind_palette[2], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[3])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[3]), side=3, line=4, adj=0, cex=0.6, col=safe_colorblind_palette[3], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[4])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[4]), side=3, line=3, adj=0, cex=0.6, col=safe_colorblind_palette[4], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[5])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[5]), side=3, line=2, adj=0, cex=0.6, col=safe_colorblind_palette[5], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[6])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[6]), side=3, line=1, adj=0, cex=0.6, col=safe_colorblind_palette[6], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[7])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[7]), side=3, line=0, adj=0, cex=0.6, col=safe_colorblind_palette[7], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[8])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[8]), side=3, line=7, adj=1, cex=0.6, col=safe_colorblind_palette[8], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[9])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[9]), side=3, line=6, adj=1, cex=0.6, col=safe_colorblind_palette[9], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[10])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[10]), side=3, line=5, adj=1, cex=0.6, col=safe_colorblind_palette[10], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[11])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[11]), side=3, line=4, adj=1, cex=0.6, col=safe_colorblind_palette[11], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[12])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[12]), side=3, line=3, adj=1, cex=0.6, col=safe_colorblind_palette[12], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[13])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[13]), side=3, line=2, adj=1, cex=0.6, col=safe_colorblind_palette[13], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[14])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[14]), side=3, line=1, adj=1, cex=0.6, col=safe_colorblind_palette[14], outer=TRUE)
      }
      if( is.na(unique(Alluvial_Plot[,plot_color_variable])[15])==FALSE)  {
        mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[15]), side=3, line=7, adj=1, cex=0.6, col=safe_colorblind_palette[15], outer=TRUE)
      }
    } 
    dev.off()# Close the pdf file
    
    
  }
}

