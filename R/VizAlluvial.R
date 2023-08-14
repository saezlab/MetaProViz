## ---------------------------
##
## Script name: Visualization
##
## Purpose of script: Data Visualisation of the MetaProViz analysis to aid biological interpretation
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
#'
#' This script allows you to perform different data visualizations using the results of the MetaProViz analysis

######################################
### ### ### Alluvial Plots ### ### ###
######################################=
#' Notes; Check the alluvial in the Metabolic Clusters. This one is a bit old but kept here to return to if nedded

VizAlluvial <- function(Input_data1,
                     Input_data2,
                     Output_Name = "Metabolic_Clusters_Output_Condition1-versus-Condition2",
                     Condition1,
                     Condition2,
                     pCutoff= 0.05 ,
                     FCcutoff=0.5,
                     test = "p.adj",
                     plot_column_names= c("class", "MetaboliteChange_Significant", "Overall_Change", "Metabolite"),
                     safe_colorblind_palette = c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882255",  "#6699CC", "#117733", "#888888","#CC6677", "#FFF", "#000"), # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
                     plot_color_variable = "Overall_Change",
                     plot_color_remove_variable = "SameDirection_NoChange",
                     Save_as_Plot = pdf
                     ){


  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "alluvial")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages)
  }
  suppressMessages(library(tidyverse))

  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/MetaProViz_Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_MetabolicCluster_folder = paste(Results_folder,"/Alluvial",  sep="")
  if (!dir.exists(Results_folder_plots_MetabolicCluster_folder)) {dir.create(Results_folder_plots_MetabolicCluster_folder)}



  #####################################################
  ### ### ### make output plot Save_as_Plot name ### ### ###
  Save_as_Plot_var <- Save_as_Plot
  Save_as_Plot= deparse(substitute(Save_as_Plot))


  C1 <- Input_data1
  C1 <- na.omit(C1)
  C1$class <- paste (Condition1)
  C2 <- Input_data2
  C2 <- na.omit(C2)
  C2$class <- paste (Condition2)

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
  names(Alluvial_DF.final) <- gsub(".x",paste(".",substr(Condition1, 1, 3), sep = ""),names(Alluvial_DF.final))
  names(Alluvial_DF.final) <- gsub(".y",paste(".",substr(Condition2, 1, 3), sep = ""),names(Alluvial_DF.final))
  names(Alluvial_DF.final) <- gsub(x = names(Alluvial_DF.final), pattern = "MetaboliteChange_Significant", replacement =  paste("MetaboliteChange_Significant_",test,pCutoff,"logFC",FCcutoff, sep = ""))

  ##Write to file
  # This is not needed fot the plots
  #  writexl::write_xlsx(Alluvial_DF.final, paste("Results_", Sys.Date(), "/MetabolicCluster_plots/","Metabolic_Clusters_Output_",Condition1,"-versus-",Condition2,Output_Name,".xlsx", sep = ""))
  # write.csv(Alluvial_DF2, paste("AlluvianDF", Output, ".csv", sep="_"), row.names= TRUE)


  # 1. Regulation:
  Alluvial_DF2  <- Alluvial_DF  %>%
    mutate(MetaboliteChange = case_when(Log2FC >= FCcutoff  ~ 'UP',
                                        Log2FC <= -FCcutoff ~ 'DOWN',
                                        TRUE ~ 'No_Change'))

  if (test == "p.val"){
    # 2. Excluded according to p-value:
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(Excluded_by_pval = case_when(p.val <= pCutoff ~ 'NO',
                                          p.val > pCutoff ~ 'YES'))
    #3. Excluded according to P-value?
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Excluded = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                   Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                   Log2FC >= FCcutoff & p.val >= pCutoff ~ 'UP_Excluded',
                                                   Log2FC <= -FCcutoff & p.val >= pCutoff ~ 'DOWN_Excluded',
                                                   TRUE ~ 'No_Change'))
    #4. After exlusion by p-value
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
  }else if(test=="p.adj"){
    # 2. Excluded according to p-value:
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(Excluded_by_pval = case_when(p.adj <= pCutoff ~ 'NO',
                                          p.adj > pCutoff ~ 'YES'))
    #3. Excluded according to P-value?
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Excluded = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                   Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                   Log2FC >= FCcutoff & p.adj >= pCutoff ~ 'UP_Excluded',
                                                   Log2FC <= -FCcutoff & p.adj >= pCutoff ~ 'DOWN_Excluded',
                                                   TRUE ~ 'No_Change'))
    #4. After exlusion by p-value
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
  }


  #5. Frequency:
  Alluvial_DF2[,"Frequency"]  <- as.numeric("1")
  #6. Safe DF
  # Alluvial_DF3 <- Alluvial_DF2[, c(1:6,12:14,7:8,10,9,11,15)]
  # Alluvial_Plot <- Alluvial_DF3[,c(6:14,2,15)]

  Alluvial_Plot <- Alluvial_DF2

  #Make the SelectionPlot:
  if(plot_color_remove_variable %in% Alluvial_Plot[,plot_color_variable] ){
    Alluvial_Plot<- Alluvial_Plot[-which(Alluvial_Plot[,plot_color_variable]==plot_color_remove_variable),]#remove the metabolites that do not change in either of the conditions
  }


  Save_as_Plot_var(paste(Results_folder_plots_MetabolicCluster_folder,"/Metabolic_Clusters_",Condition1,"-versus-",Condition2,"_", Output_Name,  ".",Save_as_Plot, sep=""), width=12, height=9)
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
  mtext(paste("",Output_Name), side=3, line=5, cex=0.8, col="black", outer=TRUE)
  mtext(paste("Underlying comparison: ",Condition1,"-versus-",Condition2), side=2, line=0, cex=0.8, col="black", outer=TRUE)
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
  dev.off()# Close the pdf file

}
