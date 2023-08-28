## ---------------------------
##
## Script name: Visualization Lollipop graph
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
#' This script allows you to perform Lollipop visualizations using the results of the MetaProViz analysis


#####################################
### ### ### Lolipop Plots ### ### ###
#####################################
#' @param Plot_Settings \emph{Optional: } Choose between "Standard" (Input_data), "Compare" (plot two comparisons together Input_data and Input_data2) or "PEA" (Pathway Enrichment Analysis) \strong{Default = "Standard"}
#' @param Plot_SettingsInfo \emph{Optional: } NULL or Named vector including at least one of those three information for Plot_Settings="Standard" or "Compare": c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile", individual="ColumnName_Plot_SettingsFile"). For Plot_Settings="PEA" a named vector with c(PEA_Pathway="ColumnNameAdditionalInput_data", PEA_score="ColumnNameAdditionalInput_data", PEA_stat= "ColumnNameAdditionalInput_data", individual="Plot_SettingsFile), optionally you can additionally include c(color="ColumnName_Plot_SettingsFile", shape= "ColumnName_Plot_SettingsFile").\strong{Default = NULL}
#' @param Plot_SettingsFile \emph{Optional: } DF with column "Metabolite" including the Metabolite names (needs to match Metabolite names of Input_data) and other columns with required PlotSettingInfo. \strong{Default = NULL}
#' @param Input_data DF with column "Metabolite" including the Metabolite names, Log2FC, pvalue/padjusted values. Can also include additional columns with metadata usable for Plot_Setting_Info. Multiple Input data frames can be added using a list ie. list(df1, df2, df3)
#' @param OutputPlotName \emph{Optional: } String which is added to the output files of the plot. \strong{Default = ""}
#' @param test \emph{Optional: } String which selects pvalue or padj for significance. \strong{Default = padj}
#' @param Comparison_name \emph{Optional: } List including those information about the datasets that are compared on the plots when choosing Plot_Settings= "Compare". \strong{Default = list("Cond1", "Cond2")}
#' @param pCutoff \emph{Optional: } Number of the desired p value cutoff for assessing significance. \strong{Default = 0.05}
#' @param FCcutoff \emph{Optional: } Number of the desired log fold change cutoff for assessing significance. \strong{Default = 0.5}
#' @param Legend \emph{Optional: } Legend=="Pie" will plot a PieChart as the legend for color or Legend="Standard, plot the standard legend for color. \strong{Default = "Standard"}
#' @param Subtitle \emph{Optional: } \strong{Default = ""}
#' @param Theme \emph{Optional: } Selection of theme for plot. \strong{Default = NULL}
#'
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#'
#' @keywords Lolipop plot, pathways
#' @export

VizLolipop<- function(Plot_Settings="Standard",
                      Plot_SettingsInfo= NULL,
                      Plot_SettingsFile= NULL, # Input_pathways = NULL,
                      Input_data, # a dataframe of list of dataframes
                      AdditionalInput_data= NULL, # used only for PEA
                      x = "Log2FC",
                      y = "Metabolite",
                      OutputPlotName= "",
                      Comparison_name= c(Input_data="Cond1", AdditionalInput_data= "Cond2"),
                      Subtitle= "",
                      Theme= NULL,
                      Save_as_Plot = "svg",
                      parameter_size="Reverse" #or default "Standard"
){
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "showtext", "cowplot")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))

  ## ------------ Check Input files ----------- ##
  # 1. The input data:
  if(inherits(Input_data,"list") ==FALSE){
    temp <- list(Input_data)
    Input_data <- temp
  }
  if(inherits(Plot_SettingsFile,"list") ==FALSE){
    temp <- list(Plot_SettingsFile)
    Plot_SettingsFile <- temp
  }
  Plot_SettingsFile_List <- Plot_SettingsFile
  if(length(Plot_SettingsFile_List)==1){
    Plot_SettingsFile <- Plot_SettingsFile_List[[1]]
  }

  Flip <- c() # Vector to check if x or y in numeric
  for(i in 1:length(Input_data)){
    # i=2
    data <- Input_data[[i]]
    if(any(duplicated(row.names(data)))==TRUE){
      stop("Duplicated row.names of Input_data, whilst row.names must be unique")
    }
    if("Metabolite" %in% colnames(data) == FALSE){
      stop("Check input. Input_data must contain a column named `Metabolite` including the metabolite names.")
    }
    # Check if the next lines work correctly in case of duplicated metabolites
    if(length(data[duplicated(data$Metabolite), "Metabolite"]) > 0){
      doublons <- as.character(data[duplicated(data$Metabolite), "Metabolite"])#number of duplications
      data <-data[!duplicated(data$Metabolite),]#remove duplications
      warning("Input_data contained duplicates based on Metabolite! Dropping duplicate IDs and kept only the first entry. You had ", length(doublons), " duplicates. Note that you should do this before running VizLolipop.")
    }
    if(paste(x) %in% colnames(data)==FALSE | paste(y) %in% colnames(data)==FALSE){
      stop("Check your input. The column name of x and/or y does not exist in Input_data.")
    }
    if (is.numeric(data[[x]]) && is.character(data[[y]])) {
      Flip[i] <-FALSE
    } else if (is.character(data[[x]]) && is.numeric(data[[y]])) {
      Flip[i] <- TRUE
    } else {
      stop("One of the x or y must by numeric and the other must be a character")
    }
  }

  if (sum(Flip) == length(Input_data)) {
    Flip=TRUE
    temp<- x
    x<-y
    y<- temp
  } else{
    Flip <- FALSE
  }

  # 2. The Plot_settings: Plot_Settings, Plot_SettingInfo and Plot_SettingFile
  Plot_options <- c("Standard", "Compare", "PEA")
  if (Plot_Settings %in% Plot_options == FALSE){
    stop("Plot_Settings option is incorrect. The allowed options are the following: ",paste(Plot_options, collapse = ", "),"." )
  }
  if(is.vector(Plot_SettingsInfo)==TRUE & is.null(Plot_SettingsFile)==TRUE){
    stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile.")
  }
  if(Plot_Settings=="Compare" & "color" %in% names(Plot_SettingsInfo)==TRUE){
    warning("When Plot_Settings='Compare' color is used to colorcode the input datasets. Color was added in the Plot_SettingsInfo but it will be ignored. Please use only the size option.")
    Plot_SettingsInfo <- Plot_SettingsInfo[!names(Plot_SettingsInfo) == "color"]
  }

  if("color" %in% names(Plot_SettingsInfo)==TRUE){
    for(i in 1:length(Plot_SettingsFile_List)){
      data <- Plot_SettingsFile_List[[i]]
      if (Plot_SettingsInfo[["color"]] %in% colnames(data) == FALSE) {
        stop("You have chosen color = ",paste(Plot_SettingsInfo[["color"]]), ", ", paste(Plot_SettingsInfo[["color"]])," does not exist in the PlotSettingsFile."   )
      }
    }
  }
  if("size" %in% names(Plot_SettingsInfo)==TRUE){
    for(i in 1:length(Plot_SettingsFile_List)){
      data <- Plot_SettingsFile_List[[i]]
      if (Plot_SettingsInfo[["size"]] %in% colnames(data) == FALSE) {
        stop("You have chosen size = ",paste(Plot_SettingsInfo[["size"]]), ", ", paste(Plot_SettingsInfo[["size"]])," does not exist in the PlotSettingsFile."   )
      }
    }
  }
  if("label_dot" %in% names(Plot_SettingsInfo)==TRUE){
    for(i in 1:length(Plot_SettingsFile_List)){
      data <- Plot_SettingsFile_List[[i]]
      if (Plot_SettingsInfo[["label_dot"]] %in% colnames(data) == FALSE) {
        stop("You have chosen label_dot = ",paste(Plot_SettingsInfo[["label_dot"]]), ", ", paste(Plot_SettingsInfo[["label_dot"]])," does not exist in the PlotSettingsFile."   )
      }
    }
  }

  if(is.vector(Plot_SettingsInfo)==TRUE){
    if("color" %in% names(Plot_SettingsInfo)==TRUE & "size" %in% names(Plot_SettingsInfo)==TRUE){
      if((Plot_SettingsInfo[["size"]] == Plot_SettingsInfo[["color"]])==TRUE){
        Plot_SettingsFile$size <- Plot_SettingsFile[,paste(Plot_SettingsInfo[["color"]])]
        Plot_SettingsFile<- Plot_SettingsFile%>%
          dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
      }
      if((Plot_SettingsInfo[["size"]] == Plot_SettingsInfo[["color"]])==FALSE & "color" %in% names(Plot_SettingsInfo)==TRUE){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
      }
      if((Plot_SettingsInfo[["size"]] == Plot_SettingsInfo[["color"]])==FALSE & "size" %in% names(Plot_SettingsInfo)==TRUE){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("size"=paste(Plot_SettingsInfo[["size"]]))
      }
      Plot_SettingsFile2 <- Plot_SettingsFile%>% select(Metabolite, color, size)########################################################
    } else if("color" %in% names(Plot_SettingsInfo)==TRUE & "size" %in% names(Plot_SettingsInfo)==FALSE){
      Plot_SettingsFile <- Plot_SettingsFile%>%
        dplyr::rename("color"=paste(Plot_SettingsInfo[["color"]]))
      Plot_SettingsFile2 <- Plot_SettingsFile%>% select(Metabolite, color)##########################################################
    } else if("color" %in% names(Plot_SettingsInfo)==FALSE & "size" %in% names(Plot_SettingsInfo)==TRUE){
      if(length(Plot_SettingsFile_List)==1){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("size"=paste(Plot_SettingsInfo[["size"]]))
        Plot_SettingsFile2 <- Plot_SettingsFile%>% select(Metabolite,size) ##############################
      }else if(length(Plot_SettingsFile_List)>1){
        for(i in 1:length(Plot_SettingsFile)){
          file <- Plot_SettingsFile[[i]]
          file <- file%>%
            dplyr::rename("size"=paste(Plot_SettingsInfo[["size"]]))
          Plot_SettingsFile_List[[i]] <- file %>% select(Metabolite, size)
          Plot_SettingsFile2 <- Plot_SettingsFile_List[[i]] # Needs fix This is not needed it is here for the code not to break
        }
      }
    }


    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      if(length(Plot_SettingsFile_List)==1){
        Plot_SettingsFile <- Plot_SettingsFile%>%
          dplyr::rename("individual"=paste(Plot_SettingsInfo[["individual"]]))
        if("color" %in% names(Plot_SettingsInfo) | "size" %in% names(Plot_SettingsInfo)){
          Plot_SettingsFile3 <- Plot_SettingsFile %>% select(individual, Metabolite)
          Plot_SettingsFile <- merge(Plot_SettingsFile2, Plot_SettingsFile3, by="Metabolite")
        }
      }else if(length(Plot_SettingsFile_List)>1){
        for(i in 1:length(Plot_SettingsFile)){
          file <- Plot_SettingsFile[[i]]
          file <- file%>%
            dplyr::rename("individual"=paste(Plot_SettingsInfo[["individual"]]))

          if("size" %in% names(Plot_SettingsInfo)){
            Plot_SettingsFile_List[[i]] <- merge(Plot_SettingsFile_List[[i]] ,file %>% select(Metabolite, individual), by = "Metabolite")
            Plot_SettingsFile2 <- Plot_SettingsFile_List[[i]] # Needs fix This is not needed it is here for the code not to break

          }else{
            Plot_SettingsFile_List[[i]] <- file %>% select(Metabolite, individual)
            Plot_SettingsFile2 <- Plot_SettingsFile_List[[i]] # Needs fix This is not needed it is here for the code not to break
          }
        }
      }
    }else{
      Plot_SettingsFile <- Plot_SettingsFile2
    }
  }else if(is.vector(Plot_SettingsInfo)==FALSE & is.null(Plot_SettingsInfo)==FALSE){
    stop("Plot_SettingsInfo must be named vector or NULL.")
  }


  if(Plot_Settings=="PEA" & is.vector(Plot_SettingsInfo)==FALSE){
    stop("You have chosen Plot_Settings=`PEA` that requires you to provide a vector for Plot_SettingsInfo.")
  }else if(Plot_Settings=="PEA" & is.null(Plot_SettingsFile)==TRUE){
    stop("You have chosen Plot_Settings=`PEA` that requires you to provide a DF Plot_SettingsFile including the pathways used for the enrichment analysis.")
  } else if(Plot_Settings=="PEA" & is.null(Plot_SettingsFile)==FALSE & is.null(Plot_SettingsFile)==FALSE){
    if("individual" %in% names(Plot_SettingsInfo)==FALSE | "PEA_score" %in% names(Plot_SettingsInfo)==FALSE | "PEA_stat" %in% names(Plot_SettingsInfo)==FALSE | "PEA_Pathway" %in% names(Plot_SettingsInfo)==FALSE){
      stop("You have chosen Plot_Settings=`PEA` that requires you to provide a vector for Plot_SettingsInfo including `individual`, `PEA_Pathway`, `PEA_stat` and `PEA_score`.")
    }
  }

  # The next lines need checks/corrections
  if(Plot_Settings=="PEA" & is.null(AdditionalInput_data)==TRUE){
    stop("If Plot_Settings=`PEA` you have to provide a DF for AdditionalInput_data including the results of an enrichment analysis.")
  } else if(Plot_Settings=="PEA" & is.null(AdditionalInput_data)==FALSE){
    AdditionalInput_data <- AdditionalInput_data%>%
      dplyr::rename("PEA_score"=paste(Plot_SettingsInfo[["PEA_score"]]),
                    "PEA_stat"=paste(Plot_SettingsInfo[["PEA_stat"]]),
                    "PEA_Pathway"=paste(Plot_SettingsInfo[["PEA_Pathway"]]))
  }

  # 4. Check other plot-specific parameters:
  if (!is.null(Save_as_Plot)) {
    Save_as_Plot_options <- c("svg","pdf","png")
    if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the following: ",paste(Save_as_Plot_options,collapse = ", "),"." )
    }
  }

  #Add the theme

  ## ------------ Create Output folders ----------- ##\
  if (!is.null(Save_as_Plot)) {
    name <- paste0("MetaProViz_Results_",Sys.Date())
    WorkD <- getwd()
    Results_folder <- file.path(WorkD, name)
    if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
    Results_folder_plots_Lolipop_folder = file.path(Results_folder, "Lolipop")  # This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
    if (!dir.exists(Results_folder_plots_Lolipop_folder)) {dir.create(Results_folder_plots_Lolipop_folder)}  # check and create folder
  }
  ##########################################################################################################################################
  #--------------Plots:
  if(Plot_Settings=="Standard"){############################################################################################################
    Input_data <- Input_data[[1]]
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){
      # Create the list of individual plots that should be made:
      IndividualPlots <- Plot_SettingsFile[!duplicated(Plot_SettingsFile$individual),]
      IndividualPlots <- IndividualPlots$individual

      PlotList <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        # i = IndividualPlots[1]
        Plot_SettingsInfo_indi <- Plot_SettingsInfo
        Plot_SettingsFile_Select <- subset(Plot_SettingsFile, individual == paste(i))
        InputLolipop  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
          na.omit()

        #Select metabolites for the cut offs selected
        loli.data <- InputLolipop %>% mutate(names=Metabolite)

        if("size" %in% names(Plot_SettingsInfo_indi)==TRUE ){
          if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
            stop("Size can take only numeric values")
          }else{# color = continuous
            keyvalssize <- loli.data$size
          }
        } else{
          Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi,size="p.adj")
          keyvalssize <- loli.data$size
        }


        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          label <-  Plot_SettingsInfo[["label_dot"]]
          loli.data[ Plot_SettingsInfo[["label_dot"]]] <- round(loli.data[[label]], digits = 3)
        }else{
          label <-  ""
        }

        #Create the plot:
        p2 <- NULL
        if("color" %in% names(Plot_SettingsInfo_indi)==TRUE ){
          if(is.numeric(loli.data$color)==FALSE){ # run if color is discrete

            col_var_name <- Plot_SettingsInfo_indi[['color']]

            position <- which(names(loli.data)=="color" )
            names(loli.data)[position]<-"plot_color_variable"


            loli.data <- loli.data %>%
              arrange(plot_color_variable, get(x),Metabolite)

            loli.data_avg <- loli.data %>%
              arrange(plot_color_variable, get(x), Metabolite) %>%
              mutate(Metab_name = row_number()) %>%
              group_by(plot_color_variable) %>%
              mutate(
                avg = mean(get(x))
              ) %>%
              ungroup() %>%
              mutate(plot_color_variable = factor(plot_color_variable))


            loli_lines <-   loli.data_avg %>%
              arrange(plot_color_variable, Metabolite) %>%
              group_by(plot_color_variable) %>%
              summarize(
                start_x = min(Metab_name) -0.5,
                end_x = max(Metab_name) + 0.5,
                y = 0#unique(avg)
              ) %>%
              pivot_longer(
                cols = c(start_x, end_x),
                names_to = "type",
                values_to = "x"
              ) %>%
              mutate(
                x_group = if_else(type == "start_x", x + .1, x - .1),
                x_group = if_else(type == "start_x" & x == min(x), x_group - .1, x_group),
                x_group = if_else(type == "end_x" & x == max(x), x_group + .1, x_group) )

            #rm(p2)
            p2 <- loli.data_avg %>%
              ggplot(aes(Metab_name, get(x))) + # names in aes ro Metab_name
              geom_hline(
                data = tibble(y = -5:5),
                aes(yintercept = y),
                color = "grey82",
                size = .5 )

            p2 <- p2 + geom_segment(
              aes(
                xend = Metab_name,          # names
                yend = 0,#avg,
                color = plot_color_variable,
                #color = after_scale(colorspace::lighten(color, .2))
              ))

            p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
              geom_line(
                data = loli_lines,
                aes( x_group, y,
                     color = plot_color_variable,
                     #  color = after_scale(colorspace::darken(color, .2))
                ), size = 2.5) +  geom_point(aes(size = keyvalssize, color = plot_color_variable)
                )

            p2 <- p2 +theme(axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()
            )

            if(Flip == TRUE){
              p2 <- p2 +theme(axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())
              lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)

              p2 <- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite,angle = 90, size = 3)

              lab_neg_metab <-  loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
              p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, angle = 90,size = 3)

              if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
                dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
                p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

                dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
                p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

              }

              p2 <- p2+Theme
              p2 <- p2+ labs(color=col_var_name)+
                labs(size=Plot_SettingsInfo[['size']])

              p2 <- p2  + labs(title = paste(OutputPlotName),subtitle = Subtitle)+
                theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                      plot.subtitle = element_text(color = "black", size=10),
                      plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

            }else{
              p2<- p2 + coord_flip()
              p2 <- p2 +theme(axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())
              lab_pos_metab <- loli.data_avg %>% filter(get(x)>0) %>% select(Metabolite, Metab_name, get(x))
              p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite, size = 3)

              lab_neg_metab <- loli.data_avg %>% filter(get(x)<0) %>% select(Metabolite, Metab_name, get(x))
              p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, size = 3)

              if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
                dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
                p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

                dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
                p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

              }

              p2 <- p2+Theme

              p2 <- p2  + labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle) +
                theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                      plot.subtitle = element_text(color = "black", size=10),
                      plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

              p2 <- p2+ labs(color=col_var_name)+
                labs(size=Plot_SettingsInfo_indi[['size']])
            }

            lolipop_plot <-  p2
            if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
              lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
            }

            # Put back the correct name in the data df
            names(loli.data)[position]<- col_var_name

          }else{# color = continuous
            keyvals <- loli.data$color
          }
        }else{
          Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi, color="p.adj")
          keyvals <- loli.data$color
        }

        if(is.null(p2)==TRUE){
          loli.data$names <- as.factor(loli.data$names)
          loli.data[[x]]<- as.numeric(loli.data[[x]])
          loli.data$names <- reorder(loli.data$names, -loli.data[[x]])

          lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names,label=!!label)) +
            geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
            geom_point(aes(colour = keyvals, size = keyvalssize ))   +
            scale_colour_gradient(low = "red", high = "blue")+
            geom_vline(xintercept = 0)+
            Theme+
            labs(color=Plot_SettingsInfo_indi[['color']])+
            labs(size=Plot_SettingsInfo_indi[['size']]) +
            ylab(y)+
            xlab(x)+
            labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle)+
            theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                  plot.subtitle = element_text(color = "black", size=10),
                  plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
            lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
          }
          if(Flip==TRUE){
            lolipop_plot <- lolipop_plot + coord_flip() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
          }
          if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
            lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
          }
        }

        #save plot and get rid of extra signs before saving
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        if (!is.null(Save_as_Plot)) {
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
          }else{
            ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
          }
        }
        ## Store the plot in the 'plots' list
        PlotList[[cleaned_i]] <- lolipop_plot
        plot(lolipop_plot)
      }
      # Return PlotList into the environment to enable the user to view the plots directly
      #assign("LolipopPlots", PlotList, envir=.GlobalEnv)
      # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
      Return <- PlotList
    }
    else if("individual" %in% names(Plot_SettingsInfo)==FALSE){#############################################################################################
      #Assign Data
      if(is.null(Plot_SettingsFile)==FALSE){
        InputLolipop  <- merge(x=Plot_SettingsFile,y=Input_data, by="Metabolite", all.x=TRUE)%>%
          na.omit()
      }else{
        InputLolipop  <- Input_data
      }

      InputLolipop<- InputLolipop %>% drop_na()
      loli.data <- InputLolipop %>% mutate(names=Metabolite)

      #Assign parameters size, label_dot, color
      if("size" %in% names(Plot_SettingsInfo)==TRUE ){
        if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
          stop("Size can take only numeric values")
        }else{# color = continuous
          keyvalssize <- loli.data$size
        }
      } else{
        Plot_SettingsInfo= c(Plot_SettingsInfo,size="p.adj")
        keyvalssize <- loli.data$size
      }

      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        label <-  Plot_SettingsInfo[["label_dot"]]
        loli.data[ Plot_SettingsInfo[["label_dot"]]] <- round(loli.data[[label]], digits = 3)
      }else{
        label <-  ""
      }

      p2 <- NULL
      if("color" %in% names(Plot_SettingsInfo)==TRUE ){
        if(is.numeric(loli.data$color)==FALSE){ # run if color is discrete
          col_var_name <- Plot_SettingsInfo[['color']]

          position <- which(names(loli.data)=="color" )
          names(loli.data)[position]<-"plot_color_variable"

          loli.data <- loli.data %>%
            arrange(plot_color_variable,get(x), Metabolite)

          loli.data_avg <- loli.data %>%
            arrange(plot_color_variable,get(x), Metabolite) %>%
            mutate(Metab_name = row_number()) %>%
            group_by(plot_color_variable) %>%
            mutate(
              avg = mean(get(x))
            ) %>%
            ungroup() %>%
            mutate(plot_color_variable = factor(plot_color_variable))

          loli_lines <-   loli.data_avg %>%
            arrange(plot_color_variable, Metabolite) %>%
            group_by(plot_color_variable) %>%
            summarize(
              start_x = min(Metab_name) -0.5,
              end_x = max(Metab_name) + 0.5,
              y = 0#unique(avg)
            ) %>%
            pivot_longer(
              cols = c(start_x, end_x),
              names_to = "type",
              values_to = "x"
            ) %>%
            mutate(
              x_group = if_else(type == "start_x", x + .1, x - .1),
              x_group = if_else(type == "start_x" & x == min(x), x_group - .1, x_group),
              x_group = if_else(type == "end_x" & x == max(x), x_group + .1, x_group) )

          #rm(p2)
          p2 <- loli.data_avg %>%
            ggplot(aes(Metab_name, get(x))) + # names in aes to Metab_name
            geom_hline(
              data = tibble(y = -5:5),
              aes(yintercept = y),
              color = "grey82",
              size = .5 )

          p2 <- p2 + geom_segment(
            aes(
              xend = Metab_name,          # names
              yend = 0,#avg,
              color = plot_color_variable,
              #color = after_scale(colorspace::lighten(color, .2))
            ))

          p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
            geom_line(
              data = loli_lines,
              aes( x_group, y,
                   color = plot_color_variable,
                   #  color = after_scale(colorspace::darken(color, .2))
              ), size = 2.5) +  geom_point(aes(size = keyvalssize, color = plot_color_variable)
              )

          if(Flip == TRUE){
            p2 <- p2 +theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)
            p2 <- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite,angle = 90, size = 3)

            lab_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, angle = 90,size = 3)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)
            }

            p2 <- p2+Theme
            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo[['size']])

            p2 <- p2  + labs(title = paste(OutputPlotName),subtitle = Subtitle)+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          }else{
            p2<- p2 + coord_flip()
            p2 <- p2 +theme(axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()
            )
            lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite, size = 3)


            lab_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, size = 3)

            # p2 <- p2+ annotate("text", x = max(lab_neg_metab$Metab_name)+ 7, y = 0, label = OutputPlotName, size = 5)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

            }


            p2 <- p2+Theme
            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo[['size']])

            p2 <- p2  + labs(title = paste(OutputPlotName),subtitle = Subtitle)+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          }
          lolipop_plot <-  p2
          if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
            lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
          }

          # Put back the correct name in the data df
          names(loli.data)[position]<- col_var_name

        }else{# color = continuous
          keyvals <- loli.data$color
        }
      }else{
        Plot_SettingsInfo= c(Plot_SettingsInfo, color="p.adj")
        keyvals <- loli.data$color
      }

      if(is.null(p2)==TRUE){
        loli.data$names <- as.factor(loli.data$names)
        loli.data[[x]]<- as.numeric(loli.data[[x]])
        loli.data$names <- reorder(loli.data$names, -loli.data[[x]])

        lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names,label=!!label)) +
          geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
          geom_point(aes(colour = keyvals, size = keyvalssize ))   +
          scale_colour_gradient(low = "red", high = "blue")+
          geom_vline(xintercept = 0)+
          Theme+
          labs(color=Plot_SettingsInfo[['color']])+
          labs(size=Plot_SettingsInfo[['size']]) +
          ylab(y)+
          xlab(x)+
          labs(title = OutputPlotName,subtitle = Subtitle)+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

        if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
          lolipop_plot <- lolipop_plot + scale_size(trans = 'reverse')
        }
        if(Flip==TRUE){
          lolipop_plot <- lolipop_plot + coord_flip() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
        }
      }

      #Add the theme
      if(is.null(Theme)==FALSE){
        lolipop_plot <- lolipop_plot+Theme
      }else{
        lolipop_plot <- lolipop_plot+theme_classic()
      }

      plot(lolipop_plot)

      if (!is.null(Save_as_Plot)) {
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop", OutputPlotName, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=10, height=10)
        }else{
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=10, height=10)
        }
      }
    }
  }else if(Plot_Settings=="Compare"){
    if("individual" %in% names(Plot_SettingsInfo)==TRUE){

      Combined_Input <- data.frame(matrix(ncol = 4, nrow = 0))
      comb.colnames <- c("Metabolite",paste(x), "Condition")
      colnames(Combined_Input) <- comb.colnames

      for (i in 1:length(Input_data)){
        Input_data[[i]]$Condition <- Comparison_name[[i]]
        data <-  Input_data[[i]]
        if(is.null(Plot_SettingsFile)==FALSE){
          data <- merge(Input_data[[i]], Plot_SettingsFile_List[[i]], by = "Metabolite")
        }
        Combined_Input <- rbind(Combined_Input, data)

      }

      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        Plot_SettingsInfo[["label_dot"]]
        Combined_Input[ Plot_SettingsInfo[["label_dot"]]] <- round(Combined_Input[stat], digits = 4)
      }
      Combined_Input <- Combined_Input[is.finite(Combined_Input[[x]]),]


      # Create the list of individual plots that should be made:
      IndividualPlots <- Combined_Input[!duplicated(Combined_Input$individual),]
      IndividualPlots <- IndividualPlots$individual

      PlotList <- list()#Empty list to store all the plots

      for (i in IndividualPlots){
        # i = IndividualPlots[1]
        Plot_SettingsInfo_indi <- Plot_SettingsInfo
        # Plot_SettingsFile_Select <- subset(Combined_Input, individual == paste(i))
        # InputLolipop  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
        #   na.omit()
        InputLolipop <- subset(Combined_Input, individual == paste(i))


        #Select metabolites for the cut offs selected
        loli.data <- InputLolipop %>% mutate(names=Metabolite)



        if("size" %in% names(Plot_SettingsInfo_indi)==TRUE ){
          if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
            stop("Size can take only numeric values")
          }else{# color = continuous
            keyvalssize <- loli.data$size
          }
        } else{
          Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi,size="p.adj")
          keyvalssize <- loli.data$size
        }

        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          Combined_Input[ Plot_SettingsInfo[["label_dot"]]] <- round(Combined_Input[stat], digits = 4)
          label <-  Plot_SettingsInfo[["label_dot"]]
        }else{
          label <-  ""
        }

        if(Flip==TRUE){
          lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names , label=get(label))) +
            geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
            geom_point(aes(colour = Condition, size = keyvalssize ))   +
            #scale_size_continuous(range = c(1,5))+
            geom_vline(xintercept = 0) +
            coord_flip()+
            Theme+
            theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                  plot.subtitle = element_text(color = "black", size=10),
                  plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5),
                  axis.text.x = element_text(angle = 90, hjust = 1))+
            ylab(y)+
            xlab(x)+
            labs(size=Plot_SettingsInfo_indi[['size']])  +
            labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle)

        }else{
          lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names , label=!!label)) +
            geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
            geom_point(aes(colour = Condition, size = keyvalssize ))   +
            # scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +
            #   scale_colour_gradient(low = "red", high = "blue")+#, limits = c(0, max(loli.data[,stat]))) +
            geom_vline(xintercept = 0) +
            Theme+
            theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                  plot.subtitle = element_text(color = "black", size=10),
                  plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))+
            ylab(x)+
            xlab(y)+
            # #  labs(color=Plot_SettingsInfo_indi[['color']]) +
            labs(size=Plot_SettingsInfo_indi[['size']])  +
            labs(title = paste(OutputPlotName,": ", i ),subtitle = Subtitle)
        }
        if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
          lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
        }
        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
        }
        lolipop_plot


        #save plot and get rid of extra signs before saving
        cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
        if (!is.null(Save_as_Plot)) {
          if(OutputPlotName ==""){
            ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
          }else{
            ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
          }
        }
        ## Store the plot in the 'plots' list
        PlotList[[cleaned_i]] <- lolipop_plot
        plot(lolipop_plot)
      }
      # Return PlotList into the environment to enable the user to view the plots directly
      #assign("LolipopPlots", PlotList, envir=.GlobalEnv)
      # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
      Return <- PlotList

    }
    else if("individual" %in% names(Plot_SettingsInfo)==FALSE){
      Combined_Input <- data.frame(matrix(ncol = 4, nrow = 0))
      comb.colnames <- c("Metabolite",paste(x), "Condition")
      colnames(Combined_Input) <- comb.colnames

      for (i in 1:length(Input_data)){
        Input_data[[i]]$Condition <- Comparison_name[[i]]
        data <-  Input_data[[i]]
        if(is.null(Plot_SettingsFile)==FALSE){
          data <- merge(Input_data[[i]], Plot_SettingsFile_List[[i]], by = "Metabolite")
        }
        Combined_Input <- rbind(Combined_Input, data)

      }

      if("size" %in% names(Plot_SettingsInfo)==TRUE ){
        if(is.numeric(Combined_Input$size)==FALSE){ # run is color is discrete
          stop("Size can take only numeric values")
        }else{# color = continuous
          keyvalssize <- Combined_Input$size
        }
      } else{
        Plot_SettingsInfo= c(Plot_SettingsInfo,size="p.adj")
        keyvalssize <- Combined_Input$size
      }

      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        Combined_Input[ Plot_SettingsInfo[["label_dot"]]] <- round(Combined_Input[stat], digits = 4)
        label <-  Plot_SettingsInfo[["label_dot"]]
      }else{
        label <-  ""
      }

      # Remove the metabolite with inf in logFC because it messes the plot
      Combined_Input <- Combined_Input[is.finite(Combined_Input[[x]]),]


      if(Flip==TRUE){
        lolipop_plot <- ggplot(Combined_Input, aes(x=reorder(Metabolite, + get(x)), y=.data[[x]], label=!!label)) +
          geom_point(stat = 'identity', aes(size = keyvalssize, col = Condition))  +
          #   geom_text(color="black", size=2) +
          ylim(((Reduce(min,Combined_Input[[x]]))-0.5),((Reduce(max,Combined_Input[[x]]))+0.5)) +
          geom_hline(yintercept = 0) +
          Theme+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5),
                axis.text.x = element_text(angle = 90, hjust = 1))+
          ylab(x)+
          xlab(y)+
          labs(size=Plot_SettingsInfo[['size']])  +
          labs(title = OutputPlotName,subtitle = Subtitle)

      }else{
        lolipop_plot <- ggplot(Combined_Input, aes(x=reorder(Metabolite,+ get(x)), y=.data[[x]], label=!!label)) +
          geom_point(stat = 'identity', aes(size = keyvalssize, col = Condition))  +
          #geom_text(color="black", size=2) +
          ylim(((Reduce(min,Combined_Input[[x]]))-0.5),((Reduce(max,Combined_Input[[x]]))+0.5)) +
          coord_flip()+
          geom_hline(yintercept = 0) +
          Theme+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))+
          ylab(y)+
          xlab(x)+
          labs(size=Plot_SettingsInfo[['size']])  +
          labs(title = OutputPlotName,subtitle = Subtitle)
      }
      if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
        lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
      }

      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
      }

      if (!is.null(Save_as_Plot)) {
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop", ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=12, height=14)
        }else{
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=12, height=14)
        }
      }
      plot(lolipop_plot)

    }


  }else if(Plot_Settings=="PEA"){# Code Missing
    Input_data <- Input_data[[1]]

    # Create the list of individual plots that should be made:
    IndividualPlots <- Plot_SettingsFile[!duplicated(Plot_SettingsFile$individual),]
    IndividualPlots <- IndividualPlots$individual

    PlotList <- list()#Empty list to store all the plots

    for (i in IndividualPlots){
      # i <- IndividualPlots[1]
      Plot_SettingsInfo_indi <- Plot_SettingsInfo
      Plot_SettingsFile_Select <- subset(Plot_SettingsFile, individual == paste(i))
      InputLolipop  <- merge(x=Plot_SettingsFile_Select,y=Input_data, by="Metabolite", all.x=TRUE)%>%
        na.omit()

      AdditionalInput_data_Select<- subset(AdditionalInput_data, PEA_Pathway == paste(i)) #Select pathway we plot and use the score and stats

      #Select metabolites for the cut offs selected
      loli.data <- InputLolipop %>% mutate(names=Metabolite)


      if("size" %in% names(Plot_SettingsInfo_indi)==TRUE ){
        if(is.numeric(loli.data$size)==FALSE){ # run is color is discrete
          stop("Size can take only numeric values")
        }else{# color = continuous
          keyvalssize <- loli.data$size
        }
      } else{
        Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi,size="p.adj")
        keyvalssize <- loli.data$size
      }


      if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
        label <-  Plot_SettingsInfo[["label_dot"]]
        loli.data[ Plot_SettingsInfo[["label_dot"]]] <- round(loli.data[[label]], digits = 3)
      }else{
        label <-  ""
      }

      p2 <- NULL
      # add color here
      if("color" %in% names(Plot_SettingsInfo_indi)==TRUE ){
        if(is.numeric(loli.data$color)==FALSE){ # run if color is discrete

          col_var_name <- Plot_SettingsInfo_indi[['color']]

          position <- which(names(loli.data)=="color" )
          names(loli.data)[position]<-"plot_color_variable"


          loli.data <- loli.data %>%
            arrange(plot_color_variable, get(x),Metabolite)

          loli.data_avg <- loli.data %>%
            arrange(plot_color_variable, get(x), Metabolite) %>%
            mutate(Metab_name = row_number()) %>%
            group_by(plot_color_variable) %>%
            mutate(
              avg = mean(get(x))
            ) %>%
            ungroup() %>%
            mutate(plot_color_variable = factor(plot_color_variable))


          loli_lines <-   loli.data_avg %>%
            arrange(plot_color_variable, Metabolite) %>%
            group_by(plot_color_variable) %>%
            summarize(
              start_x = min(Metab_name) -0.5,
              end_x = max(Metab_name) + 0.5,
              y = 0#unique(avg)
            ) %>%
            pivot_longer(
              cols = c(start_x, end_x),
              names_to = "type",
              values_to = "x"
            ) %>%
            mutate(
              x_group = if_else(type == "start_x", x + .1, x - .1),
              x_group = if_else(type == "start_x" & x == min(x), x_group - .1, x_group),
              x_group = if_else(type == "end_x" & x == max(x), x_group + .1, x_group) )

          #rm(p2)
          p2 <- loli.data_avg %>%
            ggplot(aes(Metab_name, get(x))) + # names in aes ro Metab_name
            geom_hline(
              data = tibble(y = -5:5),
              aes(yintercept = y),
              color = "grey82",
              size = .5 )

          p2 <- p2 + geom_segment(
            aes(
              xend = Metab_name,          # names
              yend = 0,#avg,
              color = plot_color_variable,
              #color = after_scale(colorspace::lighten(color, .2))
            ))

          p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
            geom_line(
              data = loli_lines,
              aes( x_group, y,
                   color = plot_color_variable,
                   #  color = after_scale(colorspace::darken(color, .2))
              ), size = 2.5) +  geom_point(aes(size = keyvalssize, color = plot_color_variable)
              )

          p2 <- p2 +theme(axis.text.y=element_blank(),
                          axis.ticks.y=element_blank()
          )

          if(Flip == TRUE){
            p2 <- p2 +theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            lab_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x)

            p2 <- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite,angle = 90, size = 3)

            lab_neg_metab <-  loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x)
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, angle = 90,size = 3)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

            }

            p2 <- p2+Theme
            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo[['size']])

            p2 <- p2  +  labs(title = paste(OutputPlotName, ": ", i, sep=""),
                              subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
                              caption = paste0("total = ", nrow(InputLolipop), " metabolites of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"))+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

          }else{
            p2<- p2 + coord_flip()
            p2 <- p2 +theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            lab_pos_metab <- loli.data_avg %>% filter(get(x)>0) %>% select(Metabolite, Metab_name, get(x))
            p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab[[x]]+1.5, label = lab_pos_metab$Metabolite, size = 3)

            lab_neg_metab <- loli.data_avg %>% filter(get(x)<0) %>% select(Metabolite, Metab_name, get(x))
            p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab[[x]]-1.5, label = lab_neg_metab$Metabolite, size = 3)

            if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
              dot_pos_metab <- loli.data_avg[loli.data_avg[x]>0,]  %>% select(Metabolite, Metab_name, x, label)
              p2<- p2+ annotate("text", x = dot_pos_metab$Metab_name, y = dot_pos_metab[[x]], label = dot_pos_metab[[label]], size = 3)

              dot_neg_metab <- loli.data_avg[loli.data_avg[x]<0,]  %>% select(Metabolite, Metab_name, x,label)
              p2<- p2+ annotate("text", x = dot_neg_metab$Metab_name, y = dot_neg_metab[[x]], label = dot_neg_metab[[label]], size = 3)

            }

            p2 <- p2+Theme

            p2 <- p2  +  labs(title = paste(OutputPlotName, ": ", i, sep=""),
                              subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
                              caption = paste0("total = ", nrow(InputLolipop), " metabolites of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"))+
              theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                    plot.subtitle = element_text(color = "black", size=10),
                    plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

            p2 <- p2+ labs(color=col_var_name)+
              labs(size=Plot_SettingsInfo_indi[['size']])
          }

          lolipop_plot <-  p2
          if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
            lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
          }

          # Put back the correct name in the data df
          names(loli.data)[position]<- col_var_name

        }else{# color = continuous
          keyvals <- loli.data$color
        }
      }else{
        Plot_SettingsInfo_indi= c(Plot_SettingsInfo_indi, color="p.adj")
        keyvals <- loli.data$color
      }

      if(is.null(p2)==TRUE){

        loli.data$names <- as.factor(loli.data$names)
        loli.data[[x]]<- as.numeric(loli.data[[x]])
        loli.data$names <- reorder(loli.data$names, -loli.data[[x]])

        lolipop_plot <- ggplot(loli.data , aes(x = get(x), y = names,label=!!label)) +
          geom_segment(aes(x = 0, xend = get(x), y = names, yend = names)) +
          geom_point(aes(colour = keyvals, size = keyvalssize ))   +
          # scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +

          scale_colour_gradient(low = "red", high = "blue")+
          geom_vline(xintercept = 0)+
          Theme+
          labs(color=Plot_SettingsInfo_indi[['color']])+
          labs(size=Plot_SettingsInfo_indi[['size']]) +
          ylab(y)+
          xlab(x)+
          labs(title = paste(OutputPlotName, ": ", i, sep=""),
               subtitle = paste(Plot_SettingsInfo[["PEA_score"]],"= ", AdditionalInput_data_Select$PEA_score, ", ",Plot_SettingsInfo[["PEA_stat"]] , "= ", AdditionalInput_data_Select$PEA_stat, sep=""),
               caption = paste0("total = ", nrow(InputLolipop), " metabolites of ", nrow(Plot_SettingsFile_Select), " metabolites in pathway"))+
          theme(plot.title = element_text(color = "black", size = 12, face = "bold"),
                plot.subtitle = element_text(color = "black", size=10),
                plot.caption = element_text(color = "black",size=9, face = "italic", hjust = 2.5))

        if(parameter_size=="Reverse" & is.null(keyvalssize)==FALSE){
          lolipop_plot <-   lolipop_plot + scale_size(trans = 'reverse')
        }
        if(Flip==TRUE){
          lolipop_plot <- lolipop_plot + coord_flip() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        if("label_dot" %in% names(Plot_SettingsInfo)==TRUE ){
          lolipop_plot <-   lolipop_plot + geom_text(color="black", size=2)
        }
        plot(lolipop_plot)
      }

      #save plot and get rid of extra signs before saving
      cleaned_i <- gsub("[[:space:],/\\\\]", "-", i)#removes empty spaces and replaces /,\ with -
      if (!is.null(Save_as_Plot)) {
        if(OutputPlotName ==""){
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
        }else{
          ggsave(file=paste(Results_folder_plots_Lolipop_folder,"/", "Lolipop_", OutputPlotName, "_",cleaned_i, ".",Save_as_Plot, sep=""), plot=lolipop_plot, width=8, height=6)
        }
      }
      ## Store the plot in the 'plots' list
      PlotList[[cleaned_i]] <- lolipop_plot
      plot(lolipop_plot)
    }
    # Return PlotList into the environment to enable the user to view the plots directly
    #assign("LolipopPlots", PlotList, envir=.GlobalEnv)
    # Combine plots into a single plot using facet_grid or patchwork::wrap_plots
    invisible(PlotList)#returns PlotList if assigned

  }
}




##############################################################
### ### ### lollipop helper function: Internal Function ### ### ###
##############################################################

#' @param Input This is the ggplot object generated within the VizLolipop function.
#'
#' @keywords lollipop helper function
#' @noRd


plotGrob_Lollipop <- function(Input){
  #------- Set the total heights and width
  #we need ggplot_grob to edit the gtable of the ggplot object. Using this we can manipulate the gtable arguments directly.
  plottable<- ggplot2::ggplotGrob(Input) # Convert the plot to a gtable
  if("color" %in% names(Plot_SettingsInfo)==FALSE & "shape" %in% names(Plot_SettingsInfo)==FALSE){
    #-----widths
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(10)] <- unit(0,"cm")#controls margins --> Figure legend
    plottable$widths[c(7,8,9,11)] <- unit(0,"cm")#controls margins --> not needed

    if(is.numeric(Input_data[,paste(x)])==TRUE){
      plottable$widths[5] <- unit(8, "cm")#controls x-axis
      plot_widths <- 11
    }else{
      Features <- nrow(Input_data)/4
      plottable$widths[5] <- unit(Features, "cm")#controls x-axis
      plot_widths <- 11+Features
    }

    #-----heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11,12)] <- unit(0,"cm")#controls margins --> not needed

    if(OutputPlotName=="" & Subtitle==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(1,2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <- 10.5
    } else{
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
      plottable$heights[c(1,2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      plot_heights <-11
    }
  }else if("color" %in% names(Plot_SettingsInfo)==TRUE & "shape" %in% names(Plot_SettingsInfo)==TRUE){
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))+(round(as.numeric(Legend$heights[5]),1))

    #-----Plot widths
    plottable$widths[5] <- unit(8, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed

    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- 11+Value

    #-----Plot heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11)] <- unit(0,"cm")#controls margins --> not needed

    if(OutputPlotName=="" & Subtitle==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed

      if(Legend_heights>10.5){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-10.5)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 10.5
      }
    } else{#If we do have Title and or subtitle
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
      plottable$heights[c(2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      if(Legend_heights>11){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-11)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 11
      }
    }
  }else if("color" %in% names(Plot_SettingsInfo)==TRUE | "shape" %in% names(Plot_SettingsInfo)==TRUE){
    #------- Legend heights
    Legend <- ggpubr::get_legend(Input) # Extract legend to adjust separately
    Legend_heights <- (round(as.numeric(Legend$heights[3]),1))

    #----- Plot widths
    plottable$widths[5] <- unit(8, "cm")#controls x-axis
    plottable$widths[c(3)] <- unit(2,"cm")#controls margins --> y-axis label is there
    plottable$widths[c(1,2,4)] <- unit(0,"cm")#controls margins --> not needed
    plottable$widths[c(6)] <- unit(1,"cm")#controls margins --> start Figure legend
    plottable$widths[c(7,8,10,11)] <- unit(0,"cm")#controls margins --> not needed

    Value <- round(as.numeric(plottable$widths[9]),1) #plottable$widths[9] is a <unit/unit_v2> object and we can extract the extract the numeric part
    plot_widths <- 11+Value

    #-----Plot heigths
    plottable$heights[7] <- unit(8, "cm")#controls x-axis
    plottable$heights[c(8)] <- unit(1,"cm")#controls margins --> x-axis label
    plottable$heights[c(10)] <- unit(1,"cm")#controls margins --> Figure caption
    plottable$heights[c(9,11)] <- unit(0,"cm")#controls margins --> not needed

    if(OutputPlotName=="" & Subtitle==""){
      plottable$heights[c(6)] <- unit(0.5,"cm")#controls margins --> Some space above the plot
      plottable$heights[c(2,3,4,5)] <- unit(0,"cm")#controls margins --> not needed

      if(Legend_heights>10.5){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-10.5)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 10.5
      }
    }else{#If we do have Title and or subtitle
      plottable$heights[c(3)] <- unit(1,"cm")#controls margins --> OutputPlotName and subtitle
      plottable$heights[c(2,4,5,6)] <- unit(0,"cm")#controls margins --> not needed
      if(Legend_heights>11){#If the legend requires more heights than the Plot
        Add <- (Legend_heights-11)/2
        plottable$heights[1] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(Add,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- Legend_heights
      }else{
        plottable$heights[1] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the top
        plottable$heights[12] <- unit(0,"cm")#controls margins --> Can be increased if Figure legend needs more space on the bottom
        plot_heights <- 11
      }
    }
  }
  #plot_param <-c(plot_heights=plot_heights, plot_widths=plot_widths)
  Output<- list(plot_heights, plot_widths, plottable)
}



# Helper function needed for adding column to pathway file defining if this metabolite is unique/multiple pathways

