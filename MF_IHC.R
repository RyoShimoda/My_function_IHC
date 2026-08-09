# My Function Immunohistochemistry cFos

mf_IHC <- function(Folder_name = "", Path, File_type = "csv", 
                         Red_name = "RED", Green_name = "GREEN", Blue_name = "BLUE",
                         Area = "Area", Luminance = "Luminance", Red_intensity = "RED intensity", Count = "Count",
                         Region_list, Analyse_color = "R/G", Aim_color = "R",
                         D_V = "T", PFC = "F", number_of_region = 9, number_of_sample = 16,
                         Font_family = "Times New Roman", graphic_refference = "SED", plot_number = "",
                         title = "", line_size = .8, errobar_width = .2, 
                         legend_text_size = 16, legend_position  = c(.85, .9),
                         point_size = 5, point_shape = 21,
                         title_size = 18, title_size_i = 18, title_position = .5,
                         axis_text_size = 16, axis_text_size_i = 6, axis_text_size_bar = 20,
                         axis_title_size = 16, axis_title_size_i = 10, 
                         Countlimits = c(0, 4000), Countbreaks = seq(0,4000, by = 1000),
                         Arealimits = c(0, 100), Areabreaks = seq(0,100, by = 20),
                         Luminancelimits = c(0, 40), Luminancebreaks = seq(0,40, by = 10),
                         REDintensitylimits = c(0, 180), REDintensitybreaks = seq(0,180, by = 30),
                         Histmaxlimits_y = c(0, 0.08), Histmaxbreaks_y = seq(0.00, 0.08, by = 0.02),
                         Histmaxlimits_x = c(0, 100), Histmaxbreaks_x = seq(0,100, by = 20),
                         Histmeanlimits_y = c(0, 0.20), Histmeanbreaks_y = seq(0.00,0.20, by = 0.05),
                         Histmeanlimits_x = c(0, 50), Histmeanbreaks_x = seq(0,50, by = 10),
                         axis_text_color = "black", axis_line_color = "black",
                         SED = "SED", LIE = "LIE", MOE = "MOE", NCS = "NCS",
                         color_SED = "grey85", color_LIE = "skyblue", color_MOE = "lightgreen", color_NCS = "goldenrod",
                         jitter_fill_color = "black", jitter_fill_color_box = "white", jitter_surrounding_color = "black",
                         jitter_shape = 21, jitter_height = 0, jitter_height_all = 0, jitter_width = .1, jitter_width_all = .2,
                         jitter_size = 3, jitter_size_all = 1.5, jitter_size_boxplot = 2, jitter_alpha = .7, 
                         barplot_posision = 'dodge', barplot_surrounding_color = "black", barplot_width = .7, barplot_width_all = .7, 
                         boxplot_surrounding_color = "black", boxplot_width = .5, boxplot_symbol = "mean",
                         boxplot_symbol_shape = 23, boxplot_symbol_size = 4, boxplot_symbol_fill = "white",
                         Result_folder = "Result", Plot_folder = "Plot", Region_folder = "Region", All_folder = "All",
                         DV_folder = "DV", Individual_folder = "Individual_data", Hist_folder = "Histgram",
                         save_plot_width = 3, save_lineplot_width = 3.5, save_plot_height = 3, save_plot_dpi = 300,
                         save_plot_width_i = 12, save_plot_height_i = 9,save_plot_width_all = 9, save_histplot_width = 12, 
                         save_histplot_height = 9,
                         plot_lang = ""
){
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  
  #Times New Roman
  windowsFonts("TNR" = windowsFont("Times New Roman"))
  #Meiryo
  windowsFonts("MEI" = windowsFont("Meiryo"))
  
  if(plot_lang == "ja"){
    Font_family = "Meiryo"
    Plot_folder = paste0(Plot_folder, "_ja")
  }
  
  # Create Result Folder
  dir.create(paste0(Folder_name,"/", Result_folder), showWarnings = F)
  dir.create(paste0(Folder_name, "/", Result_folder, "/", Plot_folder), showWarnings = F)
  dir.create(paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", Region_folder), showWarnings = F)
  dir.create(paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", All_folder), showWarnings = F)
  dir.create(paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", DV_folder), showWarnings = F)
  dir.create(paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", Individual_folder), showWarnings = F)
  dir.create(paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", Hist_folder), showWarnings = F)
  
  
  #Import from AUTOPLOT
  #Data retrieval----
  namestmp <- list.files(path = Path,
                         full.names = F,
                         pattern = paste0("\\.", File_type, "$")) %>% 
    gsub(paste0(".", File_type), "",.) %>% 
    gsub("Macro", "",.)
  
  paths <- list.files(path = Path,
                      full.names = T,
                      pattern = paste0(".", File_type, "$"))
  #Modifying Data----
  dataset <- function(rawdata){
    name <- gsub(paste0("_", Region), "", namestmp[i])
    tmp <- read.csv(rawdata, skip = 3, fileEncoding = "Shift-JIS") %>%
      dplyr::slice(5) %>%
      # dplyr::select(2,4,8) %>% 
      dplyr::select(2) %>% 
      rename("Area" = "–ÊÏ") %>% 
      # rename("LI" = "‹P“x..ÏŽZ.") %>% 
      # rename("RI" = paste0(Aim_color, "..ÏŽZ.")) %>% 
      mutate(Area = as.numeric(Area)) %>% 
      # mutate(LI = as.numeric(LI)) %>% 
      # mutate(RI = as.numeric(RI)) %>% 
      mutate(No = name) %>% 
      mutate(Region = Region) %>% 
      mutate(Group = if_else(str_detect(No, pattern = SED),SED,
                             if_else(str_detect(No, pattern = LIE),LIE,
                             if_else(str_detect(No, pattern = MOE),MOE,NCS)))) %>% 
      # dplyr::select(4,6,5,1,2,3)
      dplyr::select(2,4,3,1)
  }
  countdata <- function(rawdata){
    name <- gsub(paste0("_", Region), "", namestmp[i])
    tmp <- read.csv(rawdata, nrows = 2, header = FALSE, fileEncoding = "Shift-JIS") %>% 
      dplyr::select(2) %>% 
      dplyr::slice(2) %>% 
      rename("Count" = "V2") %>% 
      mutate(Count = as.numeric(Count)) %>% 
      mutate(No = name) %>% 
      mutate(Region = Region) %>% 
      mutate(Group = if_else(str_detect(No, pattern = SED),SED,
                             if_else(str_detect(No, pattern = LIE),LIE,
                             if_else(str_detect(No, pattern = MOE),MOE,NCS)))) %>% 
      
      dplyr::select(2,3,4,1)
  }
  luminancedata <- function(rawdata){
    name <- gsub(paste0("_", Region), "", namestmp[i])
    tmp <- read.csv(rawdata, skip = 3, fileEncoding = "Shift-JIS") %>%
      dplyr::select(7, 11) %>% 
      dplyr::slice(1) %>% 
      rename("LI" = "‹P“x..•½‹Ï.") %>% 
      rename("RI" = paste0(Aim_color, "..•½‹Ï.")) %>% 
      mutate(LI = as.numeric(LI)) %>% 
      mutate(RI = as.numeric(RI)) %>% 
      mutate(No = name) %>% 
      mutate(Region = Region) %>% 
      mutate(Group = if_else(str_detect(No, pattern = SED),SED,
                             if_else(str_detect(No, pattern = LIE),LIE,
                             if_else(str_detect(No, pattern = MOE),MOE,NCS)))) %>% 
      dplyr::select(3,5,4,1,2)
  }
  luminancehistdata <- function(rawdata){
    name <- gsub(paste0("_", Region), "", namestmp[i])
    tmp <- read.csv(rawdata, skip = 10, fileEncoding = "Shift-JIS") %>% 
      dplyr::select(5, 7) %>% 
      rename("Luminance Max" = "‹P“x..Å‘å.") %>% 
      rename("Luminance Mean" = "‹P“x..•½‹Ï.") %>% 
      mutate("Luminance Max" = as.numeric(`Luminance Max`)) %>% 
      mutate("Luminance Mean" = as.numeric(`Luminance Mean`)) %>% 
      mutate(No = name) %>% 
      mutate(Region = Region) %>% 
      mutate(Group = if_else(str_detect(No, pattern = SED),SED,
                             if_else(str_detect(No, pattern = LIE),LIE,
                             if_else(str_detect(No, pattern = MOE),MOE,NCS)))) %>% 
      dplyr::select(3,4,5,1,2)
  }
  
  #Main Loop----
  RED_1 <- data.frame()
  RED_2 <- data.frame()
  RED_3 <- data.frame()
  RED_4 <- data.frame()
  RED_5 <- data.frame()
  RED_6 <- data.frame()
  RED_7 <- data.frame()
  RED_8 <- data.frame()
  RED_9 <- data.frame()
  CountRED_1 <- data.frame()
  CountRED_2 <- data.frame()
  CountRED_3 <- data.frame()
  CountRED_4 <- data.frame()
  CountRED_5 <- data.frame()
  CountRED_6 <- data.frame()
  CountRED_7 <- data.frame()
  CountRED_8 <- data.frame()
  CountRED_9 <- data.frame()
  LuminanceRED_1 <- data.frame()
  LuminanceRED_2 <- data.frame()
  LuminanceRED_3 <- data.frame()
  LuminanceRED_4 <- data.frame()
  LuminanceRED_5 <- data.frame()
  LuminanceRED_6 <- data.frame()
  LuminanceRED_7 <- data.frame()
  LuminanceRED_8 <- data.frame()
  LuminanceRED_9 <- data.frame()
  LuminancehistRED_1 <- data.frame()
  LuminancehistRED_2 <- data.frame()
  LuminancehistRED_3 <- data.frame()
  LuminancehistRED_4 <- data.frame()
  LuminancehistRED_5 <- data.frame()
  LuminancehistRED_6 <- data.frame()
  LuminancehistRED_7 <- data.frame()
  LuminancehistRED_8 <- data.frame()
  LuminancehistRED_9 <- data.frame()
  GREEN_1 <- data.frame()
  GREEN_2 <- data.frame()
  GREEN_3 <- data.frame()
  GREEN_4 <- data.frame()
  GREEN_5 <- data.frame()
  GREEN_6 <- data.frame()
  GREEN_7 <- data.frame()
  GREEN_8 <- data.frame()
  GREEN_9 <- data.frame()
  
  if(PFC == "F"){
    for(i in 1:length(paths)){
      for (Region in Region_list){
        if((Region == "dDGsp") || (Region == "dDGip") || (Region == "dCA1") || (Region == "dCA3")){
          if(str_detect(namestmp[i], pattern = "_R")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              RED_1 <- dplyr::bind_rows(RED_1, assign(namestmp[i], dataset(paths[i])))
              CountRED_1 <- dplyr::bind_rows(CountRED_1, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_1 <- dplyr::bind_rows(LuminanceRED_1, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_1 <- dplyr::bind_rows(LuminancehistRED_1, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              RED_2 <- dplyr::bind_rows(RED_2, assign(namestmp[i], dataset(paths[i])))
              CountRED_2 <- dplyr::bind_rows(CountRED_2, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_2 <- dplyr::bind_rows(LuminanceRED_2, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_2 <- dplyr::bind_rows(LuminancehistRED_2, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_3"))){
              RED_3 <- dplyr::bind_rows(RED_3, assign(namestmp[i], dataset(paths[i])))
              CountRED_3 <- dplyr::bind_rows(CountRED_3, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_3 <- dplyr::bind_rows(LuminanceRED_3, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_3 <- dplyr::bind_rows(LuminancehistRED_3, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_4"))){
              RED_4 <- dplyr::bind_rows(RED_4, assign(namestmp[i], dataset(paths[i])))
              CountRED_4 <- dplyr::bind_rows(CountRED_4, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_4 <- dplyr::bind_rows(LuminanceRED_4, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_4 <- dplyr::bind_rows(LuminancehistRED_4, assign(namestmp[i], luminancehistdata(paths[i])))
            }
          }
          else if(str_detect(namestmp[i], pattern = "_G")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              GREEN_1 <- dplyr::bind_rows(GREEN_1, assign(namestmp[i], dataset(paths[i])))
              #CountGREEN_1 <- dplyr::bind_rows(CountGREEN_1, assign(namestmp[i], countdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              GREEN_2 <- dplyr::bind_rows(GREEN_2, assign(namestmp[i], dataset(paths[i])))
              #CountGREEN_1 <- dplyr::bind_rows(CountGREEN_2, assign(namestmp[i], countdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_3"))){
              GREEN_3 <- dplyr::bind_rows(GREEN_3, assign(namestmp[i], dataset(paths[i])))
              #CountGREEN_1 <- dplyr::bind_rows(CountGREEN_2, assign(namestmp[i], countdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_4"))){
              GREEN_4 <- dplyr::bind_rows(GREEN_4, assign(namestmp[i], dataset(paths[i])))
              #CountGREEN_1 <- dplyr::bind_rows(CountGREEN_2, assign(namestmp[i], countdata(paths[i])))
            }
          }
          else if(str_detect(namestmp[i], pattern = "_B")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              BLUE_1 <- dplyr::bind_rows(BLUE_1, assign(namestmp[i], dataset(paths[i])))
              #CountBLUE_1 <- dplyr::bind_rows(CountBLUE_1, assign(namestmp[i], countdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              BLUE_2 <- dplyr::bind_rows(BLUE_2, assign(namestmp[i], dataset(paths[i])))
              #CountBLUE_2 <- dplyr::bind_rows(CountBLUE_2, assign(namestmp[i], countdata(paths[i])))
            }
            
          }
        }
        else{
          if(str_detect(namestmp[i], pattern = "_R")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              RED_5 <- dplyr::bind_rows(RED_5, assign(namestmp[i], dataset(paths[i])))
              CountRED_5 <- dplyr::bind_rows(CountRED_5, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_5 <- dplyr::bind_rows(LuminanceRED_5, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_5 <- dplyr::bind_rows(LuminancehistRED_5, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              RED_6 <- dplyr::bind_rows(RED_6, assign(namestmp[i], dataset(paths[i])))
              CountRED_6 <- dplyr::bind_rows(CountRED_6, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_6 <- dplyr::bind_rows(LuminanceRED_6, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_6 <- dplyr::bind_rows(LuminancehistRED_6, assign(namestmp[i], luminancehistdata(paths[i])))
            }
          }
          else if(str_detect(namestmp[i], pattern = "_G")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              GREEN_5 <- dplyr::bind_rows(GREEN_5, assign(namestmp[i], dataset(paths[i])))
              #CountGREEN_1 <- dplyr::bind_rows(CountGREEN_1, assign(namestmp[i], countdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              GREEN_6 <- dplyr::bind_rows(GREEN_6, assign(namestmp[i], dataset(paths[i])))
              #CountGREEN_1 <- dplyr::bind_rows(CountGREEN_2, assign(namestmp[i], countdata(paths[i])))
            }
          }
          else if(str_detect(namestmp[i], pattern = "_B")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              BLUE_1 <- dplyr::bind_rows(BLUE_1, assign(namestmp[i], dataset(paths[i])))
              #CountBLUE_1 <- dplyr::bind_rows(CountBLUE_1, assign(namestmp[i], countdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              BLUE_2 <- dplyr::bind_rows(BLUE_2, assign(namestmp[i], dataset(paths[i])))
              #CountBLUE_2 <- dplyr::bind_rows(CountBLUE_2, assign(namestmp[i], countdata(paths[i])))
            }
          }
        }
      }
    }
    
    assign(paste0(Red_name, "_1"), RED_1, envir = .GlobalEnv)
    assign(paste0(Red_name, "_2"), RED_2, envir = .GlobalEnv)
    assign(paste0(Red_name, "_3"), RED_3, envir = .GlobalEnv)
    assign(paste0(Red_name, "_4"), RED_4, envir = .GlobalEnv)
    assign(paste0(Red_name, "_5"), RED_5, envir = .GlobalEnv)
    assign(paste0(Red_name, "_6"), RED_6, envir = .GlobalEnv)
    assign(paste0("Count_", Red_name, "_1"), CountRED_1, envir = .GlobalEnv)
    assign(paste0("Count_", Red_name, "_2"), CountRED_2, envir = .GlobalEnv)
    assign(paste0("Count_", Red_name, "_3"), CountRED_3, envir = .GlobalEnv)
    assign(paste0("Count_", Red_name, "_4"), CountRED_4, envir = .GlobalEnv)
    assign(paste0("Count_", Red_name, "_5"), CountRED_5, envir = .GlobalEnv)
    assign(paste0("Count_", Red_name, "_6"), CountRED_6, envir = .GlobalEnv)
    assign(paste0("Luminance_", Red_name, "_1"), LuminanceRED_1, envir = .GlobalEnv)
    assign(paste0("Luminance_", Red_name, "_2"), LuminanceRED_2, envir = .GlobalEnv)
    assign(paste0("Luminance_", Red_name, "_3"), LuminanceRED_3, envir = .GlobalEnv)
    assign(paste0("Luminance_", Red_name, "_4"), LuminanceRED_4, envir = .GlobalEnv)
    assign(paste0("Luminance_", Red_name, "_5"), LuminanceRED_5, envir = .GlobalEnv)
    assign(paste0("Luminance_", Red_name, "_6"), LuminanceRED_6, envir = .GlobalEnv)
    assign(paste0(Green_name, "_1"), GREEN_1, envir = .GlobalEnv)
    assign(paste0(Green_name, "_2"), GREEN_2, envir = .GlobalEnv)
    assign(paste0(Green_name, "_3"), GREEN_3, envir = .GlobalEnv)
    assign(paste0(Green_name, "_4"), GREEN_4, envir = .GlobalEnv)
    assign(paste0(Green_name, "_5"), GREEN_5, envir = .GlobalEnv)
    assign(paste0(Green_name, "_6"), GREEN_6, envir = .GlobalEnv)
    
    GatheringData_Dorsal <- function(Data1, Data2, Data3, Data4, Data5, Data6, Data7, Data8, 
                                     CountData1, CountData2, CountData3, CountData4, 
                                     luminancedata1, luminancedata2, luminancedata3, luminancedata4){
      
      Data <-  dplyr::bind_cols(Data1, Data2, Data3, Data4, Data5, Data6, Data7, Data8, 
                                CountData1, CountData2, CountData3, CountData4, 
                                luminancedata1, luminancedata2, luminancedata3, luminancedata4) %>% 
        # mutate(Area_R = (Area...4 + Area...10)) %>% 
        # mutate(Area_G = (Area...16 + Area...22)) %>% 
        # mutate(Count_R = (Count...28 + Count...32)) %>% 
        # mutate(LI_R = (LI...5 + LI...11)) %>%
        # mutate(RI_R = (RI...6 + RI...12)) %>%
        mutate(Area_R = (Area...4 + Area...8 + Area...12 + Area...16)) %>% 
        mutate(Area_G = (Area...20 + Area...24 + Area...28 + Area...32)) %>% 
        mutate(Count_R = (Count...36 + Count...40 + Count...44 + Count...48)) %>% 
        mutate(LI_R = (LI...52 + LI...57 + LI...62 + LI...67)) %>%
        mutate(RI_R = (RI...53 + RI...58 + RI...63 + RI...68)) %>%
        mutate(LI = (LI_R/4)) %>% 
        mutate(RI = (RI_R/4)) %>% 
        mutate(Area = (Area_R/Area_G)*100) %>%
        mutate(Count = (Count_R/Area_G)) %>% 
        # mutate(LI = (LI_R/Area_G)) %>% 
        # mutate(RI = (RI_R/Area_G)) %>%
        rename("No" = "No...1") %>% 
        rename("Group" = "Group...2") %>% 
        rename("Region" = "Region...3") %>% 
        mutate(Group = as.factor(Group)) %>%
        mutate(Region = as.factor(Region)) %>%
        mutate(No = as.factor(No)) %>% 
        mutate(Group = relevel(Group, ref = graphic_refference)) %>% 
        dplyr::select("No", "Group", "Region", 
                      "Area_R", "Area_G", "Count_R", 
                      "LI_R", "RI_R", "LI", "RI", "Area", "Count")
    }
    
    GatheringData_Ventral <- function(Data1, Data2, Data3, Data4, CountData1, CountData2, luminancedata1, luminancedata2){
      
      Data <-  dplyr::bind_cols(Data1, Data2, Data3, Data4, CountData1, CountData2, luminancedata1, luminancedata2) %>% 
        # mutate(Area_R = (Area...4 + Area...10)) %>% 
        # mutate(Area_G = (Area...16 + Area...22)) %>% 
        # mutate(Count_R = (Count...28 + Count...32)) %>% 
        # mutate(LI_R = (LI...5 + LI...11)) %>%
        # mutate(RI_R = (RI...6 + RI...12)) %>%
        mutate(Area_R = (Area...4 + Area...8)) %>% 
        mutate(Area_G = (Area...12 + Area...16)) %>% 
        mutate(Count_R = (Count...20 + Count...24)) %>% 
        mutate(LI_R = (LI...28 + LI...33)) %>%
        mutate(RI_R = (RI...29 + RI...34)) %>%
        mutate(LI = (LI_R/2)) %>% 
        mutate(RI = (RI_R/2)) %>% 
        mutate(Area = (Area_R/Area_G)*100) %>%
        mutate(Count = (Count_R/Area_G)) %>% 
        # mutate(LI = (LI_R/Area_G)) %>% 
        # mutate(RI = (RI_R/Area_G)) %>%
        rename("No" = "No...1") %>% 
        rename("Group" = "Group...2") %>% 
        rename("Region" = "Region...3") %>% 
        mutate(Group = as.factor(Group)) %>%
        mutate(Region = as.factor(Region)) %>%
        mutate(No = as.factor(No)) %>% 
        mutate(Group = relevel(Group, ref = graphic_refference)) %>% 
        dplyr::select("No", "Group", "Region", 
                      "Area_R", "Area_G", "Count_R", 
                      "LI_R", "RI_R", "LI", "RI", "Area", "Count")
    }
    
    if(Analyse_color == "R/G"){
      GatheringData_RG_d <- GatheringData_Dorsal(Data1 = RED_1, Data2 = RED_2, Data3 = RED_3, Data4 = RED_4,
                                                 Data5 = GREEN_1, Data6 = GREEN_2, Data7 = GREEN_3, Data8 = GREEN_4,
                                                 CountData1 =  CountRED_1, CountData2 = CountRED_2, CountData3 = CountRED_3, CountData4 = CountRED_4,
                                                 luminancedata1 = LuminanceRED_1, luminancedata2 = LuminanceRED_2, 
                                                 luminancedata3 = LuminanceRED_3, luminancedata4 = LuminanceRED_4)
      GatheringData_RG_v <- GatheringData_Ventral(Data1 = RED_5, Data2 = RED_6, Data3 = GREEN_5, Data4 = GREEN_6,
                                                  CountData1 = CountRED_5, CountData2 = CountRED_6,
                                                  luminancedata1 = LuminanceRED_5, luminancedata2 = LuminanceRED_6)
      GatheringData_RG <- bind_rows(GatheringData_RG_d, GatheringData_RG_v)
      HistData_RG <- bind_rows(LuminancehistRED_1, LuminancehistRED_2, LuminancehistRED_3, LuminancehistRED_4,
                               LuminancehistRED_5, LuminancehistRED_6) %>% 
        mutate(Group = as.factor(Group)) %>% 
        mutate(Group = relevel(Group, ref = graphic_refference))
      assign("GatheringData_Hipp", GatheringData_RG, envir = .GlobalEnv)
      assign("HistData_Hipp", HistData_RG, envir = .GlobalEnv)
      write.csv(GatheringData_RG, paste0(Folder_name, "/", Result_folder,"/GatheringData.csv"))
      write.csv(HistData_RG, paste0(Folder_name, "/", Result_folder,"/HistData.csv"))
    }
  }
  else{
    for (Region in Region_list) {
      if(Region == "IL"){
        for(i in 1:length(paths)){
          if(str_detect(namestmp[i], pattern = "_R")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              RED_1 <- dplyr::bind_rows(RED_1, assign(namestmp[i], dataset(paths[i])))
              CountRED_1 <- dplyr::bind_rows(CountRED_1, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_1 <- dplyr::bind_rows(LuminanceRED_1, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_1 <- dplyr::bind_rows(LuminancehistRED_1, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              RED_2 <- dplyr::bind_rows(RED_2, assign(namestmp[i], dataset(paths[i])))
              CountRED_2 <- dplyr::bind_rows(CountRED_2, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_2 <- dplyr::bind_rows(LuminanceRED_2, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_2 <- dplyr::bind_rows(LuminancehistRED_2, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_3"))){
              RED_3 <- dplyr::bind_rows(RED_3, assign(namestmp[i], dataset(paths[i])))
              CountRED_3 <- dplyr::bind_rows(CountRED_3, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_3 <- dplyr::bind_rows(LuminanceRED_3, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_3 <- dplyr::bind_rows(LuminancehistRED_3, assign(namestmp[i], luminancehistdata(paths[i])))
            }
          }
          else if(str_detect(namestmp[i], pattern = "_G")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              GREEN_1 <- dplyr::bind_rows(GREEN_1, assign(namestmp[i], dataset(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              GREEN_2 <- dplyr::bind_rows(GREEN_2, assign(namestmp[i], dataset(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_3"))){
              GREEN_3 <- dplyr::bind_rows(GREEN_3, assign(namestmp[i], dataset(paths[i])))
            }
          }
        }
        GatheringData <- function(Data1, Data2, Data3, Data4, Data5, Data6, 
                                  CountData1, CountData2, CountData3,
                                  LuminanceData1, LuminanceData2, LuminanceData3){
          
          Data <-  dplyr::bind_cols(Data1, Data2, Data3, Data4, Data5, Data6, 
                                    CountData1, CountData2, CountData3,
                                    LuminanceData1, LuminanceData2, LuminanceData3) %>% 
            mutate(Area_R = (Area...4 + Area...8 + Area...12)) %>% 
            mutate(Area_G = (Area...16 + Area...20 + Area...24)) %>% 
            mutate(Count_R = (Count...28 + Count...32 + Count...36)) %>% 
            mutate(LI_R = (LI...40 + LI...45 + LI...50)) %>% 
            mutate(RI_R = (RI...41 + RI...46 + RI...51)) %>% 
            mutate(LI = (LI_R/3)) %>% 
            mutate(RI = (RI_R/3)) %>%
            mutate(Area = (Area_R/Area_G)*100) %>%
            mutate(Count = (Count_R/Area_G)* 1000000) %>% 
            rename("No" = "No...1") %>% 
            rename("Group" = "Group...2") %>% 
            rename("Region" = "Region...3") %>% 
            mutate(Group = as.factor(Group)) %>%
            mutate(Region = as.factor(Region)) %>%
            mutate(No = as.factor(No)) %>% 
            mutate(Group = relevel(Group, ref = graphic_refference)) %>% 
            dplyr::select("No", "Group", "Region", 
                          "Area_R", "Area_G", "Count_R", 
                          "LI_R", "RI_R", "LI", "RI", "Area", "Count")
        }
        
        GatheringData_IL <- GatheringData(Data1 = RED_1, Data2 = RED_2, Data3 = RED_3,
                                          Data4 = GREEN_1, Data5 = GREEN_2, Data6 = GREEN_3,  
                                          CountData1 =  CountRED_1, CountData2 = CountRED_2, CountData3 = CountRED_3,
                                          LuminanceData1 = LuminanceRED_1, LuminanceData2 = LuminanceRED_2, LuminanceData3 = LuminanceRED_3)
      }
      else if(Region == "PL"){
        for(i in 1:length(paths)){
          if(str_detect(namestmp[i], pattern = "_R")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              RED_4 <- dplyr::bind_rows(RED_4, assign(namestmp[i], dataset(paths[i])))
              CountRED_4 <- dplyr::bind_rows(CountRED_4, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_4 <- dplyr::bind_rows(LuminanceRED_4, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_4 <- dplyr::bind_rows(LuminancehistRED_4, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              RED_5 <- dplyr::bind_rows(RED_5, assign(namestmp[i], dataset(paths[i])))
              CountRED_5 <- dplyr::bind_rows(CountRED_5, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_5 <- dplyr::bind_rows(LuminanceRED_5, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_5 <- dplyr::bind_rows(LuminancehistRED_5, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_3"))){
              RED_6 <- dplyr::bind_rows(RED_6, assign(namestmp[i], dataset(paths[i])))
              CountRED_6 <- dplyr::bind_rows(CountRED_6, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_6 <- dplyr::bind_rows(LuminanceRED_6, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_6 <- dplyr::bind_rows(LuminancehistRED_6, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_4"))){
              RED_7 <- dplyr::bind_rows(RED_7, assign(namestmp[i], dataset(paths[i])))
              CountRED_7 <- dplyr::bind_rows(CountRED_7, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_7 <- dplyr::bind_rows(LuminanceRED_7, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_7 <- dplyr::bind_rows(LuminancehistRED_7, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_5"))){
              RED_8 <- dplyr::bind_rows(RED_8, assign(namestmp[i], dataset(paths[i])))
              CountRED_8 <- dplyr::bind_rows(CountRED_8, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_8 <- dplyr::bind_rows(LuminanceRED_8, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_8 <- dplyr::bind_rows(LuminancehistRED_8, assign(namestmp[i], luminancehistdata(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_6"))){
              RED_9 <- dplyr::bind_rows(RED_9, assign(namestmp[i], dataset(paths[i])))
              CountRED_9 <- dplyr::bind_rows(CountRED_9, assign(namestmp[i], countdata(paths[i])))
              LuminanceRED_9 <- dplyr::bind_rows(LuminanceRED_9, assign(namestmp[i], luminancedata(paths[i])))
              LuminancehistRED_9 <- dplyr::bind_rows(LuminancehistRED_9, assign(namestmp[i], luminancehistdata(paths[i])))
            }
          }
          else if(str_detect(namestmp[i], pattern = "_G")){
            if(str_detect(namestmp[i], pattern = paste0(Region, "_1"))){
              GREEN_4 <- dplyr::bind_rows(GREEN_4, assign(namestmp[i], dataset(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_2"))){
              GREEN_5 <- dplyr::bind_rows(GREEN_5, assign(namestmp[i], dataset(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_3"))){
              GREEN_6 <- dplyr::bind_rows(GREEN_6, assign(namestmp[i], dataset(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_4"))){
              GREEN_7 <- dplyr::bind_rows(GREEN_7, assign(namestmp[i], dataset(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_5"))){
              GREEN_8 <- dplyr::bind_rows(GREEN_8, assign(namestmp[i], dataset(paths[i])))
            }
            else if(str_detect(namestmp[i], pattern = paste0(Region, "_6"))){
              GREEN_9 <- dplyr::bind_rows(GREEN_9, assign(namestmp[i], dataset(paths[i])))
            }
          }
        }
        GatheringData <- function(Data1, Data2, Data3, Data4, Data5, Data6,
                                  Data7, Data8, Data9, Data10, Data11, Data12,
                                  CountData1, CountData2, CountData3, CountData4, CountData5, CountData6,
                                  LuminanceData1, LuminanceData2, LuminanceData3, LuminanceData4, LuminanceData5, LuminanceData6){
          
          Data <-  dplyr::bind_cols(Data1, Data2, Data3, Data4, Data5, Data6,
                                    Data7, Data8, Data9, Data10, Data11, Data12,
                                    CountData1, CountData2, CountData3, CountData4, CountData5, CountData6,
                                    LuminanceData1, LuminanceData2, LuminanceData3, LuminanceData4, LuminanceData5, LuminanceData6) %>% 
            mutate(Area_R = (Area...4 + Area...8 + Area...12 + Area...16 + Area...20 + Area...24)) %>% 
            mutate(Area_G = (Area...28 + Area...32 + Area...36 + Area...40 + Area...44 + Area...48)) %>% 
            mutate(Count_R = (Count...52 + Count...56 + Count...60 + Count...64 + Count...68 + Count...72)) %>% 
            mutate(LI_R = (LI...76 + LI...81 + LI...86 + LI...91 + LI...96 + LI...101)) %>% 
            mutate(RI_R = (RI...77 + RI...82 + RI...87 + RI...92 + RI...97 + RI...102)) %>% 
            mutate(LI = (LI_R/6)) %>% 
            mutate(RI = (RI_R/6)) %>%
            mutate(Area = (Area_R/Area_G)*100) %>%
            mutate(Count = (Count_R/Area_G)* 1000000) %>% 
            rename("No" = "No...1") %>% 
            rename("Group" = "Group...2") %>% 
            rename("Region" = "Region...3") %>% 
            mutate(Group = as.factor(Group)) %>%
            mutate(Region = as.factor(Region)) %>%
            mutate(No = as.factor(No)) %>% 
            mutate(Group = relevel(Group, ref = graphic_refference)) %>% 
            dplyr::select("No", "Group", "Region", 
                          "Area_R", "Area_G", "Count_R", 
                          "LI_R", "RI_R", "LI", "RI", "Area", "Count")
        }
        
        GatheringData_PL <- GatheringData(Data1 = RED_4, Data2 = RED_5, Data3 = RED_6, Data4 = RED_7, Data5 = RED_8, Data6 = RED_9, 
                                          Data7 = GREEN_4, Data8 = GREEN_5, Data9 = GREEN_6, Data10 = GREEN_7, Data11 = GREEN_8, Data12 = GREEN_9,
                                          CountData1 =  CountRED_4, CountData2 = CountRED_5, CountData3 = CountRED_6, 
                                          CountData4 =  CountRED_7, CountData5 = CountRED_8, CountData6 = CountRED_9,
                                          LuminanceData1 = LuminanceRED_4, LuminanceData2 = LuminanceRED_5, 
                                          LuminanceData3 = LuminanceRED_6, LuminanceData4 = LuminanceRED_7,
                                          LuminanceData5 = LuminanceRED_8, LuminanceData6 = LuminanceRED_9)
        
      }
    }
    GatheringData_RG <- dplyr::bind_rows(GatheringData_IL, GatheringData_PL)
    HistData_RG <- bind_rows(LuminancehistRED_1, LuminancehistRED_2, LuminancehistRED_3, LuminancehistRED_4,
                             LuminancehistRED_5, LuminancehistRED_6, LuminancehistRED_7, LuminancehistRED_8, LuminancehistRED_9) %>% 
      mutate(Group = as.factor(Group)) %>% 
      mutate(Group, relevel(Group, ref = graphic_refference))
    assign("GatheringData_PFC", GatheringData_RG, envir = .GlobalEnv)
    assign("HistData_PFC", HistData_RG, envir = .GlobalEnv)
    write.csv(GatheringData_RG, paste0(Folder_name, "/", Result_folder,"/GatheringData.csv")) 
    write.csv(HistData_RG, paste0(Folder_name, "/", Result_folder,"/HistData.csv"))
  }
  
  #Plot ----
  if(Font_family == "Times New Roman"){
    Font = "TNR"
  }
  else if(Font_family == "Meiryo"){
    Font = "MEI"
  }
  theme_bar <- function(){
    theme_classic(base_family = Font) +
      theme(plot.title = element_text(size = title_size, hjust = title_position),
            legend.position = "none",
            axis.text.x = element_text(size = axis_text_size_bar, colour = axis_text_color),
            axis.text.y = element_text(size = axis_text_size, colour = axis_text_color),
            axis.line = element_line(colour = axis_line_color),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = axis_title_size)) 
  }
  theme_bar_all <- function(){
    theme_classic(base_family = Font) +
      theme(plot.title = element_text(size = title_size, hjust = title_position),
            legend.position = legend_position,
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = legend_text_size),
            axis.text.x = element_text(size = axis_text_size, colour = axis_text_color),
            axis.text.y = element_text(size = axis_text_size, colour = axis_text_color),
            axis.line = element_line(colour = axis_line_color),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = axis_title_size)) 
  }
  theme_bar_DV <- function(){
    theme_classic(base_family = Font) +
      theme(plot.title = element_text(size = title_size, hjust = title_position),
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = legend_text_size),
            axis.text.x = element_text(size = axis_text_size_bar, colour = axis_text_color),
            axis.text.y = element_text(size = axis_text_size, colour = axis_text_color),
            axis.line = element_line(colour = axis_line_color),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = axis_title_size)) 
  }
  theme_bar_i <- function(){
    
    theme_classic(base_family = Font) +
      theme(plot.title = element_text(size = title_size_i, hjust = title_position),
            legend.position = "none",
            axis.text.x = element_text(size = axis_text_size_i, colour = axis_text_color),
            axis.text.y = element_text(size = axis_title_size_i, colour = axis_text_color),
            axis.line = element_line(colour = axis_line_color),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = axis_title_size_i))
  }
  theme_hist <- function(){
    theme_classic(base_family = Font) +
      theme(plot.title = element_text(size = title_size, hjust = title_position),
            legend.position = legend_position,
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = legend_text_size),
            axis.text.x = element_text(size = axis_text_size_bar, colour = axis_text_color),
            axis.text.y = element_text(size = axis_text_size, colour = axis_text_color),
            axis.line = element_line(colour = axis_line_color),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = axis_title_size)) 
  }
  
  Barplot <- function(dataset, datajitter, purpose, region, titlename){
    if(purpose == Count){
      g <- ggplot(dataset, aes(x = Group, y = meanCount, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision,
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanCount - seCount,
                          ymax = meanCount + seCount),
                      width = errobar_width) +
        geom_jitter(data = datajitter, aes(x = Group, y = Count),
                    height = jitter_height, width = jitter_width, 
                    size = jitter_size, alpha = jitter_alpha,
                    fill = jitter_fill_color, color = jitter_surrounding_color, 
                    shape = jitter_shape) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "& NeuN"^{"+"} ~ "/ NeuN"^{"+"} ~ "(# /" ~ mm^2 ~")"))) +
        scale_y_continuous(expand = c(0, 0), limits = Countlimits, breaks = Countbreaks) + 
        scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                     MOE = color_MOE, NCS = color_NCS)) +
        theme_bar()
    }
    else if(purpose == Area){
      g <- ggplot(dataset, aes(x = Group, y = meanArea, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanArea - seArea,
                          ymax = meanArea + seArea),
                      width = errobar_width) +
        geom_jitter(data = datajitter, aes(x = Group, y = Area),
                    height = jitter_height, width = jitter_width, 
                    size = jitter_size, alpha = jitter_alpha,
                    fill = jitter_fill_color, color = jitter_surrounding_color, 
                    shape = jitter_shape) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "area (% NeuN"^{"+"} ~ ")"))) +
        scale_y_continuous(expand = c(0, 0), limits = Arealimits, breaks = Areabreaks) + 
        scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                     MOE = color_MOE, NCS = color_NCS)) +
        theme_bar()
    }
    else if(purpose == Luminance){
      g <- ggplot(dataset, aes(x = Group, y = meanLI, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanLI - seLI,
                          ymax = meanLI + seLI),
                      width = errobar_width) +
        geom_jitter(data = datajitter, aes(x = Group, y = LI),
                    height = jitter_height, width = jitter_width, 
                    size = jitter_size, alpha = jitter_alpha,
                    fill = jitter_fill_color, color = jitter_surrounding_color, 
                    shape = jitter_shape) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "Luminance\n(Region average)"))) +
        scale_y_continuous(expand = c(0, 0), limits = Luminancelimits, breaks = Luminancebreaks) +
        scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                     MOE = color_MOE, NCS = color_NCS)) +
        theme_bar()
    }
    else if(purpose == Red_intensity){
      g <- ggplot(dataset, aes(x = Group, y = meanRI, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanRI - seRI,
                          ymax = meanRI + seRI),
                      width = errobar_width) +
        geom_jitter(data = datajitter, aes(x = Group, y = RI),
                    height = jitter_height, width = jitter_width, 
                    size = jitter_size, alpha = jitter_alpha,
                    fill = jitter_fill_color, color = jitter_surrounding_color, 
                    shape = jitter_shape) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "intensity\n(Region average)"))) +
        scale_y_continuous(expand = c(0, 0), limits = REDintensitylimits, 
                           breaks = REDintensitybreaks) +
        scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                     MOE = color_MOE, NCS = color_NCS)) +
        theme_bar()}
    plot(g)
    
    assign(paste0("Plot_", region, "_", purpose), g, envir = .GlobalEnv)
    
    if(Folder_name == ""){
      ggsave(filename = paste0(Folder_name, Result_folder, "/", Plot_folder, "/",
                               Region_folder, "/", name, plot_number, ".png"),
             width = save_plot_width, height = save_plot_height, dpi = save_plot_dpi)
    }
    else{
      ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", 
                               Region_folder, "/", name, plot_number, ".png"),
             width = save_plot_width, height = save_plot_height, dpi = save_plot_dpi)
    }
  }
  BarplotAll <- function(dataset, datajitter, purpose, titlename){
    if(purpose == Count){
      g <- ggplot(dataset, aes(x = Region, y = meanCount, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanCount - seCount,
                          ymax = meanCount + seCount),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = Region, y = Count, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "& NeuN"^{"+"} ~ "/ NeuN"^{"+"} ~ "(# /" ~ mm^2 ~")"))) +
        scale_y_continuous(expand = c(0, 0), limits = Countlimits, 
                           breaks = Countbreaks) + 
        scale_x_discrete(limits = Region_list) +
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_all()
    }
    else if(purpose == Area){
      g <- ggplot(dataset, aes(x = Region, y = meanArea, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanArea - seArea,
                          ymax = meanArea + seArea),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = Region, y = Area, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "area (% NeuN"^{"+"} ~ ")"))) +
        scale_y_continuous(expand = c(0, 0), limits = Arealimits, 
                           breaks = Areabreaks) +
        scale_x_discrete(limits = Region_list) +
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_all()
    }
    else if(purpose == Luminance){
      g <- ggplot(dataset, aes(x = Region, y = meanLI, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanLI - seLI,
                          ymax = meanLI + seLI),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = Region, y = LI, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "Luminance\n(Region average)"))) +
        scale_y_continuous(expand = c(0, 0), limits = Luminancelimits, 
                           breaks = Luminancebreaks) +
        scale_x_discrete(limits = Region_list) +
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_all()
    }
    else if(purpose == Red_intensity){
      g <- ggplot(dataset, aes(x = Region, y = meanRI, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanRI - seRI,
                          ymax = meanRI + seRI),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = Region, y = RI, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "intensity\n(Region average)"))) +
        scale_y_continuous(expand = c(0, 0), limits = REDintensitylimits, 
                           breaks = REDintensitybreaks) +
        scale_x_discrete(limits = Region_list) +
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_all()
    }
    plot(g)
    
    assign(paste0("Plot_all", purpose), g, envir = .GlobalEnv)
    
    if(Folder_name == ""){
      ggsave(filename = paste0(Folder_name, Result_folder, "/", Plot_folder, "/", 
                               All_folder, "/", name, plot_number, ".png"),
             width = save_plot_width_all, height = save_plot_height, dpi = save_plot_dpi)
    }
    else{
      ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", 
                               All_folder, "/", name, plot_number, ".png"),
             width = save_plot_width_all, height = save_plot_height, dpi = save_plot_dpi)
    }
  }
  BarplotD_V <- function(dataset, datajitter, purpose, titlename){
    if(purpose == Count){
      g <- ggplot(dataset, aes(x = D_V, y = meanCount, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanCount - seCount,
                          ymax = meanCount + seCount),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = D_V, y = mCount, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "& NeuN"^{"+"} ~ "/ NeuN"^{"+"} ~ "(# /" ~ mm^2 ~")"))) +
        scale_y_continuous(expand = c(0, 0), limits = Countlimits, 
                           breaks = Countbreaks) + 
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_DV()
    }
    else if(purpose == Area){
      g <- ggplot(dataset, aes(x = D_V, y = meanArea, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanArea - seArea,
                          ymax = meanArea + seArea),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = D_V, y = mArea, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "area (% NeuN"^{"+"} ~ ")"))) +
        scale_y_continuous(expand = c(0, 0), limits = Arealimits, 
                           breaks = Areabreaks) +
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_DV()
    }
    else if(purpose == Luminance){
      g <- ggplot(dataset, aes(x = D_V, y = meanLI, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision,
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanLI - seLI,
                          ymax = meanLI + seLI),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = D_V, y = mLI, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "Luminance\n(Region average)"))) +
        scale_y_continuous(expand = c(0, 0), limits = Luminancelimits, 
                           breaks = Luminancebreaks) +
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_DV()
    }
    else if(purpose == Red_intensity){
      g <- ggplot(dataset, aes(x = D_V, y = meanRI, fill = Group)) +
        geom_bar(stat = 'identity', position = barplot_posision, 
                 width = barplot_width, colour = barplot_surrounding_color) + 
        geom_errorbar(aes(ymin = meanRI - seRI,
                          ymax = meanRI + seRI),
                      width = errobar_width, position = position_dodge(barplot_width_all)) +
        geom_jitter(data = datajitter, aes(x = D_V, y = mRI, color = Group),
                    size = jitter_size_all, alpha = jitter_alpha, 
                    fill = jitter_fill_color, shape = jitter_shape,
                    position = position_jitterdodge(jitter.width = jitter_width_all, 
                                                    jitter.height = jitter_height_all)) +
        labs(title = titlename, y = bquote(paste("c-Fos"^{"+"} ~ "intensity\n(Region average)"))) +
        scale_y_continuous(expand = c(0, 0), limits = REDintensitylimits, 
                           breaks = REDintensitybreaks) +
        scale_color_manual(labels = c(SED = SED, LIE = LIE),
                           values = c(SED = jitter_surrounding_color, 
                                      LIE = jitter_surrounding_color)) +
        scale_fill_manual(labels = c(SED = SED, LIE = LIE),
                          values = c(SED = color_SED, LIE = color_LIE)) +
        theme_bar_DV()
    }
    plot(g)
    
    assign(paste0("Plot_DV", purpose), g, envir = .GlobalEnv)
    
    if(Folder_name == ""){
      ggsave(filename = paste0(Folder_name, Result_folder, "/", Plot_folder, "/", 
                               DV_folder, "/", name, plot_number, ".png"),
             width = save_lineplot_width, height = save_plot_height, dpi = save_plot_dpi)
    }
    else{
      ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", 
                               DV_folder, "/", name, plot_number, ".png"),
             width = save_lineplot_width, height = save_plot_height, dpi = save_plot_dpi)
    }
  }
  BarplotIndividual <- function(dataset, purpose, titlename){
    if(purpose == Count){
      assign(paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
             ggplot(dataset, aes(x = Region, y = Count, fill = Group)) +
               geom_bar(stat = 'identity', position = barplot_posision,
                        width = barplot_width, colour = barplot_surrounding_color) + 
               labs(title = paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
                    y = bquote(paste("c-Fos"^{"+"} ~ "& NeuN"^{"+"} ~ "/ NeuN"^{"+"} ~ "(# /" ~ mm^2 ~")")))+
               scale_y_continuous(expand = c(0, 0), limits = Countlimits, 
                                  breaks = Countbreaks) +
               scale_x_discrete(limits = Region_list) +
               scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                            MOE = color_MOE, NCS = color_NCS)) +
               theme_bar_i())
    }
    else if(purpose == Area){
      
      #g <-
      assign(paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
             ggplot(dataset, aes(x = Region, y = Area, fill = Group)) +
               geom_bar(stat = 'identity', position = barplot_posision, 
                        width = barplot_width, colour = barplot_surrounding_color) + 
               labs(title = paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
                    y = bquote(paste("c-Fos"^{"+"} ~ "area (% NeuN"^{"+"} ~ ")"))) +
               scale_y_continuous(expand = c(0, 0), limits = Arealimits, 
                                  breaks = Areabreaks) +
               scale_x_discrete(limits = Region_list) +
               scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                            MOE = color_MOE, NCS = color_NCS)) +
               theme_bar_i())
    }
    else if(purpose == Luminance){
      # g <- 
      assign(paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
             ggplot(dataset, aes(x = Region, y = LI, fill = Group)) +
               geom_bar(stat = 'identity', position = barplot_posision, 
                        width = barplot_width, colour = barplot_surrounding_color) + 
               labs(title = paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
                    y = bquote(paste("c-Fos"^{"+"} ~ "Luminance\n(Region average)"))) +
               scale_y_continuous(expand = c(0, 0), limits = Luminancelimits, 
                                  breaks = Luminancebreaks) +
               scale_x_discrete(limits = Region_list) +
               scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                            MOE = color_MOE, NCS = color_NCS)) +
               theme_bar_i())
    }
    else if(purpose == Red_intensity){
      # g <- 
      assign(paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
             ggplot(dataset, aes(x = Region, y = RI, fill = Group)) +
               geom_bar(stat = 'identity', position = barplot_posision, 
                        width = barplot_width, colour = barplot_surrounding_color) + 
               labs(title = paste0(gsub("_1_R", "", dataset$No)," ", titlename), 
                    y = bquote(paste("c-Fos"^{"+"} ~ "intensity\n(Region average)"))) +
               scale_y_continuous(expand = c(0, 0), limits = REDintensitylimits, 
                                  breaks = REDintensitybreaks) +
               scale_x_discrete(limits = Region_list) +
               scale_fill_manual(values = c(SED = color_SED, LIE = color_LIE, 
                                            MOE = color_MOE, NCS = color_NCS)) +
               theme_bar_i())  
    }
  }
  HistPlot <- function(dataset, titlename, luminance, region){
    if(luminance == "Luminance Max"){
      g <- ggplot(dataset, aes(x = `Luminance Max`, y = ..density.., 
                               fill = Group, colour = Group)) +
        geom_histogram(position = "identity", alpha = .4) +
        geom_density(alpha = 0) +
        labs(title = paste0(region, " ", luminance)) +
        scale_y_continuous(expand = c(0, 0), limits = Histmaxlimits_y, 
                           breaks = Histmaxbreaks_y) +
        scale_x_continuous(limits = Histmaxlimits_x, 
                           breaks = Histmaxbreaks_x) +
        scale_fill_manual(values = c(SED = color_SED, 
                                     LIE = color_LIE)) +
        scale_color_manual(values = c(SED = color_SED, 
                                      LIE = color_LIE)) +
        theme_hist()
    }
    else if(luminance == "Luminance Mean"){
      g <- ggplot(dataset, aes(x = `Luminance Mean`, y = ..density.., 
                               fill = Group, colour = Group)) +
        geom_histogram(position = "identity", alpha = .4) +
        geom_density(alpha = 0) +
        labs(title = paste0(region, " ", luminance)) +
        scale_y_continuous(expand = c(0, 0), limits = Histmeanlimits_y, 
                           breaks = Histmeanbreaks_y) +
        scale_x_continuous(limits = Histmeanlimits_x, 
                           breaks = Histmeanbreaks_x) +
        scale_fill_manual(values = c(SED = color_SED, 
                                     LIE = color_LIE)) +
        scale_color_manual(values = c(SED = color_SED, 
                                      LIE = color_LIE)) +
        theme_hist()
    }
    assign(paste0("Plot_", region, "_", luminance), g, envir = .GlobalEnv)
    # if(Folder_name == ""){
    #   ggsave(filename = paste0(Folder_name, Result_folder, "/", Plot_folder, "/", Hist_folder, "/", name, plot_number, ".png"),
    #          width = save_plot_width * 1.5, height = save_plot_height, dpi = save_plot_dpi)
    # }
    # else{
    #   ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", Plot_folder, "/", Hist_folder, "/", name, plot_number, ".png"),
    #          width = save_plot_width * 1.5, height = save_plot_height, dpi = save_plot_dpi)
    # }
  }
  
  DataList <- c()
  AllData <- data.frame()
  LuminanceMax <- list()
  LuminanceMean <- list()
  purpose <- c(Count, Area, Luminance, Red_intensity)
  for (REGION in Region_list) {
    assign(paste0("Data_", REGION), GatheringData_RG[GatheringData_RG$Region == REGION,] %>% 
             # dplyr::select(1,2,3,9,10,11,12) , envir = .GlobalEnv)
             dplyr::select(1,2,3,7,8,9,10,11,12) , envir = .GlobalEnv)
    
    Sumdata <- GatheringData_RG[GatheringData_RG$Region == REGION,] %>% 
      group_by(Group, Region) %>% 
      summarise(meanCount = mean(Count),
                meanArea = mean(Area),
                meanLI = mean(LI),
                meanRI = mean(RI),
                seCount = sd(Count)/sqrt(n()-1),
                seArea = sd(Area)/sqrt(n()-1),
                seLI = sd(LI)/sqrt(n()-1),
                seRI = sd(RI)/sqrt(n()-1)) %>% 
      mutate(Group = as.factor(Group)) %>%
      mutate(Region = as.factor(Region)) %>%
      mutate(Group = relevel(Group, ref = graphic_refference))
    
    AllData <- dplyr::bind_rows(AllData, Sumdata)
    
    for (Purpose in purpose) {
      name <- paste0(REGION, "_", Purpose)
      Barplot(dataset = Sumdata, datajitter = GatheringData_RG[GatheringData_RG$Region == REGION,], 
              purpose = Purpose, region = REGION, titlename = Purpose)
    }
    luminance <- c("Luminance Max", "Luminance Mean")
    for (LuminanceHist in luminance) {
      name <- paste0(REGION, "_", LuminanceHist)
      HistPlot(dataset = HistData_RG[HistData_RG$Region == REGION,],
               titlename = LuminaceHist, luminance = LuminanceHist, region = REGION)
      if(LuminanceHist == "Luminance Max"){
        LuminanceMax <- c(LuminanceMax, 
                          list(HistPlot(dataset = HistData_RG[HistData_RG$Region == REGION,],
                                        titlename = LuminaceHist, luminance = LuminanceHist, region = REGION)))
      }
      else if(LuminanceHist == "Luminance Mean"){
        LuminanceMean <- c(LuminanceMean, 
                           list(HistPlot(dataset = HistData_RG[HistData_RG$Region == REGION,],
                                         titlename = LuminaceHist, luminance = LuminanceHist, region = REGION)))
      }
    }
    wrap_plots(LuminanceMax)
    ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", 
                             Plot_folder, "/", Hist_folder, "/Luminance_Max", plot_number, ".png"),
           width = save_histplot_width, height = save_histplot_height, 
           dpi = save_plot_dpi, device = "png")
    wrap_plots(LuminanceMean)
    ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", 
                             Plot_folder, "/", Hist_folder, "/Luminance_Mean", plot_number, ".png"),
           width = save_histplot_width, height = save_histplot_height, 
           dpi = save_plot_dpi, device = "png")
  }
  
  for (Purpose in purpose) {
    name <- paste0("All_", Purpose)
    BarplotAll(dataset = AllData, datajitter = GatheringData_RG, 
               purpose = Purpose, titlename = Purpose)
  }
  
  if(D_V == "T"){
    D_Vdata <- GatheringData_RG %>% 
      mutate(ID = rep(1:number_of_sample, each = number_of_region)) %>% 
      mutate(D_V = if_else(str_detect(Region, pattern = "d"), "Dorsal", "Ventral")) %>% 
      mutate(D_V = as.factor(D_V))
    
    dvdata <- data.frame()
    DV <- c("Dorsal", "Ventral")
    for (id in 1:number_of_sample) {
      for (dv in DV) {
        d <- D_Vdata[D_Vdata$D_V == dv & D_Vdata$ID == id, ] %>% 
          group_by(No, Group, D_V) %>% 
          summarise(mCount = (sum(Count_R)/sum(Area_G))*1000000,
                    mArea = (sum(Area_R)/sum(Area_G))*100,
                    #mLI = sum(LI_R)/sum(Area_G),
                    mLI = sum(LI)/(length(LI)),
                    # mRI = sum(RI_R)/sum(Area_G)) %>% 
                    mRI = sum(RI)/(length(RI))) %>%
          mutate(Group = as.factor(Group)) %>%
          mutate(Group = relevel(Group, ref = graphic_refference))
        dvdata <- dplyr::bind_rows(dvdata, d)
      }
    }
    
    assign("DVdata", dvdata, envir = .GlobalEnv)
    write.csv(dvdata, paste0(Folder_name, "/", Result_folder,"/DV_Data.csv"))
    
    Sumdata <- dvdata %>% 
      group_by(Group, D_V) %>% 
      summarise(meanCount = mean(mCount),
                meanArea = mean(mArea),
                meanLI = mean(mLI),
                meanRI = mean(mRI),
                seCount = sd(mCount)/sqrt(n()-1),
                seArea = sd(mArea)/sqrt(n()-1),
                seLI = sd(mLI)/sqrt(n()-1),
                seRI = sd(mRI)/sqrt(n()-1)) %>% 
      mutate(Group = as.factor(Group)) %>%
      mutate(Group = relevel(Group, ref = graphic_refference))
    
    for (Purpose in purpose) {
      name <- paste0("DV_", Purpose)
      BarplotD_V(dataset = Sumdata, datajitter = dvdata, purpose = Purpose, 
                 titlename = Purpose)
    }
  }
  
  if(PFC == "F"){
    IData <- GatheringData_RG %>% 
      mutate(ID = rep(1:number_of_sample, each = number_of_region))
  }
  else{
    IData <- GatheringData_RG %>% 
      mutate(ID = rep(1:number_of_sample, times = number_of_region))
  }
  
  
  plotlistCount <- list()
  plotlistArea <- list()
  plotlistLuminance <- list()
  plotlistREDintensity <- list()
  for (p in purpose) {
    if(p == Count){
      for (id in 1:number_of_sample) {
        plotlistCount <- c(plotlistCount, 
                           list(BarplotIndividual(dataset = IData[IData$ID == id,], 
                                                  purpose = p, titlename = p)))
      }
      g <- wrap_plots(plotlistCount)
      plot(g)
      ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", 
                               Plot_folder, "/", Individual_folder, "/", Count, plot_number, ".png"),
             width = save_plot_width_i, height = save_plot_height_i, dpi = save_plot_dpi, device = "png")
    }
    else if(p == Area){
      for (id in 1:number_of_sample) {
        plotlistArea <- c(plotlistArea, 
                          list(BarplotIndividual(dataset = IData[IData$ID == id,], 
                                                 purpose = p, titlename = p)))
      }
      g <- wrap_plots(plotlistArea)
      plot(g)
      ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", 
                               Plot_folder, "/", Individual_folder, "/", Area, plot_number, ".png"),
             width = save_plot_width_i, height = save_plot_height_i, dpi = save_plot_dpi, device = "png")
    }
    else if(p == Luminance){
      for (id in 1:number_of_sample) {
        plotlistLuminance <- c(plotlistLuminance, 
                               list(BarplotIndividual(dataset = IData[IData$ID == id,],
                                                      purpose = p, titlename = p)))
      }
      g <- wrap_plots(plotlistLuminance)
      plot(g)
      ggsave(filename = paste0(Folder_name, "/", Result_folder, "/", 
                               Plot_folder, "/", Individual_folder, "/", Luminance, plot_number, ".png"),
             width = save_plot_width_i, height = save_plot_height_i, dpi = save_plot_dpi, device = "png")
    }
    else if(p == Red_intensity){
      for (id in 1:number_of_sample) {
        plotlistREDintensity <- c(plotlistREDintensity, 
                                  list(BarplotIndividual(dataset = IData[IData$ID == id,], 
                                                         purpose = p, titlename = p)))
      }
      g <- wrap_plots(plotlistREDintensity)
      plot(g)
      ggsave(filename = paste0(Folder_name, "/", Result_folder, "/",
                               Plot_folder, "/", Individual_folder, "/", Red_intensity, plot_number, ".png"),
             width = save_plot_width_i, height = save_plot_height_i, dpi = save_plot_dpi, device = "png")
    }
    assign("PlotListCount", plotlistCount, envir = .GlobalEnv)
    assign("PlotListArea", plotlistArea, envir = .GlobalEnv)
    assign("PlotListLuminace", plotlistLuminance, envir = .GlobalEnv)
    assign("PlotListREDintensity", plotlistREDintensity, envir = .GlobalEnv)
  }
  rm(g)
}
