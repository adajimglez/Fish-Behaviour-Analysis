# Prerequesites
# 1. Make sure all metadata follows a consistent formatting
# 2.Create a folder per experiment containing one sub folder for metadata and one for behaviour
# 3. Save all metadata and behaviour as .csv in their respective folders (making sure is pure csv and not csv UTF-8)

library(tidyverse)
library(ggforce)
library(ggpubr)
library(rstatix)
library(ggprism)
library(hgutils)
library(plotrix)

########################################## Define the working directories ##############################################
# Define your behaviour data directory

viewpoint_dir <- "G:/Shared drives/Project/Embryo_Swimming_Data/"

# Start by working from the metadata directory
setwd(viewpoint_dir)


# Read the necessary CSV files (we need to combine viewpoint and metadata)
viewpoint <- read.csv("1 second bins 20251205-172003.xls.csv")



######################################### Tidy the data #######################################

### 1. Create a new output folder and file name for this part of the analysis ### Name it with strain_drugs_repnumber
InputFolder<-viewpoint_dir
OutputFolder <- "AB_Zn" # This is where the replicates for the experiment will go. Label strain and compound used and other data if necessary. 
dir.create(paste0 (InputFolder, "/", OutputFolder)) #where the final data per replicate will go to be used in step 2

#### Set up a name for the final output file following a similar structure to the directory (strain_drug)
StrainDrug <- "AB_Zn" #label the name of strain and compound used. This will be used as label for file names later on
Replicate <- "R1"     #Keep the format of the naming consistent



## 2.  Adjust the analysis to your experimental design 

############# For mutant lines assign the wells to the genotype (in this example all A-D wells are wt) ##########


#### Indicate the position of each experimental condition using well number here (following viewpoint's "aname" column)
T1<- c("01|07")
T2 <- c("02|08")
T3 <- c("03|09")
T4 <- c("04|10")
T5 <- c("05|11")
T6 <- c("06|12")




#### Describe the names per group. This is well name-independent. e.g. group 1 will always be associated to the wells that you labelled as "control" in the previous step
group1<- "Control"
group2<- "200uM EDTA"
group3<- "1000uM EDTA"
group4<- "64uM Zn"
group5<- "64uM Zn + 200uM EDTA"
group6<- "64uM Zn + 1000uM EDTA"





## 3. Prepare the data

#Reassign your behaviour + metadata merged file to a new working variable for this part of the analysis
data <- viewpoint
head(data)


exclude <- c("B03", "C03", "H04", "06", "12", "B09", "G09", "H09") 
#Creates a vector with those wells that need removing if larvae show anatomical defects (e.g no inflation of swim bladder) or have been left empty. 
#Give specific well names here following the same format that viewpoint aname column. 

	


datafiltered<-data[grep(paste0(exclude, collapse = "|"), data$aname, invert = TRUE), ] # Creates a new dataframe that does not contain abnormal larvae



#### Calculate total distance and duration by adding small and large data

alldist <- datafiltered$smldist+datafiltered$lardist #Calculate the sum of small and large distance
alldur <-  datafiltered$smldur+datafiltered$lardur #Calculate the sum of small and large duration

datafiltered<- cbind(datafiltered, alldist) #Adds a new column to the dataframe with all distance
datafiltered<- cbind(datafiltered, alldur) #Adds a new column to the dataframe with all duration




#### Calculate the mean per larva of duration or distance keeping the separation between on/off stimuli

filteredDur<-aggregate(alldur~aname+stimuli_name+sttime+pn, data=datafiltered, FUN=mean) #This create a new dataframe that only contains the needed columns

filteredDur <- filteredDur %>% 
  filter(!(stimuli_name %in% c("on", "off")))


filteredDur <- filteredDur %>% 
  mutate(stimuli_name = case_when(
    pn >= 1 & pn <= 300 ~ "on",
    TRUE ~ stimuli_name
  ))

filteredDur <- filteredDur %>% 
  mutate(stimuli_name = case_when(
    pn >= 301 & pn <= 600 ~ "off",
    TRUE ~ stimuli_name
  ))


filteredDur <- filteredDur %>% 
  mutate(stimuli_name = case_when(
    pn >= 601 & pn <= 900 ~ "on",
    TRUE ~ stimuli_name
  ))


filteredDur <- filteredDur %>% 
  mutate(stimuli_name = case_when(
    pn >= 901 & pn <= 1200 ~ "off",
    TRUE ~ stimuli_name
  ))

filteredDur <- filteredDur %>% 
  mutate(stimuli_name = case_when(
    pn >= 1201 & pn <= 1500 ~ "on",
    TRUE ~ stimuli_name
  ))

filteredDur <- filteredDur %>% 
  mutate(stimuli_name = case_when(
    pn >= 1501 & pn <= 1800 ~ "off",
    TRUE ~ stimuli_name
  ))



filteredDist<-aggregate(alldist~aname+stimuli_name+sttime+pn, data=datafiltered, FUN=mean)

filteredDist <- filteredDist %>% 
  filter(!(stimuli_name %in% c("on", "off")))


filteredDist <- filteredDist %>% 
  mutate(stimuli_name = case_when(
    pn >= 1 & pn <= 300 ~ "on",
    TRUE ~ stimuli_name
  ))

filteredDist <- filteredDist %>% 
  mutate(stimuli_name = case_when(
    pn >= 301 & pn <= 600 ~ "off",
    TRUE ~ stimuli_name
  ))


filteredDist <- filteredDist %>% 
  mutate(stimuli_name = case_when(
    pn >= 601 & pn <= 900 ~ "on",
    TRUE ~ stimuli_name
  ))


filteredDist <- filteredDist %>% 
  mutate(stimuli_name = case_when(
    pn >= 901 & pn <= 1200 ~ "off",
    TRUE ~ stimuli_name
  ))

filteredDist <- filteredDist %>% 
  mutate(stimuli_name = case_when(
    pn >= 1201 & pn <= 1500 ~ "on",
    TRUE ~ stimuli_name
  ))

filteredDist <- filteredDist %>% 
  mutate(stimuli_name = case_when(
    pn >= 1501 & pn <= 1800 ~ "off",
    TRUE ~ stimuli_name
  ))



cbindDurDistMean <- cbind(aname = filteredDur$aname, time = filteredDur$sttime, bin = filteredDist$pn, stimuli = filteredDur$stimuli_name, alldur = filteredDur$alldur, alldist = filteredDist$alldist)
cbindDurDistMean <- as.data.frame(cbindDurDistMean) ### These two lines create a new dataframe that holds the average of the 3 reads per larva for both distance and duration 


cbindDurDistMean<-cbindDurDistMean %>% 
  mutate(Group = case_when(str_detect(aname, T1 ) ~ group1))
cbindDurDistMean<-cbindDurDistMean %>% 
  mutate(Group = case_when(str_detect(aname, T2 ) ~ group2,  TRUE ~ Group))
cbindDurDistMean<-cbindDurDistMean %>% 
  mutate(Group = case_when(str_detect(aname, T3 ) ~ group3, TRUE ~ Group)) ##### "True ~ Group" retains previous assigned names that don't match the pattern 
cbindDurDistMean<-cbindDurDistMean %>% 
  mutate(Group = case_when(str_detect(aname, T4 ) ~ group4, TRUE ~ Group))
cbindDurDistMean<-cbindDurDistMean %>% 
  mutate(Group = case_when(str_detect(aname, T5 ) ~ group5, TRUE ~ Group))
cbindDurDistMean<-cbindDurDistMean %>% 
  mutate(Group = case_when(str_detect(aname, T6 ) ~ group6, TRUE ~ Group))




## 4. Export clean data

outputPath<- (paste0 (InputFolder,"/", OutputFolder, "/"))
write.csv(cbindDurDistMean, paste(outputPath, StrainDrug, Replicate, "DurDist", ".csv", sep = "" ))


################################################### Step 4: Integrate replicates, plot, analyse ############################################

## 1. Set up your directory and working data
directory2 <- "G:/Shared drives/Project/Embryo_Swimming_Data/AB_Zn" #This line should take you straightaway to the folder with the filtered files produced in the previous step
setwd(directory2)

data <- list.files(path = directory2 , pattern = "*.csv")

dataDurDist <- readr::read_csv(data, id = "file_name")
dataDurDist <- do.call(cbind, dataDurDist) #do call applies the function to the entire list (not each individual element like lapply)
#rbind takes row names. cbind takes column names instead



dataDurDist <- as.data.frame(dataDurDist)
names(dataDurDist)[1]<- "FileID" # renames the new column



##### Organize the order of your factors per variable #############

dataDurDist$Group <- factor(dataDurDist$Group,levels=c(group1, group2, group3, group4, group5, group6))

##### This part takes a fraction of each on/off phase 
dataDurDist <- dataDurDist %>%
  filter(!(bin %in% c(15:294, 311:594, 611:894, 911:1194, 1211:1494, 1511:1800)))


dataDurDist$alldist<- as.numeric(dataDurDist$alldist)
dataDurDist$bin<- as.numeric(dataDurDist$bin)

##### adjust this to your experiment setup (may have more/less on-off switches)
dataDurDist <- dataDurDist %>%
  mutate(bin_group = case_when(
    bin >= 1 & bin <= 15 ~ "ON 1",
    bin >= 295 & bin <= 310 ~ "OFF 1",
    bin >= 595 & bin <= 610 ~ "ON 2",
    bin >= 895 & bin <= 910 ~ "OFF 2",
    bin >= 1195 & bin <= 1210 ~ "ON 3",
    bin >= 1495 & bin <= 1510 ~ "OFF 3",
    TRUE ~ "Other" 
  ))

dataDurDist$bin_group <- factor(dataDurDist$bin_group,      
                            levels = c("ON 1", "OFF 1", "ON 2", "OFF 2", "ON 3", "OFF 3"))




vline_data <- data.frame(
  bin_group = c( "ON 1", "OFF 1", "ON 2", "OFF 2", "ON 3", "OFF 3"),
  xintercept = c(1, 300, 600, 900, 1200, 1500) 
) ###Change this based on your experiment's timings - Each vline is the time 
  ###at which the light switches in the machine

vline_data$bin_group <- factor(vline_data$bin_group,  
                                levels = c("ON 1", "OFF 1", "ON 2", "OFF 2", "ON 3", "OFF 3"))


ggplot(data = dataDurDist, aes(x= bin, y= alldist, fill= Group)) +
  #stat_summary(geom = "bar", fun = mean, position = position_dodge(0.85), size = 1, width = 0.7) + 
  stat_summary(
    fun = mean,
    geom = 'line', 
    aes(group=Group, colour = Group)) +
  #geom_sina()+
  stat_summary(fun.data = "mean_se",
               shape=18, size=0.25, aes(colour = Group))+ 
  guides(alpha = "none", fill = "none")+  
  scale_y_continuous(expand = c(0, 0))+
  labs(x="Time (s)", y="Distance (mm)")+
  facet_wrap(~bin_group, scales = "free_x", ncol = 6) +
  geom_vline(data = vline_data, aes(xintercept = xintercept), colour = "darkblue", linetype = 2,size = 1) +
  theme(aspect.ratio = 5/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        text = element_text(size = 10),
        legend.text=element_text(size=20),
        panel.spacing = unit(0.5, "cm", data = NULL),
        strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"))

ggsave("AB_Zn_distance.png",  width=15, height=12)

controlandZnMn <- subset(
  dataDurDist,
  Group %in% c("Control", "Zn_1", "Zn_2")
)

ggplot(data = controlandZnMn, aes(x= bin, y= alldist, fill= Group)) +
  #stat_summary(geom = "bar", fun = mean, position = position_dodge(0.85), size = 1, width = 0.7) + 
  stat_summary(
    fun = mean,
    geom = 'line',
    # position = position_dodge(width = 0.9), 
    aes(group=Group, colour = Group)) +
  #geom_sina()+
  stat_summary(fun.data = "mean_se",
               shape=18, size=0.25, aes(colour = Group))+ 
  #position = position_dodge(width = 0.9))+
  guides(alpha = "none", fill = "none")+  
  #scale_y_continuous(guide = "prism_minor", minor_breaks = seq(0, 100, 1), expand = c(0, 0), limits = c(0, 100))+
  labs(x="5 seconds Time bins", y="Distance (mm)")+
  facet_wrap(~bin_group, scales = "free_x", ncol = 6) +
  geom_vline(data = vline_data, aes(xintercept = xintercept), colour = "darkblue", linetype = 2,size = 1) +
  theme(aspect.ratio = 5/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        text = element_text(size = 10),
        legend.text=element_text(size=20),
        panel.spacing = unit(0.5, "cm", data = NULL),
        strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"))

ggsave("AB_Zn_subset.png",  width=15, height=12)
