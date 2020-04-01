## Gathering supercooling temperature data from thermal camera output
## March 6, 2020
## LML


## Notes from Jon:
#This is a folder with a separate file fore each video. You'll have to combine them. I've been using this code if you have the output files in the "data" folder of your R project. You need tidyverse installed for this to work.

#files <- list.files(path = "data", full.names = T, pattern = "output")
#coldtol <- files %>% map_dfr(~read.csv(.))

#The "sample" column refers to the well that the sample was in. For the supercooling value, I would use the med_supercooling, as that is the median and should be less sensitive to potentially weird measurements. I would filter out any sample that has a "n.sample" less than 5 or so and be a bit weary of samples that have an sc_var higher than 0.1 or so.

#Jon
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

#Load packges
#install.packages("TDPanalysis")
library(tidyverse)
#library(nlme)
library(lme4)
#library(foreign)
#library(agricolae)
library(TDPanalysis) #date to DOY conversion

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
##              Importing data 
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Supercooling data. Importing and merging into one file
files <- list.files(path = "data/tissue_test_results", full.names = T, pattern = "output")
coldtol_orig <- files %>% map_dfr(~read.csv(.))

sp_id <- read_csv("data/FreezeTrial.csv") # links sample ID to plate position on supercooling run
date <- read_csv("data/Video_tracking.csv") #link date to video name
plant_date <- read_csv("data/PlantingDate_TissueColdTol.csv") # Planting date
emergence <- read_csv("data/SeedlingEmergence_TissueColdTol.csv") # Seedling emergence


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
## Question 1: Does cold tolerance vary with regard to species or tissue type?         
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

# Supercooling data
## NOTE this current version doesn't have all the dry biomass yet. Ben G is currently helping to finish enter that. march 11, 2020 lml
coldtol <- dplyr::filter(coldtol_orig, n.samples > 5) # toss low sample sizes, if n.sample < 5
coldtol <- dplyr::filter(coldtol, sc.var < 0.1) # Toss highly variable measurements, if sc_var > 0.1

# fOR the first few trials, two pixel points were selected in the software to calculate super cooling. The code below averages the readings from those two points
sp_id <- sp_id %>% 
  select(-Trial_order, -position, -Date_Weighted, -Initials, -Notes, -old_Tissue) %>% 
  mutate(Plate_Num = as.numeric(Plate_Num)) %>% 
  filter(`Experiment Type` == "Seedling Freeze Trial") %>% 
  group_by(Date_Froze, Trial_num, Species, Species_Rep, Tissue, Tissue_Rep, Dry_Mass_mg) %>% 
  filter(!duplicated(Plate_Num))

# merge video id onto the species information sheet
b <- left_join(sp_id, date, by = c("Date_Froze" = "Trial_date", "Trial_num" = "Trial_num"))

# Merge species id to supercooling by mergeing by video id and plate position. Also correcting spelling errors and grouping cotyledons with seedlings
coldtol$sample <- as.numeric(coldtol$sample)
cold_all <- full_join(coldtol, b, by = c("video" = "File_name", "sample" = "Plate_Num")) %>% 
  mutate(Species = plyr::mapvalues(Species, from = c("PPREALB", "SYMPOBL", "ALLACER"), to = c("PREALB", "SYMOBL", "ALLCER"))) %>% 
  filter(Species != "empty", !is.na(Species)) %>% 
  mutate(Tissue_combined = plyr::mapvalues(Tissue, from = c("Seedling", "Seedling ", "Cotyledon"), to = c("Seedling", "Seedling", "Seedling"))) %>% 
  filter(Tissue_combined != "First_seedling_leaf") %>% 
  mutate(Tissue_combined = factor(Tissue_combined, levels = c("Seedling", "Leaf", "Root")))

## For each root and leaf, several measurements were taken. This code averages those multiple measurements. I make a separate dataset fore each tissue type and merge them all back together. It seems like the clunky way of doing things, but I think it still works
cold_all_l <- cold_all %>%
  filter(Tissue_combined == "Leaf") %>%
  group_by(Species, Species_Rep, Tissue_combined, 
           video, Date_Froze, Tissue, Dry_Mass_mg) %>%
  summarise(med_supercooling = mean(med_supercooling, na.rm = TRUE)) %>%
  ungroup()
  
cold_all_r <-cold_all %>%
  filter(Tissue_combined == "Root") %>%
  group_by(Species, Species_Rep, Tissue_combined, 
           video, Date_Froze, Tissue, Dry_Mass_mg) %>%
  summarise(med_supercooling = mean(med_supercooling, na.rm = TRUE)) %>%
  ungroup()

cold_all_s <- cold_all %>%
  filter(Tissue_combined == "Seedling") %>%
  select(Species, Species_Rep, Tissue_combined, 
         video, Date_Froze, Tissue, Dry_Mass_mg, med_supercooling)

## Merge tissue types back together 
blah <- dplyr::bind_rows(cold_all_l, cold_all_r)
cold_all_short <- dplyr::bind_rows(blah, cold_all_s)

# summarize supercooling data by species and tissue type just because I was curious and wanted to see the numbers
cold_summ <- cold_all_short %>% 
  group_by(Tissue_combined, Species) %>% 
  summarise(reps = n(),
            ave = mean(med_supercooling, na.rm = TRUE)) 

#making a total rep column to add to full dataset so I can filter by it
reps <- cold_all_short %>% 
  drop_na(med_supercooling) %>%
  group_by(Tissue_combined, Species) %>% 
  summarise(tot_reps = n())

cold_all_short <- left_join(cold_all_short, reps)

## ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# GRAPHING
## ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

## Unhash these next two lines when you want to save a graph. Possibly change the name if you don't want to overwrite the last one
#dev.off() #cleaning R, just in case
#pdf("/Users/laura/Desktop/Writing Projects/tissue cold tol/R/tissue_cold_tol/output/ColdTol_byspp.pdf", width = 10, height = 5) #This saves the pdf
## basic graph of the data. I don't really care for hte box plot format (Changes: fix colors; species codes or full spcies names; shift graph so y legend fits; center justify yaxis so cold tolerance is centered over the middle of the second line of text)
ggplot(data = cold_all_short %>% filter(tot_reps>4), 
       aes(x = Species, y = med_supercooling, fill = factor(Tissue_combined))) +
  geom_boxplot(notch = FALSE, varwidth = FALSE) +
  ylab("Cold tolerance\n(median supercooling temp "*~degree~"C)") +
  xlab("Species") +
  guides(fill=guide_legend(title="Tissue")) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1),
        plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
        axis.title.y=element_text(size=13, hjust=0.5)) 
dev.off()

## Graphing tissues separately
#dev.off() #cleaning R, just in case
#pdf("/Users/laura/Desktop/Writing Projects/tissue cold tol/R/tissue_cold_tol/output/ColdTol_byTissue.pdf", width = 10, height = 5) #This saves the pdf
ggplot(data = cold_all_short %>% filter(tot_reps>4), 
       aes(x = Species, y = med_supercooling, fill = factor(Tissue_combined))) +
  geom_boxplot(notch = FALSE, varwidth = FALSE) +
  ylab("Cold tolerance\n(median supercooling temp "*~degree~"C)") +
  xlab("Species") +
  guides(fill=guide_legend(title="Tissue")) +
  theme_classic() +
  facet_wrap(vars(Tissue_combined)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))
dev.off()

# basic graph of cold tolerance of different species
ggplot(data = cold_all_short %>% filter(tot_reps>4), 
       aes(x = Species, y = med_supercooling)) +
  geom_point(aes(color = Tissue_combined)) +
  xlab("Species") +
  theme_classic()

# basic graph of cold tolerance across the different videos
ggplot(data = cold_all_short, 
       aes(x = video, y = med_supercooling)) +
  geom_point(aes(color = Tissue_combined)) +
  xlab("Species") +
  theme_classic()

## *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Statistical test
## *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Does cold tolerance differ with species, tissue, and the interaction? Video was also included to account for the variation of different runs. Video could (should?) be a random effect but in this model it is fixed
mod_tiss  <- lm (med_supercooling ~
        Species +
          video +
          Tissue_combined +
          Species:Tissue_combined,
        data = cold_all_short %>% filter(tot_reps>4))
        
print(mod_tiss)
summary(mod_tiss)
anova(mod_tiss)

## This is the posthoc tests to see where the differences lie
posthoc <- TukeyHSD(aov(med_supercooling ~
                          video +
                          Species +
                          Tissue_combined +
                          Species:Tissue_combined, 
                        data = cold_all_short%>% filter(tot_reps>4)))
print(posthoc)

#The results above are too much to take in at once. This cleans and organizes the results so only looking at differences between tissues of the same species
## NOTE: If the number or reps is changed then make sure to check the indexing in the code below 
## These results don't make much scence. From what I understand, there are only sig diffs with DESILL leaf-seedling, HELHEL leaf-seedling and DESIL root-leaf, but when looking at the graph it doesn't look right
k <- do.call(rbind.data.frame, posthoc) %>% #turns it into a dataframe
  mutate(comp = rownames(.)) %>%
  slice(-c(1:253)) %>% #tossing all video comparisions
  slice(-c(1:108)) %>% #now taking out all spp comparisions and tissue
  separate(comp, c("first", "second"), sep = "\\.") %>% # split by period
  separate(second, c("spp1", "tiss1", "spp2", "tiss2"), sep = "([\\:\\-])") %>%
  mutate(same_spp = if_else(spp1 == spp2, 1, 0)) %>% #finding species matches
  filter(same_spp == 1)
  

  
  


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Question 2: Is cold tolerance related to seedling emergence time?
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

# Prepping emergence data for merge
emerge <- emergence %>%
  filter(!is.na(emergence_date)) %>% #take out reps that didn't sprout
  separate(emergence_date, c("year", "month", "day"), sep = "\\-") %>% # split by hash
  unite(emerg_date, c("day", "month", "year"), sep = "/") %>% #merged date back together
  mutate(emerge_doy = date.to.DOY(emerg_date, format = "dd/mm/yyyy")) %>% #create doy 
  mutate(Species = plyr::mapvalues(Species, from = c("OXACIO", "SYMONOV"), to = c("OXAVIO", "SYMNOV"))) #correct a few spelling errors

#prepping planting date data to merge
plant <- plant_date %>%
  unite(Rep, c("rep", "run"), sep = "") %>% #creates mergable value
  mutate(spp = toupper(spp)) #turn lowercase codes into all uppercase

# prepping cold tolerance data for merge  
cold_all_seed <- cold_all_short %>% filter(Tissue_combined == "Seedling") 

#Merging data and calculating emergence time (sprouttime). Note, based on the final merge it looks like a decent amount of data was tossed. Might consider loosening the varience constraints in the first clip of the cold tolerance data
timing <- full_join(plant, emerge, by = c("spp" = "Species", "Rep")) %>%
  mutate(sprouttime = emerge_doy - plant_doy) %>% #calculate time to sprout
  full_join(., cold_all_seed, by = c("spp" = "Species", "Rep" = "Species_Rep")) #join in cold
  
# making a summary to graph
timing_sum <- timing %>%
  group_by(spp) %>%
  summarise(cold = mean(med_supercooling, na.rm = TRUE), 
            sprout = mean(sprouttime, na.rm = TRUE))

## A box plot that shows variation along both yand x axes would be good.
ggplot(data = timing_sum %>% filter(spp != "AQUCAN"),
       aes(x = sprout, y = cold))+
  geom_point() +
  ylim(-20, -7) +
  xlim(3, 20) +
  geom_smooth(method = "lm") +
  theme_classic ()

q2_mod <- lm(cold ~ sprout, data = timing_sum %>% filter(spp != "AQUCAN"))       
q2_mod <- lm(cold ~ log(sprout), data = timing_sum)  
print(q2_mod)
summary(q2_mod)

q2_cor <- cor(x = timing_sum$sprout, y = timing_sum$cold)
print(q2_cor)
summary(q2_cor)

# graphing to see how time to emergence relates to cold tolerance of seedligns
# ISSUES: it'g graphing species that don't have data...
ggplot(data = timing, #%>% filter(spp != "AQUCAN"), 
       aes(x = log(sprouttime), y = med_supercooling)) +
  geom_point(aes(color = spp)) +
  #stat_summary(
  #  geom = "point",
  #  fun.y = "mean", fun.x = "mean", #this is wrong. Doesn't do mean for each spp
  #  col = "black", #fill = spp,
  #  size = 3, shape = 24
  #) +
  xlab("log Time to seedling emergence (days)") +
  ylab("Seedling cold tolerance (C)") +
  geom_smooth(method = "lm") +
  theme_classic()  

#*~*~*~*~*~*~*~*~*~*~*~
# Statistical analysis
#~*~*~*~*~*~*~*~*~*~*~

q2_all_mod <-lm(med_supercooling ~ sprouttime + spp + video, data = timing %>% filter(spp != "AQUCAN"))

summary(q2_all_mod)
anova(q2_all_mod)



#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Question 3: Sources of variation
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
#Seedlings appear to have a lot more varience in cold tolerance than leaves and roots. Is this a true pattern or is it related to sample size (a lot more seedligns were measured)





#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Question 4: Does size of seedling relate to cold tolerance?
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
#There's a lot of variation in seedling cold tolerance, so we were tyring to understand some of this variation. Seedling size mgiht be one factor.
# For this we want to use only seedlings, not cotyleadons, since their mass would be different

##****WARNING**** AS OF THIS MORNING, MISSING SOME DRY MASSES FOR SEEDLINGS, ESPECIALLY TRAOHI. LOTS OF MASSES MISSING FOR LEAVES AND ROOTS, TOO
## Update: These data have been entered, need to look on google drive to get the most uptodate version

size <- cold_all_short %>%
  filter(Tissue == "Seedling") # we don't want cotyledons so using "Tissue"
  

ggplot(data = size, #%>% filter(spp != "AQUCAN"), 
       aes(x = Dry_Mass_mg, y = med_supercooling)) +
  geom_point(aes(color = Species)) +
  xlab("Seedling Dry Mass (mg)") +
  ylab("Seedling cold tolerance (C)") +
  #geom_smooth(method = "auto") +
  xlim(0, 15) +
  facet_grid(Species) +
  theme_classic()  






#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Question 5: does age of the seedling matter?
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# calculating seedling age so I can add it to the model

#***INCOMPLETE*** NEED TO ALSO CONNECT SEEDLINGS BASED ON TISSUE REP, BUT THERE ISN'T A SIMILAR COLUMN IN THE EMERGENCE DATA....
r <- cold_all_short %>%
  filter(Tissue_combined == "Seedling") %>% #just getting seedlings
  mutate(froze_doy = date.to.DOY(Date_Froze, format = "mm/dd/yyyy")) %>% #create doy 

ce <- full_join(cold_all_short, emerge, by = c("Species", "Species_Rep" = "Rep"))







