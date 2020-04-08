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
library(stringr)

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
##              Importing data 
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Supercooling data. Importing and merging into one file
files <- list.files(path = "data/tissue_test_results", full.names = T, pattern = "output")
coldtol_orig <- files %>% map_dfr(~read.csv(.))

sp_id <- read_csv("data/FreezeTrial.csv") # links sample ID to plate position on supercooling run
sp_id2 <-read_csv("data/FreezeTrial_20200403.csv")#updated file that has tissue mass entered for more samples. Ideally we'd use this updated dataset instead of sp_id becaues it has more complete tissue weights, but sp_id2 doesn't work nice in the code below - something must not merge correctly and I havn't figured out what yet.  
#One issue resolved: For some of the dates in the updated file and 0 was added to the front of the month in the date so it doesn't merge well with other datasets. Fucking excel and it's horrible use of dates!!!
#formating Date so it matches the other datasets
sp_id2$Date_Froze = str_remove(sp_id2$Date_Froze, "^0+")
# but still drops some reps when merged and I'm not sure why

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
r <- sp_id %>% 
  select(-Trial_order, -position, -Date_Weighted, -Initials, -Notes, -old_Tissue) %>% 
  mutate(Plate_Num = as.numeric(Plate_Num)) %>% 
  filter(`Experiment Type` == "Seedling Freeze Trial") %>% 
  group_by(Date_Froze, Trial_num, Species, Species_Rep, Tissue, Tissue_Rep, Dry_Mass_mg) %>% 
  filter(!duplicated(Plate_Num))

# merge video id onto the species information sheet
b <- left_join(r, date, by = c("Date_Froze" = "Trial_date", "Trial_num" = "Trial_num"))

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
#pdf("/Users/laura/Desktop/Writing Projects/tissue cold tol/R/tissue_cold_tol/output/ColdTol_byspp_wletters.pdf", width = 15, height = 7) #This saves the pdf

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
        axis.title.y=element_text(size=13, hjust=0.5),
        panel.border = element_rect(colour = "black", fill = NA, size =1)) +
  annotate("text", x = 2.75, y = -14.5, label = "B") + #ASCSYR
  annotate("text", x = 3, y = -12, label = "AB", size = 3) +  #ASCSYR
  annotate("text", x = 3.25, y = -9, label = "A") + #ASCSYR
  annotate("text", x = 4.75, y = -13, label = "B") + #DESILL
  annotate("text", x = 5, y = -12.2, label = "AB", size = 3) + #DESILL
  annotate("text", x = 5.25, y = -9, label = "A") + #DESILL
  annotate("text", x = 5.75, y = -11.5, label = "AB", size = 3) + #HELHEL
  annotate("text", x = 6, y = -13.7, label = "B") + #HELHEL
  annotate("text", x = 6.25, y = -9.05, label = "A") + #HELHEL
  annotate("text", x = 7.75, y = -15.5, label = "B") + #OXAVIO
  annotate("text", x = 8, y = -10.5, label = "A") + #OXAVIO
  annotate("text", x = 8.25, y = -9.22, label = "A") + #OXAVIO
  annotate("text", x = 8.75, y = -16, label = "B") + #PREALB
  annotate("text", x = 9, y = -11, label = "A") + #PREALB
  annotate("text", x = 9.25, y = -9.22, label = "A") + #PREALB
  annotate("text", x = 9.75, y = -12, label = "B") + #RATPIN
  annotate("text", x = 10, y = -9.3, label = "AB", size = 3) + #RATPIN
  annotate("text", x = 10.25, y = -8.95, label = "A") + #RATPIN
  annotate("text", x = 10.75, y = -16, label = "B") + #SOLGRA
  annotate("text", x = 11, y = -11, label = "A") + #SOLGRA
  annotate("text", x = 11.25, y = -9.7, label = "A") + #SOLGRA
  annotate("text", x = 13.75, y = -18, label = "B") + #SYMPIL
  annotate("text", x = 14, y = -10.4, label = "A") + #SYMPIL
  annotate("text", x = 14.25, y = -8.5, label = "A")  #SYMPIL
dev.off()

## Graphing tissues separately: This helps show how seedling were more cold tolerant but also more variable. Edits: get rid of boxes around tissue types; get rid of legend; Possibly come up with a better way to show the differences in variation. 
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
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()

## raw data graphed by tissue type. not that compelling
ggplot(data = cold_all_short %>% filter(tot_reps>4), 
       aes(x = Species, y = med_supercooling, fill = factor(Tissue_combined))) +
  geom_point(aes(color = Tissue_combined)) +
  ylab("Cold tolerance\n(median supercooling temp "*~degree~"C)") +
  xlab("Species") +
  guides(fill=guide_legend(title="Tissue")) +
  theme_classic() +
  facet_wrap(vars(Tissue_combined)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# basic graph of cold tolerance of different species
ggplot(data = cold_all_short %>% filter(tot_reps>4), 
       aes(x = Species, y = med_supercooling)) +
  geom_point(aes(color = Tissue_combined)) +
  xlab("Species") +
  theme_classic()

# basic graph of cold tolerance across the different videos
ggplot(data = cold_all_short %>% filter(Tissue == "Seedling"), 
       aes(x = Date_Froze, y = med_supercooling, color = Species)) +
  geom_point() +
  xlab("Species") +
  geom_smooth(method = "lm", se = FALSE) +
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

#This code gives me believable differences but I'm not sure if it's the correct way to do it becaues it only looks at Species:Tissue comparisons. Those are the only comparissons I'm interested in anyways, but it feels wierd not to have the other factors in the model. 
posthoc2 <- TukeyHSD(aov(med_supercooling ~
                           Species:Tissue_combined, 
                         data = cold_all_short%>% filter(tot_reps>4)))
print(posthoc2)

# adjusting the data a bit so it's easier to look at the comparisions I want  
m <- do.call(rbind.data.frame, posthoc2) %>% #turns it into a dataframe
  mutate(comp = rownames(.)) %>%
  separate(comp, c("first", "second"), sep = "\\.") %>% # split by period
  separate(second, c("spp1", "tiss1", "spp2", "tiss2"), sep = "([\\:\\-])") %>%
  mutate(same_spp = if_else(spp1 == spp2, 1, 0)) %>% #finding species matches
  filter(same_spp == 1)

#saving posthoc results so I can look at them later
write_csv(m, "/Users/laura/Desktop/Writing Projects/tissue cold tol/R/tissue_cold_tol/output/tissue_spp_posthoc.csv")

## This is my first try at the posthoc tests to see where the differences lie. It doesn't seem to be working here and I'm not sure why. One difference between this and the test above is that this uses the model object aov while the model above was a lm - not sure if that matters too much
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
# Using sp_id2 because it has more complete mass measurments than sp_id

# fOR the first few trials, two pixel points were selected in the software to calculate super cooling. The code below averages the readings from those two points
p <- sp_id2 %>% 
  select(-Trial_order, -position, -Date_Weighted, -Initials, -Notes, -old_Tissue) %>% 
  mutate(Plate_Num = as.numeric(Plate_Num)) %>% 
  filter(`Experiment Type` == "Seedling Freeze Trial") %>% 
  group_by(Date_Froze, Trial_num, Species, Species_Rep, Tissue, Tissue_Rep, Dry_Mass_mg) %>% 
  filter(!duplicated(Plate_Num))

# merge video id onto the species information sheet
j <- left_join(p, date, by = c("Date_Froze" = "Trial_date", "Trial_num" = "Trial_num"))

# Merge species id to supercooling by mergeing by video id and plate position. Also correcting spelling errors and grouping cotyledons with seedlings
coldtol$sample <- as.numeric(coldtol$sample)
cold_all_b <- full_join(coldtol, j, by = c("video" = "File_name", "sample" = "Plate_Num")) %>% 
  mutate(Species = plyr::mapvalues(Species, from = c("PPREALB", "SYMPOBL", "ALLACER"), to = c("PREALB", "SYMOBL", "ALLCER"))) %>% 
  filter(Species != "empty", !is.na(Species)) %>% 
  mutate(Tissue_combined = plyr::mapvalues(Tissue, from = c("Seedling", "Seedling ", "Cotyledon"), to = c("Seedling", "Seedling", "Seedling"))) %>% 
  filter(Tissue_combined != "First_seedling_leaf") %>% 
  mutate(Tissue_combined = factor(Tissue_combined, levels = c("Seedling", "Leaf", "Root")))

size <- cold_all_short %>%
  filter(Tissue == "Seedling") # we don't want cotyledons so using "Tissue"

dev.off() #cleaning R, just in case
pdf("/Users/laura/Desktop/Writing Projects/tissue cold tol/R/tissue_cold_tol/output/seedling_ColdTol_vs_mass.pdf", width = 10, height = 5) #This saves the pdf

## Mass doesn't seem to be tighly correlated with seedling cold tolerance
ggplot(data = cold_all_b %>% filter(Tissue == "Seedling", Species !="ALLCER", Species != "LATVEN", Species != "LIACYL"), 
       aes(x = Dry_Mass_mg, y = med_supercooling, color = Species)) +
  geom_point() +
  xlab("Seedling Dry Mass (mg)") +
  ylab("Seedling cold tolerance (C)") +
  #geom_smooth(method = "auto") +
  xlim(0, 15) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(vars(Species)) +
  theme_classic()  
dev.off()


## date froze... code not really working
ggplot(data = cold_all_b %>% filter(Tissue == "Seedling", Species !="ALLCER", Species != "LATVEN", Species != "LIACYL"), 
       aes(x = factor(Date_Froze), y = med_supercooling, color = Species)) +
  geom_boxplot() +
  xlab("Seedling Dry Mass (mg)") +
  ylab("Seedling cold tolerance (C)") +
  #geom_smooth(method = "auto") +
  xlim(0, 15) +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(vars(Species)) +
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







