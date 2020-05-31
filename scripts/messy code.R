library(tidyverse)
library(magrittr)
library(corrplot)

data <- 
  read.csv("data/filledDatabaseNUMONLY_051820.csv") %>% 
  clean_data() %>% 
  filter(GroupCat %in% c(2,3,4,5,6))

#### ----------------   AREN vs Pauling EN   ---------------- ####

electron <- data[, names(data) %>% str_detect("AREN|PaulingEN")]
electron[str_ends(names(electron), "A")] %>% cor() 
electron[str_ends(names(electron), "Aprime")] %>% cor()
electron[str_ends(names(electron), "B")] %>% cor()
electron[str_ends(names(electron), "Bprime")] %>% cor()


#### ----------------------   Radius   ---------------------- ####

radius <- data[, names(data) %>% 
                 str_detect(paste(
                   c("CrystalRadius",
                     "IonicRadii",
                     "ZungerRadii",
                     "r_XII",
                     "r_VI",
                     "r_VII"),
                   collapse = "|"))]

radius[str_ends(names(radius), "A")] %>% cor() %>% corrplot.mixed()
radius[str_ends(names(radius), "Aprime")] %>% cor() %>% corrplot.mixed()
radius[str_ends(names(radius), "B")] %>% cor() %>% corrplot.mixed()
radius[str_ends(names(radius), "Bprime")] %>% cor() %>% corrplot.mixed()



