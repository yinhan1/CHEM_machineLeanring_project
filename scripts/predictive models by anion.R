

library(tidyverse)
library(magrittr)
library(ggsci)
library(kableExtra)

library(plotly)

source("scripts/functions.R")

data <- 
  read.csv("data/filledDatabaseNUMONLY_051820.csv") %>% 
  clean_data() %>% 
  filter(GroupCat %in% c(2,3:6))

table(data$X) %>% sort()

#### =====================  Anion O  ===================== #### 

subset <- data %>% filter(X == "O")
X <- subset[,-c(1:4)]
table(subset$GroupCat) %>% sort()

#### -------------   step 0: EDA   ------------- #### 

get_bars <- function(index){
  data.frame(
    group = subset$GroupCat,
    feature = X[,index]
  ) %>% 
    ggplot(aes(x = feature, fill = as.factor(group))) +
    geom_histogram(stat = "bin", bins = 30) +
    scale_fill_jco() +
    theme_minimal() + 
    labs(x = names(X)[index], 
         fill = 'Group')
}

get_bars(1)
get_bars(2)
get_bars(3)
get_bars(4)
get_bars(5) # 4 low, 2 high
get_bars(6)
get_bars(7) # 4 low, 2 high
get_bars(8) # 4 low, 2 high
get_bars(9) # 4 high, 2 low
get_bars(10) # 4 high, 2 low

get_bars(11)
get_bars(12) # 4 high, 2 low
get_bars(13) # 4 high, 2 low ## 
get_bars(14) # 4 high, 2 low
get_bars(15)
get_bars(16)
get_bars(17)  # 4 high, 2 low ## 
get_bars(18)  # 4 high, 2 low ## 
get_bars(19)
get_bars(20)  # 4 high, 2 low 

get_bars(21)  # 4 high, 2 low 
get_bars(22)
get_bars(23)
get_bars(24)
get_bars(25)
get_bars(26)
get_bars(27)
get_bars(28)
get_bars(29)
get_bars(30)  # 4 high, 2 low ## 

get_bars(31)
get_bars(32) ##### 
get_bars(33) #####
get_bars(34) 
get_bars(35) 
get_bars(36)
get_bars(37)
get_bars(38) #
get_bars(39)
get_bars(40)

get_bars(41)
get_bars(42) # 
get_bars(43) 
get_bars(44)
get_bars(45)
get_bars(46)
get_bars(47)
get_bars(48)
get_bars(49)
get_bars(50)

get_bars(51)
get_bars(52)
get_bars(53)
get_bars(54)
get_bars(55)
get_bars(56)
get_bars(57)
get_bars(58)
get_bars(59)
get_bars(60)

get_bars(61)
get_bars(62)
get_bars(63)
get_bars(64) #### 
get_bars(65)


get_bars2 <- function(index){
  data.frame(
    group = subset$GroupCat,
    feature = X[,index]
  ) %>% 
    ggplot(aes(x = feature, fill = as.factor(group))) +
    geom_histogram(stat = "bin", bins = 70) +
    scale_fill_jco() +
    theme_minimal() + 
    labs(x = names(X)[index], 
         fill = 'Group')
}

get_bars(17) + facet_grid(group~.)
get_bars2(32) + facet_grid(group~.)
get_bars2(64) + facet_grid(group~.)

data.frame(
  compound = subset$Compound,
  group = as.factor(subset$GroupCat),
  x17 = X[,17],
  x32 = X[,32],
  x64 = X[,64]
) %>% 
  ggplot(aes(x = x17, y = x32, color = as.factor(group))) +
  geom_point(size = 2.5) +
  scale_color_jco() +
  theme_minimal() +
  labs(x = names(X)[17], y = names(X)[32], color = "Group")
  

#### -------------   step 1: PCA   ------------- #### 

pca <- prcomp(X, scale = TRUE)
X_PC <- pca$x[,1:17]

data.frame(Compound = subset$Compound, 
           Cluster = as.character(subset$GroupCat), 
           X_PC) %>%
  plot_ly(x = ~PC1, 
          y = ~PC3, 
          z = ~PC5, 
          color = ~Cluster, 
          colors = "Dark2", 
          text = ~Compound) %>%
  add_markers(marker=list(size = 7, opacity = 0.6)) %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))







