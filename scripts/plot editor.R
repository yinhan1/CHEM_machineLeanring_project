library(tidyverse)
library(magrittr)
library(fastDummies)
library(ggalt)
library(reshape2)
library(RColorBrewer)
library(ggsci)

library(cluster)
library(factoextra)
library(mclust)
library(neuralnet)
library(dendextend) 

data = read.csv("./data/filledDatabase10042019NUMonlyREDUCED.csv")
source("./R script/functions.R")
data = clean_data(data)

compound = data$Compound
group_cat = data$GroupCat
group_cat_text = paste("Grp", group_cat)
data = select(data, -c("Compound","GroupCat"))

data_biplot = data
data_pca = get_pc_space(data, k = 17) %>% scale()



####----------   compounds in PC space   ----------#### 

set_color = c("#0071C3","#DE501A","#EEB020","#7E2E8E","#79AC2C","#A51331","#4DBDF7") %>% 
  rep(10)

data.frame(Compound = compound, GroupCat = group_cat_text, data_pca) %>%
  ggplot(aes(x=PC1, y=PC2, color = GroupCat)) +
  # geom_point(aes(color = GroupCat), size = 3, alpha = 0.4) +
  geom_text(aes(label=Compound, color=GroupCat), size = 4) +
  scale_color_manual(values=set_color) +
  scale_fill_manual(values=set_color) +
  scale_shape_manual(values=1:11) +
  theme_minimal() +
  labs(color="") +
  ylim(-3,4.5) +
  xlim(-2.5,2.5)



####-----------------   Biplot   ----------------#### 

rownames(data_biplot) = make.names(compound, unique=TRUE)
fit <- princomp(data_biplot, cor=TRUE)
fviz_pca_biplot(fit, aesx = c(1,2),
                # individual
                label = "var", labelsize = 4,
                geom = c("point","text"), 
                fill.ind = group_cat_text, alpha.ind = 0.7,
                pointsize = 1, 
                pointshape = 21, 
                palette = set_color[1:11],
                # variable
                col.var = "contrib", 
                gradient.cols = c("#00AFBB","orange1","#FC4E07"), 
                repel=TRUE) +
  labs(title="",
       fill = "", 
       color = "Contribution",
       x = "PC1",
       y = "PC2") 



####------------   Cubic vs Tilted   -----------#### 

set_color2 = c("#A51331","#0071C3")
data.frame(Compound = compound, GroupCat = group_cat_text, data_pca) %>%
  filter(GroupCat %in% c("Grp 3","Grp 5")) %>%
  ggplot(aes(x=PC1, y=PC2, color = GroupCat)) +
  # geom_point(aes(color = GroupCat), size = 4, alpha = 0.4) +
  geom_text(aes(label=Compound, color=GroupCat), size = 4, alpha = 0.7) +
  scale_color_manual(values=set_color2) +
  scale_fill_manual(values=set_color2) +
  scale_shape_manual(values=1:11) +
  theme_minimal() +
  labs(color="") +
  ylim(-3,4.5) +
  xlim(-2.5,2.5)
