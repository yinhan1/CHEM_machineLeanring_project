library(tidyverse)
library(magrittr)
library(ggsci)
library(kableExtra)

library(glmnet)
library(caret)

library(plotly)

source("scripts/functions.R")

data <- 
  read.csv("data/filledDatabaseNUMONLY_051820.csv") %>%
  clean_data() %>% 
  filter(GroupCat %in% c(2,3,4,5,6)) %>% 
  mutate(GroupCat ,
         GroupCat = recode(GroupCat, 
                           "2" = "LiNb03", 
                           "4" = "NCOT",
                           "3" = "Cubic", 
                           "5" = "Tilted",
                           "6" = "Hexagonal"),
         GroupCat = factor(GroupCat, levels = c("Cubic","Tilted","Hexagonal","LiNb03","NCOT")))

table(data$X) %>% sort(decreasing = TRUE)

#### =====================  Anion O  ===================== #### 

subset <- data %>% filter(X == "O")
table(subset$GroupCat) %>% sort(decreasing = TRUE)

subset <- data %>% filter(X == "O") %>% filter(GroupCat != "NCOT") %>% droplevels()

X <- subset[,-c(1:4)] %>% remove_identical_cal() %>% as.matrix()
Y <- subset$GroupCat %>% droplevels() %>% as.matrix()
folds <- caret::createFolds(1:nrow(X), k = 5, list = TRUE, returnTrain = FALSE)

#### -------------   step 0: PCA   ------------- #### 

pca <- prcomp(X, scale = TRUE)
summary(pca) 
X_PC <- pca$x[,1:17] %>% as.matrix()

PC_point <- data.frame(Compound = subset$Compound, 
           Cluster = as.character(subset$GroupCat), 
           X_PC)

PC <- pca$rotation %>%
  as.data.frame() %>% 
  select(PC1,PC2,PC3) %>% 
  rownames_to_column("variable") %>% 
  mutate(tag = "end",
         contribution = PC1^2 + PC2^2 + PC3^2)
PC[,2:4] <- PC[,2:4]*20
PC_initial <- PC %>% mutate(tag = "start")
PC_initial[,2:4] = 0
PC_arrow <- bind_rows(PC, PC_initial)

PC_arrow %>% 
  group_by(variable) %>% 
  plot_ly() %>% 
  add_trace(
    x = ~PC1, 
    y = ~PC2, 
    z = ~PC3,
    color = ~contribution,
    text = ~variable,
    type = 'scatter3d', mode = 'lines', opacity = 1, 
    line = list(width = 6, reverscale = FALSE)
  ) %>% 
  add_markers(
    data = PC_point,
    x = ~PC1, 
    y = ~PC2, 
    z = ~PC3, 
    color = ~Cluster, 
    colors = "Dark2", 
    text = ~Compound,
    opacity = 0.9)


#### -------------   step 1: ridge   ------------- #### 

ridge_cv = cv.glmnet(x = X, y = Y, alpha = 0, nfolds = 5, type.measure = "deviance", family = "multinomial")
tb_ridge = prediction_table(alpha = 0, lambda = ridge_cv$lambda.min) 

tb_ridge$r %>% print_accurate_tb()

tb_ridge$t %>% highlight_tb_count()
tb_ridge$t[,-5] %>% highlight_tb_percent()


#### -------------   step 2: lasso   ------------- #### 
lasso_cv = cv.glmnet(x = X, y = Y, alpha = 1, nfolds = 5, type.measure = "deviance", family = "multinomial")
tb_lasso = prediction_table(alpha = 0, lambda = lasso_cv$lambda.min) 

tb_lasso$r %>% print_accurate_tb()

tb_lasso$t %>% highlight_tb_count()
tb_lasso$t[,-5] %>% highlight_tb_percent()


#### -------------   step 3: elastic net   ------------- #### 
elastic_cv <-
  train(GroupCat ~., data = data.frame(X, GroupCat=Y), method = "glmnet",
        trControl = trainControl("cv", number = 5),
        tuneLength = 10)
tb_elastic = prediction_table(alpha = elastic_cv$bestTune[[1]], lambda = elastic_cv$bestTune[[2]])

tb_elastic$r %>% print_accurate_tb()

tb_elastic$t %>% highlight_tb_count()
tb_elastic$t[,-5] %>% highlight_tb_percent()









