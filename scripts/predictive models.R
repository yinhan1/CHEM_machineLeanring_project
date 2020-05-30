library(tidyverse)
library(magrittr)
library(ggsci)
library(kableExtra)
library(plotly)

library(glmnet)
library(caret)

library(gbm)

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
         GroupCat = factor(GroupCat, 
                           levels = c("Cubic",
                                      "Tilted",
                                      "Hexagonal",
                                      "LiNb03",
                                      "NCOT")))
table(data$X) %>% sort(decreasing = TRUE)



#### ================================  Anion O  ================================ #### 


subset <- data %>% filter(X == "O")
table(subset$GroupCat) %>% sort(decreasing = TRUE)

X <- subset[,-c(1:4)] %>% 
  remove_identical_cal() %>% 
  as.matrix()
Y <- subset$GroupCat %>% 
  droplevels() %>% 
  as.matrix()


#### ---------------------------   section 0: PCA   --------------------------- #### 


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
PC[,2:4] <- PC[,2:4]*30

PC_initial <- PC %>% mutate(tag = "start")
PC_initial[,2:4] = 0

bind_rows(PC, PC_initial) %>% 
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
  add_trace(
    data = PC_point,
    x = ~PC1, 
    y = ~PC2, 
    z = ~PC3, 
    color = ~Cluster, 
    colors = "Dark2", 
    text = ~Compound,
    type = 'scatter3d', mode = 'markers',
    opacity = 0.9)



#### --------------------   section 1: multinomial reg   -------------------- #### 

subset2 <- data %>% 
  filter(X == "O") %>% 
  filter(GroupCat != "NCOT") %>% 
  droplevels()

X <- subset2[,-c(1:4)] %>% remove_identical_cal() %>% as.matrix()
Y <- subset2$GroupCat %>% droplevels() %>% as.matrix()

set.seed(2020)
folds <- createFolds(1:nrow(X), k = 5, list = TRUE, returnTrain = FALSE)

#### ridge

ridge_cv = cv.glmnet(x = X, y = Y, alpha = 0, 
                     nfolds = 5, 
                     type.measure = "deviance", 
                     family = "multinomial")

tb_ridge = prediction_table(alpha = 0, lambda = ridge_cv$lambda.min) 
tb_ridge$r %>% print_accurate_tb()
tb_ridge$t[,-5] %>% highlight_tb_count()
tb_ridge$t[,-5] %>% highlight_tb_percent()

#### lasso 

lasso_cv = cv.glmnet(x = X, y = Y, alpha = 1, 
                     nfolds = 5, 
                     type.measure = "deviance", 
                     family = "multinomial")

tb_lasso = prediction_table(alpha = 1, lambda = lasso_cv$lambda.min) 
tb_lasso$r %>% print_accurate_tb()
tb_lasso$t[,-5] %>% highlight_tb_count()
tb_lasso$t[,-5] %>% highlight_tb_percent() 

#### elastic net  

elastic_cv <-
  train(GroupCat ~., data = data.frame(X, GroupCat = Y), 
        method = "glmnet",
        trControl = trainControl("cv", number = 5), 
        tuneLength = 10)

tb_elastic = prediction_table(alpha = elastic_cv$bestTune[[1]], 
                              lambda = elastic_cv$bestTune[[2]])
tb_elastic$r %>% print_accurate_tb()
tb_elastic$t[,-5] %>% highlight_tb_count()
tb_elastic$t[,-5] %>% highlight_tb_percent()


#### bootstrap coefficients

# B <- 2000
# dummy_coef <- data.frame(feature = NULL, group_cat = NULL, value = NULL)
# for (i in 1:B){
#   rows_to_take <- get_samples(Y)
#   ridge_cv_B <- cv.glmnet(x = X[rows_to_take,],
#                          y = Y[rows_to_take],
#                          alpha = 0,
#                          nfolds = 5,
#                          type.measure = "deviance",
#                          family = "multinomial")
#  
#   dummy_coef <- rbind(dummy_coef,
#                       ridge_cv_B %>% 
#                         get_coef(tuning_parameter = ridge_cv_B$lambda.min) %>% 
#                         reshape2::melt(id.vars = "feature", variable.name = "group_cat"))
# }

#### plot CI 

# full_join(
#   dummy_coef %>% 
#     group_by(feature, group_cat) %>% 
#     summarise(lower = quantile(value, probs = 0.025),
#               upper = quantile(value, probs = 0.975)),
#   ridge_cv %>% 
#     get_coef(tuning_parameter = ridge_cv$lambda.min) %>% 
#     reshape2::melt(id.vars = "feature", value.name = "est", variable.name = "group_cat")
#   ) %>% 
#   write.csv(file = "data/CI_anionO.csv", row.names = FALSE)
  
read.csv("data/CI_anionO_ridge.csv") %>% 
  filter(!(lower<0 & upper>0)) %>% 
  filter(feature != "(Intercept)") %>% 
  droplevels() %>% 
  ggplot(aes(x = feature, y = est, color = group_cat)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.5, width = 0.5) +
  geom_hline(yintercept = 0, size = 0.2, color = "grey30", linetype = "dashed") +
  scale_color_nejm() +
  facet_wrap(group_cat~., nrow=1) +
  labs(x = "", y = "Coefficient", color = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()
  


#### -------------------------   section 2: GBM   ------------------------- #### 

gbm_cv <- gbm(GroupCat~., data = subset2[,-c(1:3)], 
              shrinkage = 0.01, distribution = "multinomial", 
              cv.folds = 5, n.trees = 3000, verbose = F)
best.iter = gbm.perf(gbm_cv, method="cv")

summary(gbm_cv) %>% View()

fitControl = trainControl(method = "cv", number = 5, returnResamp = "all")
model2 = train(GroupCat~., data = subset2[,-c(1:3)], 
               method = "gbm",
               distribution = "multinomial", 
               trControl = fitControl, verbose = F, 
               tuneGrid = data.frame(.n.trees = best.iter, 
                                     .shrinkage = 0.01, 
                                     .interaction.depth = 1, 
                                     .n.minobsinnode = 1))
tb = confusionMatrix(model2)$table %>% as.matrix()
tb_sum = colSums(tb)  
tb / tb_sum



#### -------------   section 3: what we predict it wrong   ------------- #### 

Y_pred = Y
for(i in 1:length(folds)){
  id = folds[[i]]
  X_test = X[id,]; X_train = X[-id,]
  Y_test = Y[id]; Y_train = Y[-id]
  model = gbm(GroupCat~., data = data.frame(X_train, GroupCat = Y_train),
              shrinkage = 0.01, distribution = "multinomial",
              n.trees = 3000, verbose = F, train.fraction = 1)
  Y_pred[id] = predict(model, n.trees = 3000, newdata = as.data.frame(X_test), type = "response") %>% 
    apply(., 1, which.max)
}

Y_pred = recode(Y_pred, 
                "1" = "Cubic",
                "2" = "Hexagonal",
                "3" = "LiNb03",
                "4" = "Tilted")


top3 <- c("ToleranceBVP", "IonizationPotentialofA", "CrystalRadiusofA")
df_plot <- 
  data.frame(
    Compound = subset$Compound, 
    Cluster = as.character(subset$GroupCat), 
    X[,top3],
    tag = ifelse(Y_pred == Y, 'correct', 'wrong')
  ) 

plot_ly() %>% 
  add_trace(
    data = df_plot %>% filter(tag == 'correct'),
    x = ~ToleranceBVP, 
    y = ~IonizationPotentialofA, 
    z = ~CrystalRadiusofA, 
    color = ~Cluster,
    text = ~Compound,
    type = 'scatter3d', mode = 'markers',
    opacity = 0.8
  ) %>% 
  add_trace(
    data = df_plot %>% filter(tag == 'wrong'),
    x = ~ToleranceBVP, 
    y = ~IonizationPotentialofA, 
    z = ~CrystalRadiusofA, 
    color = ~Cluster,
    text = ~Compound,
    type = 'scatter3d', mode = 'markers',
    marker = list(line = list(color = "red", width = 2, opacity = 0.5)),
    opacity = 0.8
  )


#### ================================  Anion F  ================================ #### 

subset <- data %>% filter(X == "F")
table(subset$GroupCat) %>% sort(decreasing = TRUE)
X <- subset[,-c(1:4)] %>% remove_identical_cal() %>% as.matrix()
Y <- subset$GroupCat %>% droplevels() %>% as.matrix()

#### -----------------------------   step 0: PCA   ----------------------------- #### 

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
PC[,2:4] <- PC[,2:4]*30
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
subset2 <- data %>% 
  filter(X == "F") %>% 
  filter(!(GroupCat %in% c("LiNb03","NCOT"))) %>% 
  droplevels()
X <- subset2[,-c(1:4)] %>% remove_identical_cal() %>% as.matrix()
Y <- subset2$GroupCat %>% droplevels() %>% as.matrix()

set.seed(2020)
folds <- caret::createFolds(1:nrow(X), k = 5, list = TRUE, returnTrain = FALSE)

ridge_cv = cv.glmnet(x = X, y = Y, alpha = 0, nfolds = 5, 
                     type.measure = "deviance", family = "multinomial")
tb_ridge = prediction_table(alpha = 0, lambda = ridge_cv$lambda.min) 
tb_ridge$r %>% print_accurate_tb()
tb_ridge$t[,-(4:5)] %>% highlight_tb_count()
tb_ridge$t[,-(4:5)] %>% highlight_tb_percent()


#### -------------   step 2: lasso   ------------- #### 
lasso_cv = cv.glmnet(x = X, y = Y, alpha = 1, nfolds = 5, 
                     type.measure = "deviance", family = "multinomial")
tb_lasso = prediction_table(alpha = 1, lambda = lasso_cv$lambda.min) 
tb_lasso$r %>% print_accurate_tb()
tb_lasso$t[,-(4:5)] %>% highlight_tb_count()
tb_lasso$t[,-(4:5)] %>% highlight_tb_percent() 

lasso_cv %>% 
  get_coef(tuning_parameter = lasso_cv$lambda.min) %>% 
  View()

#### -------------   step 3: elastic net   ------------- #### 
elastic_cv <-
  train(GroupCat ~., data = data.frame(X, GroupCat = Y), method = "glmnet",
        trControl = trainControl("cv", number = 5),
        tuneLength = 10)
tb_elastic = prediction_table(alpha = elastic_cv$bestTune[[1]], 
                              lambda = elastic_cv$bestTune[[2]])
tb_elastic$r %>% print_accurate_tb()
tb_elastic$t[,-(4:5)] %>% highlight_tb_count()
tb_elastic$t[,-(4:5)] %>% highlight_tb_percent()

elastic_cv$finalModel %>% 
  get_coef(tuning_parameter = elastic_cv$bestTune$lambda) %>% 
  View()

#### -------------   step 4: GBM   ------------- #### 
gbm_cv <- gbm(GroupCat~., data = subset2[,-c(1:3)], 
              shrinkage = 0.01, distribution = "multinomial", 
              cv.folds = 5, n.trees = 3000, verbose = F)
best.iter = gbm.perf(gbm_cv, method="cv")
summary(gbm_cv) %>% View()

fitControl = trainControl(method = "cv", number = 5, returnResamp = "all")
model2 = train(GroupCat~., data = subset2[,-c(1:3)], method = "gbm",
               distribution = "multinomial", trControl = fitControl, verbose=F, 
               tuneGrid = data.frame(.n.trees = best.iter, 
                                     .shrinkage=0.01, 
                                     .interaction.depth=1, 
                                     .n.minobsinnode=1))
model2
tb = confusionMatrix(model2)$table %>% as.matrix()
tb_sum = colSums(tb)  
tb / tb_sum

#### -------------   step 5: what we predict wrong  ------------- #### 
X <- subset[,-c(1:4)] %>% remove_identical_cal() %>% as.matrix()
Y <- subset$GroupCat %>% droplevels() %>% as.matrix()

top3 <- c("ToleranceBVP","DensityatSpecofA","CrystalRadiusofBprime")
X_top3 <- X[,top3]

data.frame(
  Compound = subset$Compound, 
  Cluster = as.character(subset$GroupCat), 
  X_top3
) %>% 
  plot_ly() %>% 
  add_trace(
    x = ~ToleranceBVP, 
    y = ~DensityatSpecofA, 
    z = ~CrystalRadiusofBprime, 
    color = ~Cluster, 
    colors = "Paired", 
    text = ~Compound,
    type = 'scatter3d', mode = 'markers',
    opacity = 0.8
  )


