
# Edited: 03/24/2020

#-----------------------------------------------------------------------------------

clean_data = function(data){
  
  # remove duplicated rows
  data = unique(data)
  
  # remove missing values 
  Gurl = data
  data[Gurl==-300] = NA
  data = na.omit(data)
}


#-----------------------------------------------------------------------------------

print_pca_importance = function(data, k){
  
  # principal component analysis
  pr.out = prcomp(data, scale=TRUE)
  
  # print out first k PC's
  summary(pr.out)$importance[,1:k]
}

#-----------------------------------------------------------------------------------

get_pc_space = function(data, k){
  
  # principal component analysis
  pr.out = prcomp(data, scale = TRUE)
  
  # choose the first k PC's to work with
  PCA_Data = pr.out$x[,1:k]
  
  # transform to new space
  data_pca = PCA_Data %>% as.data.frame()
  
  return(data_pca)
}

# prcomp(data, scale = TRUE)$rotation %>%
#   as.data.frame() %>%
#   rownames_to_column("Predictor") %>%
#   write.csv("./result/pca rotation.csv", row.names = FALSE)

#-----------------------------------------------------------------------------------

get_kmeans_results = function(data_pca, k){
  
  # build a model for k_means
  model = kmeans(data_pca, centers = k, nstart = 25)
  
  # generate clustering result
  clusters = model$cluster
  
  # sort cluster labels based on size
  orders = clusters %>% table() %>% sort(decreasing = TRUE) %>% names()
  
  # re-label clusters based on size
  clusters = sapply(clusters, function(x) which(x == orders))
  
  return(clusters)
}

#-----------------------------------------------------------------------------------

get_kmeans_visual_2D = function(clusters){
  
  # plot the clustering results
  data.frame(Compound = compound, Cluster = as.factor(clusters), data_pca[,1:2]) %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Cluster, shape = Cluster), size = 5, alpha = 0.5) +
    geom_encircle(aes(fill = Cluster), s_shape=1, expand=0, alpha=0.2, show.legend=FALSE) +
    geom_text(aes(label = Compound), size = 4, color="black") + 
    scale_color_manual(values=set_color) + 
    scale_fill_manual(values=set_color) +
    scale_shape_manual(values=1:nlevels(as.factor(clusters))) +
    theme_bw()

}

#-----------------------------------------------------------------------------------

get_kmeans_histogram = function(clusters){
  
  data.frame(Cluster = as.factor(clusters), GroupCat = group_cat_text) %>% 
    table() %>% 
    melt() %>%
    ggplot(aes(x = GroupCat, y = value, fill = as.factor(Cluster))) +
    geom_histogram(stat="identity", position="dodge", color="black", size = 0.1) +
    scale_fill_manual(values = set_color) +
    labs(y="Count", fill="Kmean cluster") +
    theme_bw()
}

#-----------------------------------------------------------------------------------

get_kmeans_visual_3D = function(clusters){
  
  # plot the clustering results
  data.frame(Compound = compound, Cluster = as.factor(clusters), data_pca[,1:3]) %>%
    plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, color = ~Cluster, colors = set_color[1:k], text = ~Compound) %>%
    add_markers(marker=list(size = 7, opacity = 0.6)) %>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')),
           title = paste("K-means in a 3D scatter plot with", k, "clusters"))
}

#-----------------------------------------------------------------------------------

save_kmeans_results = function(data, data_pca, k){
  
  sapply(k, function(x) get_kmeans_results(data_pca,x)) %>% 
    as.data.frame() %>% 
    set_colnames(paste("kmeans_", k, sep="")) %>% 
    cbind(Compound=compound, GroupCat=group_cat, data) %>% 
    select(Compound, GroupCat, everything()) %>% 
    set_rownames(NULL) %>% 
    write.csv("../result/kmeans result.csv", row.names = FALSE)
}

#-----------------------------------------------------------------------------------

get_mclust_model = function(data_pca){
  
  # choose a model based on bic
  BIC = mclustBIC(data_pca)
  
  # build a model 
  model = Mclust(data_pca, x = BIC) 
  
  return(model)
}

#-----------------------------------------------------------------------------------

get_mclust_visual_2D = function(model){

  # generate mclust results
  clusters = model$classification
  
  # sort cluster labels based on size
  orders = clusters %>% table() %>% sort(decreasing = TRUE) %>% names()
  
  # re-label clusters based on size
  clusters = sapply(clusters, function(x) which(x == orders))

  # data.frame(Compound = compound, Cluster = as.factor(clusters), data_pca[,1:2]) %>%
  #   ggplot(aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) +
  #   geom_point(alpha = 0.8) + 
  #   # geom_encircle(aes(fill = Cluster), s_shape=1, expand=0, alpha=0.2, show.legend=FALSE) +
  #   scale_color_manual(values=set_color) + 
  #   scale_fill_manual(values=set_color) +
  #   scale_shape_manual(values=1:nlevels(as.factor(clusters))) +
  #   theme_bw()
  # 
  # data.frame(Compound = compound, Cluster = as.factor(clusters), data_pca[,c(1,3)]) %>%
  #   ggplot(aes(x = PC1, y = PC3, color = Cluster, shape = Cluster)) +
  #   geom_point(alpha = 0.8) + 
  #   geom_encircle(aes(fill = Cluster), s_shape=1, expand=0, alpha=0.2, show.legend=FALSE) +
  #   scale_color_manual(values=set_color) + 
  #   scale_fill_manual(values=set_color) +
  #   scale_shape_manual(values=1:nlevels(as.factor(clusters))) +
  #   theme_bw()
  # 
  # data.frame(Compound = compound, Cluster = as.factor(clusters), data_pca[,c(2,3)]) %>%
  #   ggplot(aes(x = PC2, y = PC3, color = Cluster, shape = Cluster)) +
  #   geom_point(alpha = 0.8) + 
  #   geom_encircle(aes(fill = Cluster), s_shape=1, expand=0, alpha=0.2, show.legend=FALSE) +
  #   scale_color_manual(values=set_color) + 
  #   scale_fill_manual(values=set_color) +
  #   scale_shape_manual(values=1:nlevels(as.factor(clusters))) +
  #   theme_bw()
  
  # plot
  plot(model, what="classification", dimens = c(1,2), col=set_color[1:9])
  plot(model, what="classification", dimens = c(2,3), col=set_color[1:9])
  plot(model, what="classification", dimens = c(1,3), col=set_color[1:9]) 
}

#-----------------------------------------------------------------------------------

get_mclust_histogram = function(model_results){
  # generate mclust results
  clusters = model$classification
  
  # sort cluster labels based on size
  orders = clusters %>% table() %>% sort(decreasing = TRUE) %>% names()
  
  # re-label clusters based on size
  clusters = sapply(clusters, function(x) which(x == orders))
  
  # histogram 
  data.frame(Cluster = as.factor(clusters), GroupCat = group_cat_text) %>% 
    table() %>% 
    melt() %>%
    ggplot(aes(x = GroupCat, y = value, fill = as.factor(Cluster))) +
    geom_histogram(stat="identity", position="dodge", color="black", size = 0.3) +
    scale_fill_manual(values = set_color) +
    labs(y="Count", fill="Mclust cluster") +
    theme_bw()
}

#-----------------------------------------------------------------------------------

get_mclust_visual_3D = function(model){
  
  # generate mclust results
  clusters = model$classification
  
  # sort cluster labels based on size
  orders = clusters %>% table() %>% sort(decreasing = TRUE) %>% names()
  
  # re-label clusters based on size
  clusters = sapply(clusters, function(x) which(x == orders))
  
  # plot
  data.frame(Compound = compound, Cluster = as.factor(clusters), data_pca) %>%
    plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, color = ~Cluster, text = ~Compound) %>%
    add_markers(marker=list(size = 4, opacity = 0.6)) %>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')),
           title = "Mclust in a 3d scatter plot with 9 components")
}
  
#-----------------------------------------------------------------------------------

save_mclust_results = function(data, data_pca, compound, group_cat){
  
  # choose a model based on bic
  BIC = mclustBIC(data_pca)
  
  # build a model 
  model = Mclust(data_pca, x = BIC)
  
  # generate mclust results
  clusters = model$classification
  
  # sort cluster labels based on size
  orders = clusters %>% table() %>% sort(decreasing = TRUE) %>% names()
  
  # re-label clusters based on size
  clusters = sapply(clusters, function(x) which(x == orders))
  
  # save results 
  data.frame(Compound=compound, GroupCat=group_cat, Mclust=as.factor(clusters), data) %>% 
    set_rownames(NULL) %>%
    # ggplot(aes(x = GroupCat, fill = Mclust)) +
    # geom_histogram(stat="count", position="dodge2") +
    # theme_classic()
    write.csv("../result/mclust result.csv", row.names = FALSE)
}

#-----------------------------------------------------------------------------------

check_missing = function(group_cat){
  
  # detect any missing
  if(sum(group_cat == -300) > 1) print("There is at least one missing GroupCat")
}

#-----------------------------------------------------------------------------------

dummy_group_cat = function(group_cat, k){
  
  # convert categories to dummy columns
  group_cat_dummy = data.frame(GroupCat = group_cat) %>% dummy_cols(select_columns = "GroupCat") %>% select(-GroupCat)
  
  # sort column names of chem group based on the size
  descending_order = group_cat_dummy %>% colSums() %>% sort(decreasing = TRUE) %>% names()
  
  # truncate for first k chem group
  group_cat_dummy %>% select(descending_order[1:k])
}

#-----------------------------------------------------------------------------------

get_nn_model = function(group_cat_dummy, data_pca){
  
  # pack a data frame 
  group_data_pca = data.frame(group_cat_dummy, data_pca)
    
  # build a model for first k largest chem groups
  formula_response = paste(names(group_cat_dummy), collapse = " + ")
  formula_predictor = paste(names(data_pca), collapse = " + ")
  set_formula = paste(formula_response, "~", formula_predictor) %>% as.formula()
  
  model = neuralnet(set_formula, data=group_data_pca, hidden=c(20,10,10,5), 
                    linear.output=FALSE, stepmax=1e6, err.fct="ce")
  
  return(model)
}
  
#-----------------------------------------------------------------------------------

get_nn_cv_rate = function(group_cat_dummy, data_pca){

  # split the data
  n = nrow(data_pca)
  folds = n
  folding = split(sample(1:n), 1:folds)
  
  # create a space to save classification rate
  rate = c()
  
  # start cv loop
  for(i in 1:folds){
    # split training and test
    test = folding[[i]]
    
    # train a model of neural nets using train data
    model = get_nn_model(group_cat_dummy = group_cat_dummy[-test,], data_pca = data_pca[-test,])
    
    # predict on test data
    group_cat_pred = 
      compute(model, data_pca[test,])$net.result %>% 
      apply(1, function(x) (x == max(x))+0) %>% 
      t() %>% 
      as.data.frame() %>% 
      set_colnames(names(group_cat_dummy)) %>% 
      set_rownames(NULL)
    
    # calculate and save classification rate
    rate[i] = mean(rowSums(abs(group_cat_pred - group_cat_dummy[test,])) == 0)
  }
  
  return(rate)
}

#-----------------------------------------------------------------------------------

save_nn_results = function(model, compound, group_cat, data, dummy_order){
  
  model$net.result[[1]] %>% 
    apply(1, function(x) which.max(x)) %>% 
    data.frame(nnGroupCat = gsub("GroupCat_","",dummy_order[.])) %>% 
    set_rownames(NULL) %>% 
    cbind(Compound=compound, GroupCat=group_cat, data) %>% 
    write.csv("../result/nn result.csv", row.names = FALSE)
}

#-----------------------------------------------------------------------------------

repeat_B = function(k){
  B = 200
  dummy = c()
  for(i in 1:B){
    group_cat_dummy = dummy_group_cat(group_cat, k = k)
    dummy[i] = get_nn_cv_rate(group_cat_dummy, data_pca) %>% mean()
  }
  return(dummy)
}

#-----------------------------------------------------------------------------------

collapse_data <- function(data){
  data %>% 
    mutate(GroupCat = fct_other(factor(GroupCat), keep = c(3,5,6), other_level = 'Others'),
           GroupCat = factor(GroupCat, 
                             levels = c("3","5","6","Others"),
                             labels = c("Cubic","Tilted","Hexagonal","Others")))
}
  
#-----------------------------------------------------------------------------------

get_coef <- function(cv_model, tuning_parameter){
  coef_packed = coef(cv_model, s = tuning_parameter)
  
  do.call(cbind, coef_packed) %>% 
    as.matrix() %>% as.data.frame() %>% 
    set_colnames(names(coef_packed)) %>% 
    rownames_to_column("feature") %>% 
    arrange(feature) %>% 
    as.data.frame()
}

#-----------------------------------------------------------------------------------

plot_coef <- function(coef_table){
  coef_table %>% 
    reshape2::melt(id.vars = c("feature")) %>%
    ggplot(aes(x = feature, y = value, color = variable)) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, size = 5, alpha = 0.3, color = "grey50") +
    scale_color_nejm() +
    facet_wrap(variable~., nrow=1) +
    labs(x = "", y = "Coefficient", color = "Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    coord_flip()
}

#-----------------------------------------------------------------------------------

plot_coef_heatmap <- function(coef_table){
  coef_table %>% 
    reshape2::melt(id.vars = "feature") %>% 
    mutate(text = paste0("predictor: ", feature, "\n", 
                         "GroupCat: ", variable, "\n", 
                         "Value: ", round(value, 4), "\n")) %>% 
    ggplot(aes(x = variable, y = feature, fill = value, text = text)) +
    geom_tile() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                         high="red", space ="Lab")
}

#-----------------------------------------------------------------------------------

prediction_table <- function(alpha, lambda){
  t = lapply(folds, function(id) {
    X_test = X[id,]; X_train = X[-id,]
    Y_test = Y[id]; Y_train = Y[-id]
    model = glmnet(x = X_train, y = Y_train, alpha = alpha, family = "multinomial")
    Y_pred = predict(model, newx = X_test, type = "class", s = lambda)
    table(factor(Y_pred, levels = levels(factor(data$GroupCat))), 
          factor(Y_test, levels = levels(factor(data$GroupCat))))
  }) 
  r = lapply(t, function(x) sum(diag(x))/sum(x)) %>% unlist()
  t = Reduce("+",t) %>% as.matrix()
  return(list(r=r, t=t))
}

#-----------------------------------------------------------------------------------

print_accurate_tb <- function(r){
  t(r) %>% data.frame() %>% cbind(Mean = mean(r)) %>% kable(escape = F, booktabs = T) %>% kable_styling()
}

#-----------------------------------------------------------------------------------

highlight_tb_count <- function(m){
  total = colSums(m)
  diag(m) = cell_spec(diag(m), 
                      background = ifelse(diag(m)>=0, "red","white"),
                      color = ifelse(diag(m)>=0, "white", "black"),
                      bold = ifelse(diag(m)>=0, T, F))
  rbind(m, Total = total) %>% 
    kable(escape = F, booktabs = T) %>%
    kable_styling()
}

#-----------------------------------------------------------------------------------

highlight_tb_percent <- function(m){
  m2 = sweep(m,2,colSums(m),`/`) %>% round(2)
  diag(m2) = cell_spec(diag(m2), 
                       background = ifelse(diag(m2)>=0, "red","white"),
                       color = ifelse(diag(m2)>=0, "white", "black"),
                       bold = ifelse(diag(m2)>=0, T, F))
  rbind(m2,
        Total = rep("100%", ncol(m))) %>% 
    kable(escape = F, booktabs = T) %>%
    kable_styling()
}

#-----------------------------------------------------------------------------------

remove_identical_cal <- function(data){
  cols_to_remove <- apply(data, 2, function(x) length(unique(x)) == 1)
  data[,!cols_to_remove]
  # data[, !str_detect(names(data), c("X"))]
}















