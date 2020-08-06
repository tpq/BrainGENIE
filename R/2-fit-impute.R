# 5. Fit LR-PCA prediction model to GTEx counts, then apply the resulting weights to the PCA scores derived from new samples
fit_lr_weights_in_gtex = function(pca_model  = NULL, tissue = NULL, gene_list = NULL, n_comps = 20){
  if(is.null(tissue)){stop("Please specify a tissue based on nomenclature available in the cross-validation table")}
  if(is.null(pca_model)){stop("Argument required for blood_PCA")}
  if(is.null(gene_list)){stop("Please specify gene IDs to include in the training process.")}
  # TODO: Linear regression: brain gene ~ blood PCA
  Y = data.frame(t(brain_expr)) # select normalized GTEx brain counts
  # filter Y by genes that are well predicted from 5-fold cross-validation
  cv_performance = gene_list
  matched_genes = intersect(colnames(Y), cv_performance)
  if(length(matched_genes) < 1){stop("Unexpected issue! No matching gene IDs between CV performance file and normalized GTEx counts.")}
  common_genes = intersect(matched_genes, pca_model$genes)
  if(length(matched_genes) < 1){stop("Unexpected issue! No genes available to run LR.")}
  message("\rTraining LR models for: ", length(matched_genes), " genes")

  Y = Y[,colnames(Y) %in% matched_genes, ]
  X = data.frame(pca_model$pca$x[,1:n_comps]) # use PCA model derived from GTEx paired blood-brain data
  fit = lm(as.matrix(Y) ~ ., data = X) # fit a LR model per gene
  return(fit) # return model
}

# 6. apply brain transcriptome prediction model to blood PCA
impute_gxp = function(pca_in_new_sample = NULL, trained_model = NULL, scale = TRUE){

  if(is.null(trained_model)){stop("Please specify object containing pre-trained models")}

  if(is.null(pca_in_new_sample)){stop("Please provide blood-based transcriptome PCs for new samples")}

  predict_vals = predict(trained_model, data.frame(pca_in_new_sample))

  if(scale == TRUE){

    predict_vals = scale(predict_vals)
  }

  return(data.frame(predict_vals))
}
