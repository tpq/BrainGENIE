# 3. principal component analysis in reference data
fit_pca = function(gene_list = NULL, autorun = FALSE){

  if(autorun == FALSE & is.null(gene_list) == T){stop("Gene list is expected!")}

  trained_pca_blood = blood_expr

  if(autorun == FALSE){
    if(is.null(gene_list)){stop("Please provide a list of genes that are present in new samples in order to run PCA correctly")}
    matched_genes = intersect(rownames(trained_pca_blood), gene_list)
    if(length(matched_genes) < 1){"No genes common to reference and new samples! Please check that IDs are in ENSEMBL gene ID format."}
    trained_pca_blood = trained_pca_blood[rownames(trained_pca_blood) %in% matched_genes, ]
  }

  return_these_genes = rownames(trained_pca_blood)

  pca_blood = prcomp(t(trained_pca_blood))

  return(list(pca = pca_blood, genes = return_these_genes))
}

# 4. predict PCA in new data
predict_pca = function(dat = NULL, pca_model = NULL, mean_imputation = FALSE){

  if(is.null(dat)){stop("Please supply data frame for new samples")}
  if(is.null(pca_model)){stop("Please supply fitted PCA model")}
  if(class(pca_model) != "list"){stop("Expecting list object for fitted PCA model")}

  new_samples = dat[,colnames(dat) %in% pca_model$genes]

  if(mean_imputation == TRUE){

    means = colMeans(new_samples,na.rm=TRUE)
    for(n in 1:ncol(new_samples)){
      new_samples[is.na(new_samples[,n]),n] = means[[n]]
    }

  }

  if(mean_imputation == FALSE){

    na_detected = colSums(is.na(new_samples))
    na_detected = na_detected[na_detected > 0]
    if(length(na_detected) > 0){stop("Warning! Missing values (NAs) detected in new samples. Switch mean_imputation to TRUE to resolve this issue.")}

  }

  preds = predict(pca_model$pca, newdata = new_samples)
  return(preds)
}
