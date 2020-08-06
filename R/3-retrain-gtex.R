# re-train PCA models using specified gene lists present in target sample:
retrain_gtex = function(gene_list = NULL, output = "", tissue = NULL, ncomps = 20, prop_for_test_set = 0.1, n_folds = 5){

  if(base::exists("brain_expr") == FALSE | base::exists("blood_expr") == FALSE){stop("Warning! Paired blood-brain data have not been imported properly. Please run `load_expr_data` first.")}
  if(is.null(n_folds) == T){stop("Specify the number of cross-validation folds")}
  if(is.null(tissue) == T){stop("Specify the name of brain region for model training")}
  if(output == ""){warning("No output directory specified!")}

  use_n_comps = ncomps

  brain_counts = brain_expr
  blood_counts = blood_expr
  blood_counts = blood_counts[rownames(brain_counts) %in% gene_list, ]
  blood_counts = blood_counts[!rowSums(is.na(blood_counts)), ]

  # optional: define test set
  if(prop_for_test_set > 0){
    test.set = sample(colnames(brain_counts), ncol(brain_counts)*prop_for_test_set)

    test.brain = brain_counts[,colnames(brain_counts) %in% test.set]
    test.blood = blood_counts[,colnames(blood_counts) %in% test.set]

    # remove test set from cross-validation samples
    brain_counts = brain_counts[,!colnames(brain_counts) %in% test.set]
    blood_counts = blood_counts[,!colnames(blood_counts) %in% test.set]
  }

  # set up fold ids for cross-validation
  folds = create.folds(nfolds = n_folds, ids=colnames(brain_counts))

  allmods = list()
  for(fold in 1:max(folds)){

    cat("\rFold: ", fold);cat("\n")

    # define training/testing splits
    x_train = t(blood_counts[,which(folds != fold)])
    y_train = t(brain_counts[,which(folds != fold)])

    x_test = t(blood_counts[,!colnames(blood_counts) %in% rownames(x_train)])
    y_test = t(brain_counts[,!colnames(brain_counts) %in% rownames(y_train)])

    # - Run PCA on blood gene expression data
    pca.blood.train = prcomp(x_train)
    n.blood.comps = pca.blood.train$sdev^2 / sum(pca.blood.train$sdev^2)

    # apply PCA model to test set
    pca.in.test = predict(pca.blood.train, x_test)

    # format PCA predictor matrices based on chosen number of components
    pca.TRAIN = data.frame(pca.blood.train$x[,1:use_n_comps])
    pca.TEST = data.frame(pca.in.test[,1:use_n_comps])


    # apply PCA model to predict gene expression in brain
    fitted.model = lm(as.matrix(y_train)  ~ ., data = pca.TRAIN)
    pred.train = predict(fitted.model)
    train.cors = lapply(1:ncol(pred.train), function(x) cor(pred.train[,x], y_train[,x]))
    names(train.cors) = colnames(pred.train)
    train.cors = ldply(train.cors, .id='gene')
    colnames(train.cors)[2] = 'train_cor'
    train.cors$train_rsq = train.cors$train_cor^2

    # apply PCA model to test set, evaluate accuracy
    pred.test = predict(fitted.model, newdata = pca.TEST)
    test.cors = lapply(1:ncol(pred.test), function(x) cor(pred.test[,x], y_test[,x]))
    names(test.cors) = colnames(pred.test)
    test.cors = ldply(test.cors, .id='gene')
    colnames(test.cors)[2] = 'test_cor'
    test.cors$test_rsq = test.cors$test_cor^2
    test.cors$N = nrow(pca.TEST)
    allmods[[fold]] = data.frame(cbind(train.cors, test.cors[,-1]))
  }

  stats = ldply(allmods)

  # deal with missing values (NAs)
  stats$test_cor[is.na(stats$test_cor)] = 0
  stats$train_cor[is.na(stats$train_cor)] = 0
  stats$test_rsq[is.na(stats$test_rsq)] = 0
  stats$train_rsq[is.na(stats$train_rsq)] = 0

  # transform correlation coefficients
  stats$fisherZ = atanh(stats$test_cor)*sqrt(stats$N-3)
  stats = data.table(stats)

  perf_train = stats[,list(Cor = mean(test_cor, na.rm=TRUE),
                           train.cor = mean(train_cor, na.rm=T),
                           SE = sd(test_cor,  na.rm=TRUE)/sqrt(length(test_cor)),
                           SD = sd(test_cor, na.rm=T),
                           Rsq=mean(test_cor)^2,
                           FisherZ = sum(fisherZ),
                           sample.size = mean(N),
                           Zscore = sum(fisherZ) / sqrt(length(unique(folds)))),by=c("gene")]

  perf_train$pval = 2*pnorm(abs(perf_train$Zscore), lower.tail = FALSE) # calculate two-tailed p-value based on combined z-score
  perf_train$fdr = p.adjust(perf_train$pval, "fdr")


  if(prop_for_test_set > 0){
    # Run final test
    # define training/testing splits
    x_train = t(blood_counts)
    y_train = t(brain_counts)

    x_test = t(test.blood)
    y_test = t(test.brain)

    # - Run PCA on blood gene expression data
    pca.blood.train = prcomp(x_train)
    n.blood.comps = pca.blood.train$sdev^2 / sum(pca.blood.train$sdev^2)
    # use_n_comps = min(which(cumsum(n.blood.comps) >= pca.prop.min))


    # apply PCA model to test set
    pca.in.test = predict(pca.blood.train, x_test)

    # format PCA predictor matrices based on chosen number of components
    pca.TRAIN = data.frame(pca.blood.train$x[,1:use_n_comps])
    pca.TEST = data.frame(pca.in.test[,1:use_n_comps])


    # apply PCA model to predict gene expression in brain
    fitted.model = lm(as.matrix(y_train)  ~ ., data = pca.TRAIN)
    pred.train = predict(fitted.model)
    train.cors = lapply(1:ncol(pred.train), function(x) cor(pred.train[,x], y_train[,x]))
    names(train.cors) = colnames(pred.train)
    train.cors = ldply(train.cors, .id='gene')
    colnames(train.cors)[2] = 'train_cor'
    train.cors$train_rsq = train.cors$train_cor^2

    # apply PCA model to test set, evaluate accuracy
    pred.test = predict(fitted.model, newdata = pca.TEST)
    test.cors = lapply(1:ncol(pred.test), function(x) cor.test(pred.test[,x], y_test[,x]))
    test.cors.est = unlist(lapply(test.cors, function(x) x$estimate))
    test.cors.pval = unlist(lapply(test.cors, function(x) x$p.value))
    names(test.cors.est) = colnames(pred.test)
    test.cors = ldply(test.cors.est, .id='gene')
    colnames(test.cors)[2] = 'test_cor'
    test.cors$test_pval = unlist(test.cors.pval)
    test.cors$test_rsq = test.cors$test_cor^2
    test.cors$N = nrow(pca.TEST)
    final.model = data.frame(cbind(train.cors, test.cors[,-1]))
    colnames(final.model)[-1] = paste("final_", colnames(final.model)[-1], sep="")


    colnames(perf_train)[-c(1)] = paste("cv_", colnames(perf_train)[-c(1)], sep="")
    all.models.merge = merge(perf_train, final.model, by='gene')
    all.models.merge$final_test_fdr = p.adjust(all.models.merge$final_test_pval ,"fdr")

  } else {
    all.models.merge = perf_train
  }

  all.models.merge$tissue = tissue

  # save output as .Rdata file
  saveRDS(all.models.merge, file=paste(output, "/", tissue, "_BrainGENIE_retrain-",ncomps,".Rdata", sep=""))

  perf <<- all.models.merge

}
