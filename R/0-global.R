# 1. Load brain data from GTEx samples
load_expr_data = function(path_to_data = NULL){
  # warning message
  if(is.null(path_to_data) == T){stop("Please provide a directory path")}
  # load data
  data_files = list.files(path_to_data, pattern=".Rdata", full.names=T)
  blood_expr = data.frame(t(readRDS(data_files[grepl("WholeBlood", data_files)])))
  brain_expr = data.frame(t(readRDS(data_files[!grepl("WholeBlood", data_files)])))
  # retain same samples
  common.sids = intersect(colnames(brain_expr), colnames(blood_expr))
  blood_expr <<- blood_expr[,colnames(blood_expr) %in% common.sids]
  brain_expr <<- brain_expr[,colnames(brain_expr) %in% common.sids]
}

# 2. load cross-validation performance
load_cv_performance = function(file_path = NULL){
  if(is.null(file_path)){stop("Please provide full path to .Rdata file with cross-validation accuracies (if relying on pre-trained models")}
  return(data.frame(readRDS(file_path)))
}

# > convert hgnc symbols to ensembl gene ids (gencode v26)
convert_hgnc_to_ensg = function(hgnc = NULL, gtf_path = NULL){

  # check if parameters are filled in
  if(is.null(gtf_path)){stop("Please specify path to the gtf file in the BrainGENIE provided repository")}
  if(is.null(hgnc)){stop("Please specify a vector of HGNC gene symbols")}

  # load GTF file
  gtf = data.frame(fread(gtf_path))
  gtf = gtf[,colnames(gtf) %in% c("gene_id", "gene_name")]
  gtf$gene_name = gsub("[-]", ".", gtf$gene_name)

  hgnc = gsub("[-]", ".", hgnc)

  symbols = data.frame(gene_name = hgnc)
  convert = merge(symbols, gtf, by='gene_name')
  return(convert)
}

convert_ensg_to_hgnc = function(ensg = NULL, gtf_path = NULL){

  # check if parameters are filled in
  if(is.null(gtf_path)){stop("Please specify path to the gtf file in the BrainGENIE provided repository")}
  if(is.null(ensg)){stop("Please specify a vector of ensembl gene ids")}

  # load GTF file
  gtf = data.frame(fread(gtf_path))
  gtf = gtf[,colnames(gtf) %in% c("gene_id", "gene_name")]
  gtf$gene_name = gsub("[-]", ".", gtf$gene_name)

  ensg = gsub("[-]", ".", ensg)

  symbols = data.frame(gene_id = ensg)
  convert = merge(symbols, gtf, by='gene_id')
  return(convert)
}


# create folds for cross-validation
create.folds = function(nfolds = 5, ids = NULL){
  if(is.null(ids)){stop("Expecting IDs column for creating folds")}
  sample(dplyr::ntile(ids, nfolds)) # random folds
}
