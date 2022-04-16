# library(snpnet)
# library(reticulate)
# source("~/TransPRS/R/snpnet_offset.R")

## basic functions
#' read a pfile and save it as a matrix
#'
#' @param new_genotype_file the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param per_idx index of selected persons
#' @param select_SNP index of selected SNPs
#' @importFrom snpnet snpnet
#' @export

read_geno <- function(new_genotype_file, per_idx = NULL, select_SNP = NULL){
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[["zstdcat.path"]], ' ', paste0(new_genotype_file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  per_name <- snpnet:::readIDsFromPsam(paste0(new_genotype_file, '.psam'))
  pvar <- pgenlibr::NewPvar(paste0(new_genotype_file, '.pvar.zst'))
  if (is.null(per_idx)) {
    per_idx <- seq_along(per_name)
  }
  chr <- pgenlibr::NewPgen(paste0(new_genotype_file, '.pgen'), pvar = pvar, sample_subset = per_idx)
  pgenlibr::ClosePvar(pvar)
  if (is.null(select_SNP)) {
    select_SNP = seq_along(vars)
  }
  buf <- pgenlibr::ReadList(chr, select_SNP, meanimpute=F)
  colnames(buf) <- vars[select_SNP]
  rownames(buf) <- per_name[per_idx]
  return(buf)
}

# logit function
logistic<-function(x){
  1/(1+exp(-x))
}

# repeat a vector to generate a matrix
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# repeat a vector to generate a matrix
rep.row<-function(x,n){
  matrix(rep(x,each=n), nrow=n, byrow=F)
}

#' count the lines
#'
#' @param f filenames using 'wc -l filename' command
#'
#' @return
#' @export
wc_R <- function(f){
  str = system(paste0("wc -l ", f), intern = T)
  wc <- as.numeric(strsplit(str, split = " ")[[1]][1])
  return(wc)
}

# set the index
ind.set <- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re <- c(ind.re,1: n.vec[1])
    }else{
      ind.re <- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}

#' Xt%*%beta
#'
#' @param new_genotype_file the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param beta beta value
#' @param p the number of SNPs in pgen file
#' @param nsplit split number of X and beta that we can product part by part to reduce memory use
#'
#' @export
Xbeta_product <- function(new_genotype_file, beta, p, nsplit = NULL){
  if(is.null(dim(beta))){
    beta = matrix(beta, ncol = 1)
  }
  record = 0
  prd = 0
  if(is.null(nsplit)){
    vecSize_seq = p
  }else{
    vecSize_seq = rep(round(p/nsplit), nsplit-1)
    vecSize_seq = c(vecSize_seq, p-sum(vecSize_seq))
  }
  for (size in vecSize_seq) {
    X_part = read_geno(new_genotype_file, select_SNP = c((record+1):(record+size)))
    X_part[is.na(X_part)] = 0
    prd_temp = X_part%*%as.matrix(beta[(record+1):(record+size),])
    prd = prd + prd_temp
    record = record + size
  }
  return(prd)
}


#' only calculate Xt%*%y of appointed samples
#'
#' @param new_genotype_file the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param y phenotypes
#' @param select_per index of person selected
#' @param phenotype the colnames of phenotypes
#'
#' @export
Xy_product <- function(new_genotype_file, y, select_per, phenotype){
  X_part = read_geno(new_genotype_file, per_idx = select_per)
  y = y[select_per, ]
  y_names = c()
  for(i in 1:nrow(y)){
    y_names = c(y_names, paste0(y[,"FID"][i], '_', y[, "IID"][i]))
  }
  rownames(y) = y_names
  y = y[rownames(X_part),phenotype]
  y = as.numeric(y)
  Xty = t(X_part)%*%as.matrix(y)
  return(Xty)
}

## file generation
# transformation of vcf file to pgen file
# id_type = "underline" or "double"
#' transformation of vcf file to pgen file
#'
#' @param out_pgen_name the output PLINK 2.0 pgen file named as out_pgen_name.{pgen,pvar.zst,psam}.
#' @param vcf.file input VCF files
#' @param id_type "underline" or "double", depend on the sample names in VCF files
#'                if the sample name is 'FID_IID', choose id_type="underline"
#'                if the sample name is 'IID' alone , choose id_type="double"
#'
#' @export
pgenGen <- function(out_pgen_name, vcf.file, id_type){
  if(id_type=="underline"){
    system(paste0(configs[['plink2.path']], " --vcf ", vcf.file, " --make-pfile --out ",
                  out_pgen_name, " --id-delim _"))
  }else if(id_type=="double"){
    system(paste0(configs[['plink2.path']], " --vcf ", vcf.file, " --make-pfile --out ",
                  out_pgen_name, "  --double-id"))
  }
  options(warn = -1)
  system(paste0("rm ", out_pgen_name, ".pvar.zst"))
  remove_anno(paste0(out_pgen_name, ".pvar"))
  system(paste0(configs[['zstd.path']], " ", out_pgen_name, ".pvar"))
  system(paste0("rm ", out_pgen_name, ".pvar"))
  system(paste0("rm ", out_pgen_name, ".log"))
}


#' transformation of phenofile with FID and IID to snpnet phenofile with split.col
#'
#' @param ori_phefile original phenofile name. With FID, IID and phenotypes in the file, with no header
#' @param phenotype defile output column names of phenotype columns. "QPHE" default
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param offset a character vector containing the names of the offset included in the lasso
#'               fitting. The names must exist in the column names of the phenotype file.
#' @param split.col defile output column names of split columns. "split" default
#' @param out_phefile output phenofile name
#' @param split_vec proportions to split the data into training set and validating set
#'                  to find the optimal lambda
#'                  c(0.9, 0.1) for "train" and "val" default
#'
#' @export
#'
pheGen <- function(ori_phefile, phenotype = "QPHE", covariates = NULL, offset = NULL, split.col = "split", out_phefile, split_vec = c(0.9, 0.1)){
  phe_out = read.delim(ori_phefile, sep = "", header = F)
  split_as = sample(c("train", "val"), nrow(phe_out), replace = T, p = split_vec)
  phe_out = cbind(phe_out, split_as)
  colnames(phe_out) = c("FID", "IID", phenotype, covariates, offset, split.col)
  write.table(phe_out, sep = ",", out_phefile, row.names = F, quote = F)
}

#' select given sample from pgen file and generate new pgen files
#'
#' @param configs a list of other config parameters.
#' @param genotype.pfile the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param keepnamefile a file with FID and IID of samples to keep names.
#' @param out.pfile the output PLINK 2.0 pgen file name that contains genotype.
#'
#' @export
keep_pfile <- function(configs, genotype.pfile, keepnamefile, out.pfile){
  cmd_keepf <- paste0(configs[['plink2.path']], " --pgen ", genotype.pfile, ".pgen --pvar ",
                      genotype.pfile, ".pvar.zst --psam ", genotype.pfile, ".psam --keep ",  keepnamefile, " --make-pgen --out ", out.pfile)
  system(cmd_keepf, intern = T)
  options(warn = -1)
  system(paste0("rm ", out.pfile, ".pvar.zst"))
  remove_anno(paste0(out.pfile, ".pvar"))
  system(paste0(configs[['zstd.path']], " ", out.pfile, ".pvar"))
}

## evaluation
#' calculate mse
#'
#' @param beta 
#' @param est estimated beta
#' @param X.test X(geno) matrix that can be added to calculate pred.err
#'
#' @export
mse.fun<- function(beta,est, X.test=NULL){
  pred.err<-NA
  est.err<- sum((beta-est)^2)
  if(!is.null(X.test)){
    pred.err <- mean((X.test%*%(beta-est))^2)
  }
  return(list(est.err=est.err, pred.err= pred.err))
}

# rsq calculation
#' return r^2 of prediction (under linear assumption)
#'
#' @param est estimated beta
#' @param X.test X(geno) matrix
#' @param y_real real y (phenotypes)
#'
#' @export
rsq_cal <- function(est, X.test, y_real){
  pred0 = X.test%*%est
  itc = mean(y_real) - mean(pred0)
  pred0 = pred0 + itc
  # rsq <- sum((pred0 - mean(y_real))^2)/sum((y_real-mean(y_real))^2)
  rsq <- (cor(pred0, y_real))^2
  return(rsq)
}


#' find the optimal lambda in linear cases
#'
#' @param genotype.pfile the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table. There should be FID (family ID) and IID (individual ID) columns
#'                       containing the identifier for each individual, and the phenotype column(s).
#'                       (optional) some covariate/offset columns and a column specifying the
#'                       training/validation split can be included in this file.
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param split.col the column name in the phenotype file that specifies the membership of individuals to
#'                  the training or the validation set. The individuals marked as "train" and "val" will
#'                  be treated as the training and validation set, respectively. When specified, the
#'                  model performance is evaluated on both the training and the validation sets.
#' @param offset a character vector containing the names of the offset included in the lasso
#'               fitting. The names must exist in the column names of the phenotype file.
#' @param configs a list of other config parameters.
#' @param nlambda the number of lambda values - default is 50
#' @param lambda.min.ratio smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value,
#'                         i.e. the smallest value for which all coefficients are zero. The default
#'                         depends on the sample size nobs relative to the number of actual variables
#'                         nvars (after QC filtering). The default is 0.01. A very small value of lambda.min.ratio
#'                         will lead to a saturated fit in the nobs < nvars case.
#' @param alpha the elastic-net mixing parameter, where the penalty is defined as
#'              alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2. alpha = 1 corresponds to the lasso penalty,
#'              while alpha = 0 corresponds to the ridge penalty.
#' 
#' @export
cv.snpnet <- function(genotype.pfile = genotype.pfile,
                      phenotype.file = phenotype.file,
                      phenotype = phenotype,
                      covariates = covariates,
                      split.col = NULL,
                      offset = NULL,
                      configs = configs, 
                      nlambda = 50, 
                      lambda.min.ratio = 0.01, 
                      alpha = 1){
  if(is.null(split.col)){
    stop("split.col cannot be NULL.\n")
  }
  fit_snpnet_train <- snpnet_offset(genotype.pfile = genotype.pfile,
                                    phenotype.file = phenotype.file,
                                    phenotype = phenotype,
                                    covariates = covariates,
                                    split.col = split.col,
                                    nlambda = nlambda,
                                    lambda.min.ratio = lambda.min.ratio,
                                    offset = offset,
                                    # mem = 128000, # amount of memory available (MB), recommended
                                    configs = configs, 
                                    alpha = alpha)
  opt_idx <- which.max(fit_snpnet_train$metric.val)
  # metric_optimal_test <- fit_snpnet_refit$metric.val[opt_idx]
  # size_optimal <- sum(fit_snpnet_refit$beta[[opt_idx]] != 0)
  lambda_optimal <- fit_snpnet_train$full.lams[opt_idx]
  beta = fit_snpnet_train$beta[[opt_idx]]
  names = names(fit_snpnet_train$var.rank)
  if(!is.null(covariates)){
    names = c(covariates, names)
  }
  beta = beta[names]
  beta <- ifelse(is.na(beta), 0, beta)
  beta = c(fit_snpnet_train$a0[[opt_idx]], beta)
  return(list(lambda = lambda_optimal, beta = beta))
}

# split the files as training set and testing set, I.til can be non-consecutive integers
#' split the files as training set and testing set, I.til can be non-consecutive integers
#'
#' @param vcf.file geno file in vcf format
#' @param genotype.pfile the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table. There should be FID (family ID) and IID (individual ID) columns
#'                       containing the identifier for each individual, and the phenotype column(s).
#'                       (optional) some covariate/offset columns and a column specifying the
#'                       training/validation split can be included in this file.
#' @param configs a list of other config parameters. 
#' @param get_path dir for saving the intermediate files and result files
#' @param I.til A vector for selected samples to be a testing set. Testing set named as outname[1]
#' @param n total sample size
#' @param outname  output filename: ${outname}.pfile, ${outname}.phe
#'                 c("test", "train") is default
#'
#' @export
split_vcf <- function(vcf.file = vcf.file,
                      genotype.pfile = genotype.pfile,
                      phenotype.file = phenotype.file,
                      configs = configs,
                      get_path = get_path,
                      I.til = I.til, n = sum(n.vec), 
                      outname = c("test", "train")){
  remove_anno(vcf.file)
  system(paste0("grep \'^#CHROM\' ", vcf.file, " | sed \'s/\\t/\\n/g\' | sed \'s/_/ /g\' > ", get_path ,"/all_name.txt"))
  ## extract selected samples
  all_name <- read.delim(paste0(get_path ,"/all_name.txt"), header=TRUE)
  all_name <- as.vector(all_name[-1:-8, 1])
  I.el = 1:n
  I.el = I.el[-match(I.til,I.el)]
  X1 = all_name[I.til]
  X2 = all_name[I.el]
  write.table(cbind(X1, X1), paste0(get_path ,"/X_", outname[1], ".txt"), sep = " ", quote = F,col.names = F, row.names = F)
  write.table(cbind(X2, X2), paste0(get_path ,"/X_", outname[2], ".txt"), sep = " ", quote = F,col.names = F, row.names = F)
  y = read.csv(phenotype.file, header = T)
  y1 = y[I.til, ]
  y2 = y[I.el, ]
  write.csv(y1, paste0(get_path ,"/y_", outname[1], ".phe"), quote = F, row.names = F)
  write.csv(y2, paste0(get_path ,"/y_", outname[2], ".phe"), quote = F, row.names = F)
  system(paste0(configs[['plink2.path']], " --pgen ", genotype.pfile, ".pgen --pvar ",
                genotype.pfile, ".pvar.zst --psam ", genotype.pfile, ".psam --keep ",
                get_path ,"/X_", outname[1], ".txt --make-pgen --out ", get_path ,"/X_", outname[1]))
  system(paste0(configs[['plink2.path']], " --pgen ", genotype.pfile, ".pgen --pvar ",
                genotype.pfile, ".pvar.zst --psam ", genotype.pfile, ".psam --keep ",
                get_path ,"/X_", outname[2], ".txt --make-pgen --out ", get_path ,"/X_", outname[2]))
  options(warn = -1)
  system(paste0("rm ", get_path ,"/X_", outname[1], ".pvar.zst"))
  system(paste0("rm ", get_path ,"/X_", outname[2], ".pvar.zst"))
  system(paste0(configs[['zstd.path']], " ", get_path ,"/X_", outname[1], ".pvar"))
  system(paste0(configs[['zstd.path']], " ", get_path ,"/X_", outname[2], ".pvar"))
  system(paste0(configs[['plink2.path']], " --pgen ", get_path ,"/X_", outname[1], ".pgen --pvar ", get_path ,"/X_", outname[1], ".pvar.zst --psam ", get_path ,"/X_", outname[1], ".psam --export vcf --out ", get_path ,"/X_", outname[1]))
  system(paste0(configs[['plink2.path']], " --pgen ", get_path ,"/X_", outname[2], ".pgen --pvar ", get_path ,"/X_", outname[2], ".pvar.zst --psam ", get_path ,"/X_", outname[2], ".psam --export vcf --out ", get_path ,"/X_", outname[2]))
  system(paste0("rm ", get_path ,"/X_", outname[1], ".txt"))
  system(paste0("rm ", get_path ,"/X_", outname[2], ".txt"))
  system(paste0("rm ", get_path ,"/*.log"))
}


#' remove all the annotations that start with '##'
#'
#' @param filename 
#'
#' @export
remove_anno <- function(filename){
  sed_command <- paste0("sed -i '/^##/d' ", filename)
  system(sed_command)
}


#' bind vcf using python
#'
#' @param nameLis We assume the existence of nameLis.{vcf}.
#' @param out output vcf file name
#' @importFrom reticulate source_python
#' @export
bindVCF_py <- function(nameLis, out){
  if(length(nameLis) == 1){
    system(paste0("cp ", nameLis, " ", out))
  }else{
    source_python(system.file('python', 'bind_vcf.py', package = 'TransPRS'))
    # source_python("../src/bind_vcf.py")
    bind_vcf(nameLis, out)
  }
}

#' bind phe
#'
#' @param nameLis nameLis of phenofile. We assume the existence of nameLis.{phe}
#' @param out output phenofile name
#' @export
bind_phe <- function(nameLis, out){
  system(paste0("head -1 ", nameLis[1], " > ", out))
  for(nm in nameLis){
    system(paste0("cat ", nm, " | sed '1d' >> ", out))
  }
}

# Aggregation rule
#' linear aggregation rule
#'
#' @param B a matrix for column-binded possible beta
#' @param Xtest.pfile genofile for test. The PLINK 2.0 pgen file that contains genotype.
#'                    We assume the existence of Xtest.pfile.{pgen,pvar.zst,psam}.
#' @param y.test phenofile for test
#' @param nsplit split number of X and beta that we can product part by part to reduce memory use. default is 4 
#' @param const 
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#'
#' @export
agg.fun<- function(B, Xtest.pfile, y.test, nsplit = 4, const=2, phenotype, covariates = NULL){
  p<-nrow(B)-1-length(covariates)
  K<-ncol(B)
  colnames(B)<-NULL
  if (is.null(covariates)){
    Xtilb <- Xbeta_product(Xtest.pfile, as.matrix(B)[-1, ], p, nsplit)+rep.row(as.matrix(B)[1, ], length(y.test[, phenotype]))
  }else{
    Xtilb <- Xbeta_product(Xtest.pfile, as.matrix(B)[-(1:(1+length(covariates))), ], p, nsplit)+rep.row(as.matrix(B)[1, ], length(y.test[, phenotype])) + as.matrix(y.test[, covariates]) %*% as.matrix(B)[2:(1+length(covariates)),]
  }
  loss.B <- apply(Xtilb, 2, function(bb){mean((y.test[, phenotype]-bb)^2)})
  eta.hat <- exp(-const*loss.B)/sum(exp(-const*(loss.B[which(!is.nan(loss.B))])))
  eta.hat <- ifelse(is.nan(eta.hat), 0, eta.hat)
  while(sum(eta.hat) == 0){
    const = const/2
    eta.hat <- exp(-const*loss.B)/sum(exp(-const*(loss.B[which(!is.nan(loss.B))])))
    eta.hat <- ifelse(is.nan(eta.hat), 0, eta.hat)
  }
  eta.hat
  # if(sum(B==0)==ncol(B)*nrow(B)){
  #   return(rep(0,nrow(B)))
  # }
  # p<-nrow(B)-1-length(covariates)
  # K<-ncol(B)
  # colnames(B)<-NULL
  # if (is.null(covariates)){
  #   temp_delta = (rep.col(y.test[, phenotype], K)-(Xbeta_product(Xtest.pfile, as.matrix(B)[-1, ], p, nsplit)+rep.row(as.matrix(B)[1, ], length(y.test[, phenotype]))))^2
  # }else{
  #   temp_delta = (rep.col(y.test[, phenotype], K)-(Xbeta_product(Xtest.pfile, as.matrix(B)[-(1:(1+length(covariates))), ], p, nsplit)+rep.row(as.matrix(B)[1, ], length(y.test[, phenotype]))))^2 + as.matrix(y.test[, covariates]) %*% as.matrix(B)[2:(1+length(covariates)),]
  # }
  # temp_delta <- ifelse(is.na(temp_delta), max(temp_delta), temp_delta)
  # if(selection){#select beta.hat with smallest prediction error
  #   khat<-which.min(colSums(temp_delta))
  #   theta.hat<-rep(0, ncol(B))
  #   theta.hat[khat] <- 1
  #   beta=as.matrix(B)[,khat]
  #   beta.ew=NULL
  # }else{#Q-aggregation
  #   theta.hat<- exp(-colSums(temp_delta)/2)
  #   if(max(theta.hat) == 0){
  #     khat<-which.min(colSums(temp_delta))
  #     theta.hat<-rep(0, ncol(B))
  #     theta.hat[khat] <- 1
  #     beta=as.matrix(B)[,khat]
  #     beta.ew=NULL
  #     return(list(theta=theta.hat, beta=beta, beta.ew=beta.ew))
  #   }
  #   theta.hat <- ifelse(is.na(theta.hat), 0, theta.hat)
  #   theta.hat=theta.hat/sum(theta.hat)
  #   theta.old=theta.hat
  #   beta<-as.numeric(as.matrix(B)%*%theta.hat)
  #   beta.ew<-beta
  #   # theta.old=theta.hat
  #   for(ss in 1:total.step){
  #     if (is.null(covariates)){
  #       theta.hat<- exp(-colSums(temp_delta)/2+colSums((rep.col((Xbeta_product(Xtest.pfile, as.matrix(beta)[-1, ], p, nsplit)+rep.row(as.matrix(beta)[1, ], nrow(temp_delta))), K)-(Xbeta_product(Xtest.pfile, as.matrix(as.matrix(B)[-1, ]),p, nsplit)+rep.row(as.matrix(B)[1, ], nrow(temp_delta))))^2)/8)
  #     }else{
  #       theta.hat<- exp(-colSums(temp_delta)/2+colSums((rep.col((Xbeta_product(Xtest.pfile, as.matrix(beta)[-(1:(1+length(covariates))), ], p, nsplit)+rep.row(as.matrix(beta)[1, ], nrow(temp_delta))+(as.matrix(y.test[, covariates]) %*% as.matrix(beta)[2:(1+length(covariates)),])), K)-(Xbeta_product(Xtest.pfile, as.matrix(as.matrix(B)[-(1:(1+length(covariates))), ]),p, nsplit)+rep.row(as.matrix(B)[1, ], nrow(temp_delta))+(as.matrix(y.test[, covariates]) %*% as.matrix(B)[2:(1+length(covariates)),])))^2)/8)
  #     }
  #     theta.hat<-theta.hat/sum(theta.hat)
  #     if(is.na(sum(theta.hat))){
  #       theta.hat = theta.old
  #       break
  #     }
  #     beta<- as.numeric(B%*%theta.hat*1/4+3/4*beta)
  #     dtheta = sum(abs(theta.hat-theta.old))
  #     if (is.na(dtheta)) {
  #       theta.hat = theta.old
  #       break
  #     }
  #     if(dtheta<10^(-3)){break}
  #     theta.old=theta.hat
  #   }
  # }
  # list(theta=theta.hat, beta=beta, beta.ew=beta.ew)
}

#' binomial aggregation rule
#'
#' @param B a matrix for column-binded possible beta
#' @param Xtest.pfile genofile for test. The PLINK 2.0 pgen file that contains genotype.
#'                    We assume the existence of Xtest.pfile.{pgen,pvar.zst,psam}.
#' @param y.test phenofile for test
#' @param nsplit split number of X and beta that we can product part by part to reduce memory use. default is 4 
#' @param const 
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#'
#' @export
Binary_agg.fun<-function(B, Xtest.pfile, y.test, nsplit = 4, const=2, phenotype, covariates = NULL){
  p<-nrow(B)-1-length(covariates)
  K<-ncol(B)
  colnames(B)<-NULL
  if (is.null(covariates)){
    Xtilb <- Xbeta_product(Xtest.pfile, as.matrix(B)[-1, ], p, nsplit)+rep.row(as.matrix(B)[1, ], length(y.test[, phenotype]))
  }else{
    Xtilb <- Xbeta_product(Xtest.pfile, as.matrix(B)[-(1:(1+length(covariates))), ], p, nsplit)+rep.row(as.matrix(B)[1, ], length(y.test[, phenotype])) + as.matrix(y.test[, covariates]) %*% as.matrix(B)[2:(1+length(covariates)),]
  }
  loss.B <- apply(Xtilb, 2, function(bb){-sum(y.test[, phenotype]*log(logistic(bb))+(1-y.test[, phenotype])*log(1-logistic(bb)))})
  eta.hat <- exp(-const*loss.B)/sum(exp(-const*(loss.B[which(!is.nan(loss.B))])))
  eta.hat <- ifelse(is.nan(eta.hat), 0, eta.hat)
  while(sum(eta.hat) == 0){
    const = const/2
    eta.hat <- exp(-const*loss.B)/sum(exp(-const*(loss.B[which(!is.nan(loss.B))])))
    eta.hat <- ifelse(is.nan(eta.hat), 0, eta.hat)
  }
  eta.hat
}


#' snpnet in binary cases
#'
#' @param genotype.pfile the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table. There should be FID (family ID) and IID (individual ID) columns
#'                       containing the identifier for each individual, and the phenotype column(s).
#'                       (optional) some covariate/offset columns and a column specifying the
#'                       training/validation split can be included in this file.
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param split.col the column name in the phenotype file that specifies the membership of individuals to
#'                  the training or the validation set. The individuals marked as "train" and "val" will
#'                  be treated as the training and validation set, respectively. When specified, the
#'                  model performance is evaluated on both the training and the validation sets.
#' @param offset a character vector containing the names of the offset included in the lasso
#'               fitting. The names must exist in the column names of the phenotype file.
#' @param configs a list of other config parameters.
#' @param nlambda the number of lambda values - default is 50
#' @param lambda.min.ratio smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value,
#'                         i.e. the smallest value for which all coefficients are zero. The default
#'                         depends on the sample size nobs relative to the number of actual variables
#'                         nvars (after QC filtering). The default is 0.01. A very small value of lambda.min.ratio
#'                         will lead to a saturated fit in the nobs < nvars case.
#' @param lambda A user supplied \code{lambda} sequence.
#' @param alpha the elastic-net mixing parameter, where the penalty is defined as
#'              alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2. alpha = 1 corresponds to the lasso penalty,
#'              while alpha = 0 corresponds to the ridge penalty.
#' @export
cv.binary_lasso <- function(genotype.pfile = genotype.pfile,
                            phenotype.file = phenotype.file,
                            phenotype = phenotype,
                            covariates = NULL,
                            offset = NULL,
                            configs = configs, 
                            split.col = NULL,
                            nlambda = 30, 
                            lambda.min.ratio = 0.01, 
                            lambda = NULL, 
                            alpha = 1){
  if (!is.null(lambda)){
    nlambda <- length(lambda)
    split.col = NULL
  }
  fit_snpnet_train <- snpnet_offset(genotype.pfile = genotype.pfile,
                                    phenotype.file = phenotype.file,
                                    phenotype = phenotype,
                                    covariates = covariates,
                                    offset = offset,
                                    split.col = split.col,
                                    family = "binomial",
                                    nlambda = nlambda,
                                    lambda.min.ratio = lambda.min.ratio, 
                                    lambda = lambda,
                                    # mem = 128000, # amount of memory available (MB), recommended
                                    configs = configs, 
                                    alpha = alpha) 
  if(is.null(split.col)){
    beta = fit_snpnet_train$beta[[1]]
    lambda_optimal = lambda
    a0 = fit_snpnet_train$a0[[1]]
  }else{
    opt_idx <- which.max(fit_snpnet_train$metric.val)
    # metric_optimal_test <- fit_snpnet_refit$metric.val[opt_idx]
    # size_optimal <- sum(fit_snpnet_refit$beta[[opt_idx]] != 0)
    lambda_optimal <- fit_snpnet_train$full.lams[opt_idx]
    beta = fit_snpnet_train$beta[[opt_idx]]
    a0 = fit_snpnet_train$a0[[opt_idx]]
  }
  names = names(fit_snpnet_train$var.rank)
  if(!is.null(covariates)){
    names = c(covariates, names)
  }
  beta <- beta[names]
  names(beta) <- names
  beta <- ifelse(is.na(beta), 0, beta)
  beta = c(a0, beta)
  return(list(lambda = lambda_optimal, beta = beta))
}

# ST.init<-function(X.tar,y.tar){ #single-task initialization
#' single-task initialization
#'
#' @param genotype.pfile the PLINK 2.0 pgen file that contains genotype.
#'                       We assume the existence of genotype.pfile.{pgen,pvar.zst,psam}.
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table. There should be FID (family ID) and IID (individual ID) columns
#'                       containing the identifier for each individual, and the phenotype column(s).
#'                       (optional) some covariate/offset columns and a column specifying the
#'                       training/validation split can be included in this file.
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param offset  a character vector containing the names of the offset included in the lasso
#'               fitting. The names must exist in the column names of the phenotype file.
#' @param split.col the column name in the phenotype file that specifies the membership of individuals to
#'                  the training or the validation set. The individuals marked as "train" and "val" will
#'                  be treated as the training and validation set, respectively. When specified, the
#'                  model performance is evaluated on both the training and the validation sets.
#' @param configs a list of other config parameters.
#' @param nlambda the number of lambda values - default is 50
#' @param lambda.min.ratio smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value,
#'                         i.e. the smallest value for which all coefficients are zero. The default
#'                         depends on the sample size nobs relative to the number of actual variables
#'                         nvars (after QC filtering). The default is 0.01. A very small value of lambda.min.ratio
#'                         will lead to a saturated fit in the nobs < nvars case.
#' @param p number of SNPs in genofile
#' @param family the type of the phenotype: "gaussian", "binomial".
#' @param alpha the elastic-net mixing parameter, where the penalty is defined as
#'              alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2. alpha = 1 corresponds to the lasso penalty,
#'              while alpha = 0 corresponds to the ridge penalty.
#' @param shortcut whether to avoid using tuned lambda to caltulate the parameters on all samples again. If the
#'                 sample size is large enough, let shortcut = T can save computational costs. #' 
#' @export
#'
ST.init<-function(genotype.pfile = genotype.pfile, 
                  phenotype.file = phenotype.file, 
                  phenotype = phenotype, 
                  covariates = NULL, 
                  offset = NULL, 
                  split.col = NULL, 
                  configs = configs, 
                  nlambda = 30, 
                  lambda.min.ratio = 0.01, 
                  p = NULL,
                  family, 
                  alpha = 1, 
                  shortcut = T){
  if(is.null(p)){
    p <- length(dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[['zstdcat.path']],
                                                                         ' ', paste0(genotype.pfile,
                                                                                     '.pvar.zst'))),
                                            'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID)
  }

  if(is.null(split.col)){
    stop("split.col cannot be NULL.\n")
  }

  n0.tar <- wc_R(phenotype.file) - 1
  # fit0.tar<-cv.glmnet(x=X.tar,y=y.tar, lambda=seq(0.25,0.05, length.out=20)*sqrt(2*log(p)/n0.tar),family='binomial')
  if(family == "binomial"){
    fit0.tar <- cv.binary_lasso(genotype.pfile = genotype.pfile,
                                phenotype.file = phenotype.file,
                                phenotype = phenotype,
                                covariates = covariates,
                                offset = offset, 
                                split.col = split.col,
                                configs = configs, 
                                nlambda = nlambda, 
                                lambda.min.ratio = lambda.min.ratio, 
                                alpha = alpha)
    
    lam0.tar<-fit0.tar$lambda
    lam.const=lam0.tar/sqrt(2*log(p)/n0.tar)
    if(shortcut == F){
      fit0.tar <- cv.binary_lasso(genotype.pfile = genotype.pfile,
                                  phenotype.file = phenotype.file,
                                  phenotype = phenotype,
                                  covariates = covariates,
                                  offset = offset, 
                                  split.col = NULL,
                                  configs = configs, 
                                  lambda = lam0.tar, 
                                  alpha = alpha)
    }
    beta0 <- fit0.tar$beta
  }else if(family == "gaussian") {
    cv_sn <- cv.snpnet(genotype.pfile = genotype.pfile,
                         phenotype.file = phenotype.file,
                         phenotype = phenotype,
                         covariates = covariates,
                         split.col = split.col,
                         nlambda = nlambda, 
                         lambda.min.ratio = lambda.min.ratio,
                         offset = offset, 
                         configs = configs, 
                         alpha = alpha)
    opt_lam <- cv_sn$lambda
    lam.const=opt_lam/sqrt(2*log(p)/n0.tar)
    if(shortcut == F){
      fit_snpnet <- snpnet_offset(genotype.pfile = genotype.pfile,
                            phenotype.file = phenotype.file,
                            phenotype = phenotype,
                            covariates = covariates,
                            lambda = opt_lam,
                            offset = offset, 
                            configs = configs, 
                            alpha = alpha)
      beta0 <- fit_snpnet$beta[[1]]
      beta0<-beta0*(abs(beta0)>=opt_lam)
      names = names(fit_snpnet$var.rank)
      if(!is.null(covariates)){
        names = c(covariates, names)
      }
      beta0 <- beta0[names]
      names(beta0) <- names
      beta0 <- ifelse(is.na(beta0), 0, beta0)
      beta0 = c(fit_snpnet$a0[[1]], beta0)
    }else{
      beta0 <- cv_sn$beta
    }
  }
  return(list(beta0=beta0, lam.const=lam.const))
}

#' TL.init
#'
#' @param X.tarLis X.tarLis = paste0(XXX, ".vcf")
#' @param y.tarLis y.tarLis = paste0(XXX, ".phe")
#' @param X.srcLis X.srcLis[[i]] = paste0(XXX, ".vcf")
#' @param y.srcLis X.srcLis[[i]] = paste0(XXX, ".vcf")
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file. 
#' @param split.col the column name in the phenotype file that specifies the membership of individuals to
#'                  the training or the validation set. The individuals marked as "train" and "val" will
#'                  be treated as the training and validation set, respectively. When specified, the
#'                  model performance is evaluated on both the training and the validation sets.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param get_path  dir for saving the intermediate files and result files
#' @param configs  a list of other config parameters.
#' @param nlambda the number of lambda values - default is 50
#' @param lambda.min.ratio smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value,
#'                         i.e. the smallest value for which all coefficients are zero. The default
#'                         depends on the sample size nobs relative to the number of actual variables
#'                         nvars (after QC filtering). The default is 0.01. A very small value of lambda.min.ratio
#'                         will lead to a saturated fit in the nobs < nvars case.
#' @param p number of SNPs in genofile
#' @param family the type of the phenotype: "gaussian", "binomial".
#' @param id_type "underline" or "double", depend on the sample names in VCF files
#'                if the sample name is 'FID_IID', choose id_type="underline"
#'                if the sample name is 'IID' alone , choose id_type="double"
#' @param rm whether to remove intermediate files after the end of the functions. default is F
#' @param split_for_agg whether to split the primary sample files for later aggregation. default is F
#' @param ValProp the proportion of validating set for aggregation when spliting the primary sample files for later aggregation. default is 1/10
#' @param nsplit split number of X and beta that we can product part by part to reduce memory use. default is 4 
#' @param alpha the elastic-net mixing parameter, where the penalty is defined as
#'              alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2. alpha = 1 corresponds to the lasso penalty,
#'              while alpha = 0 corresponds to the ridge penalty.
#' @param shortcut whether to avoid using tuned lambda to caltulate the parameters on all samples again. If the
#'                 sample size is large enough, let shortcut = T can save computational costs.
#' @export
TL.init <- function(X.tarLis, y.tarLis, X.srcLis, y.srcLis,
                    phenotype = phenotype, split.col = NULL,
                    covariates = covariates, get_path = get_path,
                    configs = configs, nlambda = 30, lambda.min.ratio = 0.01, p = NULL,
                    family, id_type = "double", rm = F, split_for_agg = F, ValProp = 1/10, nsplit = 4, alpha = 1, shortcut = T){
  lam.const = list()
  w0 = list()
  delta0 = list()
  if(is.null(split.col)){
    stop("split.col cannot be NULL.\n")
  }
  for(i in 1:length(X.srcLis)){
    bindVCF_py(X.srcLis[[i]], paste0(get_path , "/TLinit_src_", i, ".vcf"))
    pgenGen(paste0(get_path , "/TLinit_src_", i), paste0(get_path , "/TLinit_src_", i, ".vcf"), id_type = id_type)
    bind_phe(y.srcLis[[i]], paste0(get_path ,"/TLinit_src_", i, ".phe"))
    if(family == "binomial"){
      fit.src <- cv.binary_lasso(genotype.pfile = paste0(get_path , "/TLinit_src_", i),
                                 phenotype.file = paste0(get_path ,"/TLinit_src_", i, ".phe"),
                                 phenotype = phenotype,
                                 covariates = covariates,
                                 offset = NULL, 
                                 split.col = split.col,
                                 configs = configs, 
                                 nlambda = nlambda, 
                                 lambda.min.ratio = lambda.min.ratio, 
                                 alpha = alpha)
      lam0.src <- fit.src$lambda
      lam.const[[i]]=lam0.src/sqrt(2*log(p)/(wc_R(paste0(get_path ,"/TLinit_src_", i, ".phe"))-1))
      if(shortcut == F){
        fit.src <- cv.binary_lasso(genotype.pfile = paste0(get_path , "/TLinit_src_", i),
                                  phenotype.file = paste0(get_path ,"/TLinit_src_", i, ".phe"),
                                  phenotype = phenotype,
                                  covariates = covariates,
                                  offset = NULL, 
                                  split.col = NULL,
                                  configs = configs, 
                                  lambda = lam0.src, 
                                  alpha = alpha)
      }
      w0[[i]] = fit.src$beta
    }else if(family == "gaussian"){
      cv_sn <- cv.snpnet(genotype.pfile = paste0(get_path , "/TLinit_src_", i),
                        phenotype.file = paste0(get_path ,"/TLinit_src_", i, ".phe"),
                        phenotype = phenotype,
                        split.col = split.col,
                        covariates = covariates,
                        nlambda = nlambda, 
                        lambda.min.ratio = lambda.min.ratio,
                        configs = configs, 
                        alpha = alpha)
      opt_lam <- cv_sn$lambda
      lam.const[[i]]=opt_lam/sqrt(2*log(p)/(wc_R(paste0(get_path ,"/TLinit_src_", i, ".phe"))-1))
      if(shortcut == F){
        fit_snpnet <- snpnet(genotype.pfile = paste0(get_path , "/TLinit_src_", i),
                            phenotype.file = paste0(get_path ,"/TLinit_src_", i, ".phe"),
                            phenotype = phenotype,
                            covariates = covariates,
                            lambda = opt_lam,
                            split.col = NULL, # split column name in phenotype.file with train/val/test labels
                            # mem = 128000, # amount of memory available (MB), recommended
                            configs = configs, 
                            alpha = alpha)
        w0[[i]] <- fit_snpnet$beta[[1]]
        w0[[i]] <- w0[[i]]*(abs(w0[[i]])>=opt_lam)
        if (sum(w0[[i]]) == 0) {
           w0[[i]] <- cv_sn$beta
        }else{
          names = names(fit_snpnet$var.rank)
          if(!is.null(covariates)){
            names = c(covariates, names)
          }
          w0[[i]] <- w0[[i]][names]
          names(w0[[i]]) <- names
          w0[[i]] <- ifelse(is.na(w0[[i]]), 0, w0[[i]])
          w0[[i]] = c(fit_snpnet$a0[[1]], w0[[i]]) 
        } 
      }else{
        w0[[i]] <- cv_sn$beta
      }
    }
  }
  bindVCF_py(X.tarLis, paste0(get_path ,"/TLinit_tar.vcf"))
  pgenGen(paste0(get_path ,"/TLinit_tar"), paste0(get_path ,"/TLinit_tar.vcf"), id_type = id_type)
  bind_phe(y.tarLis, paste0(get_path ,"/TLinit_tar.phe"))
  if(split_for_agg){
    all_phe = read.csv(paste0(get_path ,"/TLinit_tar.phe"))
    idx_val = sort(sample(1:nrow(all_phe), ceiling(nrow(all_phe)*ValProp), replace = F))
    write.csv(all_phe[idx_val, ], paste0(get_path ,"/TLinit_val.phe"), row.names = F, quote = F)
    write.csv(all_phe[-idx_val, ], paste0(get_path ,"/TLinit_tar.phe"), row.names = F, quote = F)
    if (id_type == "underline"){
      write.table(all_phe[idx_val, 1:2], paste0(get_path ,"/keep_val.txt"), sep = ' ', col.names = F, row.names = F, quote = F)
      write.table(all_phe[-idx_val, 1:2], paste0(get_path ,"/keep_tar.txt"), sep = ' ', col.names = F, row.names = F, quote = F)
      system(paste0(configs[['plink2.path']], " --vcf ", get_path ,"/TLinit_tar.vcf --keep ",
                    get_path ,"/keep_val.txt --export vcf --out ", get_path ,"/TLinit_val --id-delim _"))
      remove_anno(paste0(get_path ,"/TLinit_val.vcf"))
      pgenGen(paste0(get_path ,"/TLinit_val"), paste0(get_path ,"/TLinit_val.vcf"), id_type = id_type)
      system(paste0(configs[['plink2.path']], " --vcf ", get_path ,"/TLinit_tar.vcf --keep ",
                    get_path ,"/keep_tar.txt --export vcf --out ", get_path ,"/TLinit_tar --id-delim _"))
      remove_anno(paste0(get_path ,"/TLinit_tar.vcf"))
      pgenGen(paste0(get_path ,"/TLinit_tar"), paste0(get_path ,"/TLinit_tar.vcf"), id_type = id_type)
    }else{
      write.table(all_phe[idx_val, 1], paste0(get_path ,"/keep_val.txt"), sep = ' ', col.names = F, row.names = F, quote = F)
      write.table(all_phe[-idx_val, 1], paste0(get_path ,"/keep_tar.txt"), sep = ' ', col.names = F, row.names = F, quote = F)
      system(paste0(configs[['plink2.path']], " --vcf ", get_path ,"/TLinit_tar.vcf --keep ",
                    get_path ,"/keep_val.txt --export vcf --out ", get_path ,"/TLinit_val"))
      remove_anno(paste0(get_path ,"/TLinit_val.vcf"))
      pgenGen(paste0(get_path ,"/TLinit_val"), paste0(get_path ,"/TLinit_val.vcf"), id_type = id_type)
      system(paste0(configs[['plink2.path']], " --vcf ", get_path ,"/TLinit_tar.vcf --keep ",
                    get_path ,"/keep_tar.txt --export vcf --out ", get_path ,"/TLinit_tar"))
      remove_anno(paste0(get_path ,"/TLinit_tar.vcf"))
      pgenGen(paste0(get_path ,"/TLinit_tar"), paste0(get_path ,"/TLinit_tar.vcf"), id_type = id_type)
    }
  }
  tar_phe = read.csv(paste0(get_path ,"/TLinit_tar.phe"), header = T)
  for(i in 1:length(X.srcLis)){
    if (is.null(covariates)){
      Xb = Xbeta_product(paste0(get_path ,"/TLinit_tar"), w0[[i]][-1], p, nsplit) + w0[[i]][1]
    }else{
      Xb = Xbeta_product(paste0(get_path ,"/TLinit_tar"), w0[[i]][-(1:(1+length(covariates)))], p, nsplit) + w0[[i]][1] + as.matrix(tar_phe[, covariates]) %*% w0[[i]][2:(1+length(covariates))]
    }
    tar_phe = cbind(tar_phe, Xb)
    colnames(tar_phe)[ncol(tar_phe)] = paste0("Xb", i)
  }
  write.csv(tar_phe, paste0(get_path, "/TLinit_tarXb.phe"), row.names = F, quote = F)
  bindVCF_py(c(paste0(get_path ,"/TLinit_tar.vcf"), unlist(X.srcLis)), paste0(get_path, "/global.vcf"))
  pgenGen(paste0(get_path, "/global"), paste0(get_path, "/global.vcf"), id_type = id_type)
  bind_phe(c(paste0(get_path ,"/TLinit_tar.phe"), unlist(y.srcLis)), paste0(get_path, "/global.phe"))
  
  if(family == "binomial"){
    for(i in 1:length(X.srcLis)){
      delta_fit <- cv.binary_lasso(genotype.pfile = paste0(get_path ,"/TLinit_tar"),
                                   phenotype.file = paste0(get_path ,"/TLinit_tarXb.phe"),
                                   phenotype = phenotype,
                                   covariates = covariates,
                                   offset = paste0("Xb", i),
                                   configs = configs, 
                                   split.col = split.col,
                                   nlambda = nlambda, 
                                   lambda.min.ratio = lambda.min.ratio, 
                                   alpha = alpha)
      lam0.tar <- delta_fit$lambda
      # lam.const[[i]]=lam0.src/sqrt(2*log(p)/(wc_R(paste0(get_path ,"/TLinit_tarXb.phe"))-1))
      if(shortcut == F){
        fit.delta <- cv.binary_lasso(genotype.pfile = paste0(get_path ,"/TLinit_tar"),
                                    phenotype.file = paste0(get_path ,"/TLinit_tarXb.phe"),
                                    phenotype = phenotype,
                                    covariates = covariates,
                                    offset = paste0("Xb", i),
                                    configs = configs, 
                                    split.col = NULL,
                                    lambda = lam0.tar, 
                                    alpha = alpha)
        delta0[[i]] <- fit.delta$beta
      }else{
        delta0[[i]] <- delta_fit$beta
      }
    }
  }else if(family == "gaussian"){
    for(i in 1:length(X.srcLis)){
      cv_sn <- cv.snpnet(genotype.pfile = paste0(get_path ,"/TLinit_tar"),
                             phenotype.file = paste0(get_path ,"/TLinit_tarXb.phe"),
                             phenotype = phenotype,
                             covariates = covariates,
                             offset = paste0("Xb", i), 
                             configs = configs, 
                             split.col = split.col,
                             nlambda = nlambda, 
                             lambda.min.ratio = lambda.min.ratio, 
                             alpha = alpha)
      lam.del <- cv_sn$lambda
      if(shortcut == F){
        fit_delta <- snpnet_offset(genotype.pfile = paste0(get_path ,"/TLinit_tar"),
                                  phenotype.file = paste0(get_path ,"/TLinit_tarXb.phe"),
                                  phenotype = phenotype,
                                  covariates = covariates,
                                  offset = paste0("Xb", i), 
                                  configs = configs, 
                                  lambda = lam.del,
                                  split.col = NULL, 
                                  alpha = alpha)
        delta0[[i]] = fit_delta$beta[[1]]
        names = names(fit_delta$var.rank)
        if(!is.null(covariates)){
          names = c(covariates, names)
        }
        delta0[[i]] <- delta0[[i]][names]
        names(delta0[[i]]) <- names
        delta0[[i]] <- ifelse(is.na(delta0[[i]]), 0, delta0[[i]])
        delta0[[i]] = c(fit_delta$a0[[1]], delta0[[i]])
      }else{
        delta0[[i]] <- cv_sn$beta
      }

    }
  }
  
  ## pooled PRS
  Z = c(rep(0, wc_R(paste0(get_path , "/TLinit_tar.phe"))-1))
  for(i in 1:length(X.srcLis)){
    if (is.null(covariates)){
      Z = c(Z, -(Xbeta_product(paste0(get_path ,"/TLinit_src_", i), delta0[[i]][-1], p, nsplit))-delta0[[i]][1])
    }else{
      src_phe = read.csv(paste0(get_path ,"/TLinit_src_", i, ".phe"), header = T)
      Z = c(Z, -(Xbeta_product(paste0(get_path ,"/TLinit_src_", i), delta0[[i]][-(1:(1+length(covariates)))], p, nsplit))-delta0[[i]][1]-(as.matrix(src_phe[, covariates]) %*% delta0[[i]][2:(1+length(covariates))]))
    }
  }
  
  tar_phe = read.csv(paste0(get_path ,"/global.phe"), header = T)
  tar_phe$Z = Z
  write.csv(tar_phe, paste0(get_path, "/global_Z.phe"), row.names = F, quote = F)
  beta0 <- ST.init(genotype.pfile = paste0(get_path , "/global"), 
                   phenotype.file = paste0(get_path , "/global_Z.phe"), phenotype = phenotype, split.col = split.col,
                   covariates = covariates, offset = "Z", configs = configs, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, family = family, alpha = alpha, shortcut = shortcut)$beta0
  
  if(rm == T){
    system(paste0("rm ", get_path, "/TLinit_tarXb.phe"))
    system(paste0("rm ", get_path ,"/TLinit_tar*"))
    system(paste0("rm ", get_path ,"/global.*"))
    for(i in 1:length(X.srcLis)){
      system(paste0("rm ", get_path ,"/TLinit_src_", i, "*"))
    }
  }
  
  return(list(w0=w0, delta0=delta0, beta0 = beta0, lam.const=lam.const))
}

#' Trans.global
#'
#' @param X.tarLis X.tarLis = paste0(XXX, ".vcf")
#' @param y.tarLis y.tarLis = paste0(XXX, ".phe")
#'                 .phe file can be generated from .txt with FID IID and phenotype columns by function `pheGen`
#' @param X.srcLis X.srcLis = paste0(XXX, ".vcf")
#' @param y.srcLis X.srcLis = paste0(XXX, ".phe")
#'                 .phe file can be generated from .txt with FID IID and phenotype columns by function `pheGen`
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file. 
#' @param split.col the column name in the phenotype file that specifies the membership of individuals to
#'                  the training or the validation set. The individuals marked as "train" and "val" will
#'                  be treated as the training and validation set, respectively. When specified, the
#'                  model performance is evaluated on both the training and the validation sets.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param get_path  dir for saving the intermediate files and result files
#' @param configs  a list of other config parameters.
#' @param nlambda the number of lambda values - default is 50
#' @param lambda.min.ratio smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value,
#'                         i.e. the smallest value for which all coefficients are zero. The default
#'                         depends on the sample size nobs relative to the number of actual variables
#'                         nvars (after QC filtering). The default is 0.01. A very small value of lambda.min.ratio
#'                         will lead to a saturated fit in the nobs < nvars case.
#' @param p number of SNPs in genofile
#' @param family the type of the phenotype: "gaussian", "binomial".
#' @param id_type "underline" or "double", depend on the sample names in VCF files. The default is "double".
#'                if the sample name is 'FID_IID', choose id_type="underline"
#'                if the sample name is 'IID' alone , choose id_type="double"
#' @param agg whether to remove aggregate results after the end of the functions. default is T
#' @param ValProp the proportion of validating set for aggregation when spliting the primary sample files for later aggregation. default is 1/10
#' @param alpha the elastic-net mixing parameter, where the penalty is defined as
#'              alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2. alpha = 1 corresponds to the lasso penalty,
#'              while alpha = 0 corresponds to the ridge penalty.
#' @param nsplit split number of X and beta that we can product part by part to reduce memory use. default is 4 
#' @param shortcut whether to avoid using tuned lambda to caltulate the parameters on all samples again. If the
#'                 sample size is large enough, let shortcut = T can save computational costs.
#' @export
Trans.global<-function(X.tarLis, y.tarLis, X.srcLis, y.srcLis,
                       phenotype = phenotype, split.col = NULL,
                       covariates = covariates, get_path = get_path,
                       configs = configs, nlambda = 30, lambda.min.ratio = 0.01, p = p,
                       family, id_type = "double", agg=T, ValProp = 1/10, alpha = 1, nsplit = 4, shortcut = T){
  if(is.null(split.col)){
    stop("split.col cannot be NULL.\n")
  }
  ## if agg, split of training sets for estimating \beta_0 and validating set for aggregation is needed
  ## then split_for_agg = agg = T
  global.re <- TL.init(X.tarLis, y.tarLis, X.srcLis, y.srcLis,
                       phenotype = phenotype, split.col = split.col,
                       covariates = covariates, get_path = get_path,
                       configs = configs, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                       family = family, id_type = id_type, p = p, rm = F, split_for_agg = agg,
                       ValProp = ValProp, nsplit = nsplit, alpha = alpha, shortcut = shortcut)
  if(agg){
    ST.re <- ST.init(genotype.pfile = paste0(get_path , "/TLinit_tar"), 
                     phenotype.file = paste0(get_path , "/TLinit_tar.phe"), phenotype = phenotype, split.col = split.col,
                     covariates = covariates, offset = NULL, configs = configs, nlambda = nlambda, lambda.min.ratio = 0.01,
                     family = family, alpha = alpha, shortcut = shortcut)
    B<-cbind(global.re$beta0, as.data.frame(global.re$w0), ST.re$beta0)
    # cat(colnames(B))
    B = as.matrix(B)
    # eta.hat<-agg.fun(B=B, X.til, y.til)
    y.test = read.csv(paste0(get_path , "/TLinit_val.phe"), header = T)
    if(family == "binomial"){
      eta.hat<-Binary_agg.fun(B=B, Xtest.pfile = paste0(get_path , "/TLinit_val"), y.test, nsplit = nsplit, phenotype = phenotype, covariates = covariates)
    }else if(family == "gaussian"){
      eta.hat<-agg.fun(B=B, Xtest.pfile = paste0(get_path , "/TLinit_val"), y.test, nsplit = nsplit, phenotype = phenotype, covariates = covariates)
    }
    beta.hat<-as.matrix(as.matrix(B)[, which(!is.nan(eta.hat))])%*%eta.hat[which(!is.nan(eta.hat))]
  }else{
    beta.hat=global.re$beta0
    eta.hat=1
  }
  # if(rm == F){
  system(paste0("rm ", get_path, "/TLinit_tarXb.phe"))
  system(paste0("rm ", get_path, "/global*"))
  ## not sparse
  system(paste0("rm ", get_path ,"/TLinit_tar*"))
  system(paste0("rm ", get_path ,"/TLinit_val*"))
  system(paste0("rm ", get_path ,"/keep_tar.txt"))
  system(paste0("rm ", get_path ,"/keep_val.txt"))
  for(i in 1:length(X.srcLis)){
    system(paste0("rm ", get_path ,"/TLinit_src_", i, "*"))
  }
  # }
  return(list(beta.hat=beta.hat, eta.hat=eta.hat))
}
