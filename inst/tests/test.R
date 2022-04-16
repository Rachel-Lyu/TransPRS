# ~/.conda/envs/sy_Env/bin/Rscript ~/TransPRS/inst/tests/test.R
# source(paste0(system.file(package = 'TransPRS'), "/tests/test.R"))
library(TransPRS)
plink2.path = "~/TransLasso+snpnet/plink2"
zstdcat.path = "~/.conda/envs/sy_Env/bin/zstdcat"
zstd.path = "~/.conda/envs/sy_Env/bin/zstd"
get_path <- paste0(system.file(package = 'TransPRS'), "/extdata")
configs <- list(
  # results.dir = "PATH/TO/SAVE/DIR",  # needed when saving intermediate results
  # save = TRUE,  # save intermediate results per iteration (default FALSE)
  # nCores = 16,  # number of cores available (default 1)
  # niter = 100,  # max number of iterations (default 50)
  # prevIter = 15,  # if we want to start from some iteration saved in results.dir
  # use.glmnetPlus = TRUE,  # recommended for faster computation
  # early.stopping = FALSE,  # whether to stop based on validation performance (default TRUE)
  plink2.path = plink2.path,   # path to plink2 program
  zstdcat.path = zstdcat.path,  # path to zstdcat program
  zstd.path = zstd.path
)
# check if the provided paths are valid
for (name in names(configs)) {
  tryCatch(system(paste(configs[[name]], "-h"), ignore.stdout = T),
           condition = function(e) cat("Please add", configs[[name]], "to PATH, or modify the path in the configs list.") ) }

# phenotype.file <- paste0(get_path, "/phe_matrix.phe")
# vcf.file <- paste0(get_path, "/geno.vcf")
tarLis = paste0(get_path, "/AFR_sel_site", 1)
srcLis = list()
srcLis[[1]] = paste0(get_path, "/EAS_sel_site", 1:2)
srcLis[[2]] = paste0(get_path, "/EUR_sel_site", 1:3)
phenotype <- "QPHE"
covariates <- c("Sex", "Age")
p = 1000
## generate phenofiles that snpnet needs
for(nm in c(tarLis, unlist(srcLis))){
  pheGen(ori_phefile = paste0(nm, ".txt"), phenotype, covariates = covariates, offset = NULL, split.col = "split", out_phefile = paste0(nm, ".phe"), split_vec = c(0.9, 0.1))
  # remove annotation at first
  remove_anno(paste0(nm, ".vcf"))
  # if no pfile
  # options(warn = -1)
  pgenGen(out_pgen_name = paste0(nm), vcf.file = paste0(nm, ".vcf"), "double")
}

X.tarLis = paste0(tarLis, ".vcf")
y.tarLis = paste0(tarLis, ".phe")
X.srcLis = list()
y.srcLis = list()
for(i in 1:length(srcLis)){
  X.srcLis[[i]] = paste0(srcLis[[i]], ".vcf")
  y.srcLis[[i]] = paste0(srcLis[[i]], ".phe")
}

TG <- Trans.global(X.tarLis, y.tarLis, X.srcLis, y.srcLis,
                    phenotype = phenotype, split.col = "split",
                    covariates = covariates, get_path = get_path,
                    configs = configs, nlambda = 10, lambda.min.ratio = 0.01, p = p, 
                    family = 'gaussian', id_type = "double", 
                    agg=T, ValProp = 1/8, alpha = 1, nsplit = 5, shortcut = F)

# TG_bn <- Trans.global(X.tarLis, y.tarLis, X.srcLis, y.srcLis,
#                     phenotype = phenotype, split.col = "split",
#                     covariates = covariates, get_path = get_path,
#                     configs = configs, nlambda = 10, lambda.min.ratio = 0.01, p = p, 
#                     family = 'binomial', id_type = "double", 
#                     agg=T, ValProp = 1/8, alpha = 0.5, nsplit = 5, shortcut = T)

metric = c()
testLis = paste0(get_path, "/AFR_test_site", 1)
for(nm in c(tarLis, testLis)){
  pheGen(ori_phefile = paste0(nm, ".txt"), phenotype, covariates = covariates, offset = NULL, split.col = "split", out_phefile = paste0(nm, ".phe"), split_vec = c(0.9, 0.1))
  # remove annotation at first
  remove_anno(paste0(nm, ".vcf"))
  # if no pfile
  pgenGen(out_pgen_name = paste0(nm), vcf.file = paste0(nm, ".vcf"), "double")
  # computeMetric <- function(pred, response, metric.type = 'auc')
  # pred = Xbeta_product(nm, TG$beta.hat[-1], p, nsplit = 4) + TG$beta.hat[1]
  # response = read.csv(paste0(nm, ".phe"), header = T)
  # response = response[, phenotype]
  # metric = c(metric, (cor(pred, response))^2)
  pred = Xbeta_product(nm, TG$beta.hat[-(1:(1+length(covariates)))], p, nsplit = 4) + TG$beta.hat[1]
  response = read.csv(paste0(nm, ".phe"), header = T)
  response = response[, phenotype] - as.matrix(response[, covariates]) %*% TG$beta.hat[2:(1+length(covariates))] 
  metric = c(metric, (cor(pred, response))^2)
  # metric = c(metric, snpnet:::computeMetric(pred, response, metric.type = 'auc'))
  system(paste0("rm ", nm, ".pgen"))
  system(paste0("rm ", nm, ".psam"))
  system(paste0("rm ", nm, ".pvar*"))
}

for(nm in c(tarLis, unlist(srcLis), testLis)){
  system(paste0("rm ", nm, ".phe"))
}

cat('\nResult:', metric, '\n')
cat(TG$eta.hat, '\n')

# write.table(metric, paste0(get_path, "/TransGlobal.txt"), append = T, col.names = F)
