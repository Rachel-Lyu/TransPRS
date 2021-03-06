% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TransPRS.R
\name{TL.init}
\alias{TL.init}
\title{TL.init}
\usage{
TL.init(
  X.tarLis,
  y.tarLis,
  X.srcLis,
  y.srcLis,
  phenotype = phenotype,
  split.col = NULL,
  covariates = covariates,
  get_path = get_path,
  configs = configs,
  nlambda = 30,
  lambda.min.ratio = 0.01,
  p = NULL,
  family,
  id_type = "double",
  rm = F,
  split_for_agg = F,
  ValProp = 1/10,
  nsplit = 4,
  alpha = 1,
  shortcut = T
)
}
\arguments{
\item{X.tarLis}{X.tarLis = paste0(XXX, ".vcf")}

\item{y.tarLis}{y.tarLis = paste0(XXX, ".phe")}

\item{X.srcLis}{X.srcLis[\link{i}] = paste0(XXX, ".vcf")}

\item{y.srcLis}{X.srcLis[\link{i}] = paste0(XXX, ".vcf")}

\item{phenotype}{the name of the phenotype. Must be the same as the corresponding column name in
the phenotype file.}

\item{split.col}{the column name in the phenotype file that specifies the membership of individuals to
the training or the validation set. The individuals marked as "train" and "val" will
be treated as the training and validation set, respectively. When specified, the
model performance is evaluated on both the training and the validation sets.}

\item{covariates}{a character vector containing the names of the covariates included in the lasso
fitting, whose coefficients will not be penalized. The names must exist in the
column names of the phenotype file.}

\item{get_path}{dir for saving the intermediate files and result files}

\item{configs}{a list of other config parameters.}

\item{nlambda}{the number of lambda values - default is 50}

\item{lambda.min.ratio}{smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value,
i.e. the smallest value for which all coefficients are zero. The default
depends on the sample size nobs relative to the number of actual variables
nvars (after QC filtering). The default is 0.01. A very small value of lambda.min.ratio
will lead to a saturated fit in the nobs < nvars case.}

\item{p}{number of SNPs in genofile}

\item{family}{the type of the phenotype: "gaussian", "binomial".}

\item{id_type}{"underline" or "double", depend on the sample names in VCF files
if the sample name is 'FID_IID', choose id_type="underline"
if the sample name is 'IID' alone , choose id_type="double"}

\item{rm}{whether to remove intermediate files after the end of the functions. default is F}

\item{split_for_agg}{whether to split the primary sample files for later aggregation. default is F}

\item{ValProp}{the proportion of validating set for aggregation when spliting the primary sample files for later aggregation. default is 1/10}

\item{nsplit}{split number of X and beta that we can product part by part to reduce memory use. default is 4}

\item{alpha}{the elastic-net mixing parameter, where the penalty is defined as
alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2. alpha = 1 corresponds to the lasso penalty,
while alpha = 0 corresponds to the ridge penalty.}

\item{shortcut}{whether to avoid using tuned lambda to caltulate the parameters on all samples again. If the
sample size is large enough, let shortcut = T can save computational costs.}
}
\description{
TL.init
}
