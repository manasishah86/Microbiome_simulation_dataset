#Microbiome simulation function from Wu et al manuscript
#Based on any real microbiome dataset, estimated OTU proportions and overdispersion is obtained as a maximum likelihood estimate
#OTU proportions generated randomly from a Dirchlet distribution
#Total count of OTUs is based on a negative binomial distribution with mean Eg:1000 and size 25


simulateData <- function(nSam=100, s=12,ncluster = 20,mu = 1000, size = 25) {
  
  data(throat.tree,envir = environment())
  data(dd,envir = environment())
  
  tree <- get("throat.tree", envir  = environment())
  dd1 = get("dd",envir  = environment())
  tree.dist <- cophenetic(tree)
  obj <- pam(as.dist(tree.dist), ncluster,diss =TRUE)
  clustering <- obj$clustering
  otu.ids <- tree$tip.label
  
  p.est = dd1$pi
  names(p.est) <- names(dd1$pi)
  theta <- dd1$theta
  gplus <- (1 - theta) / theta
  p.est <- p.est[otu.ids]
  g.est <- p.est * gplus
  p.clus <- sort(tapply(p.est, clustering, sum), decreasing=T)
  scale2 = function(x)as.numeric(scale(x))
  
  
  comm <- matrix(0, nSam, length(g.est))
  rownames(comm) <- 1:nrow(comm)
  colnames(comm) <- names(g.est)
  # comm.p hold the underlying proportions
  comm.p <- comm
  nSeq <- rnbinom(nSam, mu = mu, size = size)
  for (i in 1:nSam) {
    comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
  }
  
  otu.ids <- names(which(clustering == s))
  
  # No additional covariates in this case.
  OTU = comm[, otu.ids]
  return(list(informative.OTU = OTU,whole.OTU = OTU))
}

#SIMULATE metadata 
# Prepare covariates 
set.seed(1)
Male = rep(0, nrow(throat.meta))
Male[throat.meta$Sex == "Male"] <- 1
Smoker = rep(0, nrow(throat.meta)) 
Smoker[throat.meta$SmokingStatus == "Smoker"] <- 1 

# Simulate outcomes 
# Here, outcome is associated with covariates but unassociated with microbiota
# 33% censoring 
SurvTime <- rexp(60, (1 + Smoker + Male))
CensTime <- rexp(60, 0.75)
Delta <- as.numeric( SurvTime <= CensTime )
ObsTime <- pmin(SurvTime, CensTime)
```

######################################################################################################################################
#Microbiome simulation from waste not want not

################################################################################
#
# Differential Abundance
#
################################################################################
library("phyloseq")
data("GlobalPatterns")
SkinReal = subset_samples(GlobalPatterns, SampleType=="Skin")
SkinTop = names(sort(taxa_sums(SkinReal), decreasing=TRUE)[1:6])
SkinMat = round(as(otu_table(prune_taxa(SkinTop, SkinReal)), "matrix")/1000, 0)
# Write to csv
write.csv(SkinMat, file="sim_diff_abund_matrix_ex_real_matrix.csv", col.names=FALSE, row.names=FALSE)
# Create the rowsum version (multinomial)
SkinMatRS = matrix(rowSums(SkinMat), nrow=nrow(SkinMat), dimnames=list(rownames(SkinMat)))
write.csv(SkinMatRS, file="sim_diff_abund_matrix_ex_real_matrix_RS.csv", col.names=FALSE, row.names=FALSE)
# Simulate by sampling from multinomial
# Example of "simulated" count matrix with 4 columns/samples each class
J = 4
NL = sample(round(sample_sums(GlobalPatterns)/10000, 0), size=J*2, replace=TRUE)
# Simulate Test and NULL
nullmat = sapply(NL, function(NL, SkinMatRS){
  table(sample(rownames(SkinMatRS), NL, replace=TRUE, prob=SkinMatRS))[rownames(SkinMatRS)]
}, SkinMatRS)
nullmat[is.na(nullmat)] <- 0L
write.csv(nullmat, file="sim_diff_abund_matrix_before_effect.csv", col.names=FALSE, row.names=FALSE)
# Apply "effect" to random (arbitrary in this case) rows
testmat = nullmat
EffectSize = 10
effectrows = c(1, 4)
effectcols = 1:J
testmat[effectrows, effectcols] <- EffectSize * testmat[effectrows, effectcols]
# write to csv
write.csv(testmat,  file="sim_diff_abund_matrix_after_effect.csv", col.names=FALSE, row.names=FALSE)
######################################################################################################################################

`r opts_chunk$set(cache=FALSE, fig.width=10, message=FALSE, warning=FALSE)`
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>
  
  
  Simulation by Random Subsampling, Comparison of Normalization
========================================================
  
  ---
  ## Load Required Packages
  
  Clear workspace prior to run.

```{r clear-workspace}
rm(list=ls())
```

```{r check-install-load-packages, warning=FALSE, message=FALSE}
# The required package list:
reqpkg = c("edgeR", "PoiClaClu", "phyloseq", "DESeq", 
           "foreach", "doParallel",
           "plyr", "reshape2", "ggplot2", "grid", "scales", "cluster")
# Check against installed packages:
inpkg = installed.packages()[, "Package"]
neededpkg = reqpkg[!reqpkg %in% inpkg]
if(length(neededpkg) > 0){
  stop(paste("Need to install the following package:", neededpkg))
}
# Load all required packages and show version
for(i in reqpkg){
  print(i)
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
  packageVersion(i)
}
# Show session info
sessionInfo()
```

Load some theme parameters for ggplot2.

```{r ggplot2-params}
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}
```


---
  
  ## Set parameters for this simulation
  
  ```{r sim-params}
set.seed(20140206) 
savedatafile = "simulation-cluster-accuracy"
# load(savedatafile)
sampletypes = c("Feces", "Ocean")

# Number of OTUs in simulation.
# This is the maximum. 
# For less-diverse template, or simulated samples with fewer reads, 
# the actual number of OTUs will be much less. 
nOTUs = 2000L 
#nOTUs = 603L

# Minimum number of reads to consider an OTU "observed" in a sample
minobs= 1L

# Samples per simulation
J = 40L 
#J = 7L 

# Effect Size. The artificial mix fraction.
mixfacs = c(1, 1.15, 1.25, 1.5, 1.75, 2, 2.5, 3.5)
#mixfacs = c(1, 1.25, 2)

ns = c(1000, 2000, 5000, 1E4)
#ns = c(1000, 50000)

# Vector of the replicate numbers to repeat for
# each comb of simulation parameters (n, etc)
reps=1:5
#reps=1:2

# The delimiter for command parameters
comdelim = "_"

# Number of cores. Only relevent for certain types of back-end registration.
Ncores = 8
#Ncores = 7

# rarefying power params
rarefy_fracs = c(0, seq(0.05, 0.25, 0.05), 0.4)
#rarefy_fracs = c(0, 0.4)

rarereps = 1:2
#rarereps = 1

# The combinations of simulation parameters, 
# used later in the script.
simparams = apply(expand.grid(ns, reps, mixfacs), 1, paste0, collapse=comdelim)
# Define date-stamp for file names
datestamp = gsub(":", "_", gsub("[[:space:]]+", "-", date()), fixed=TRUE)
print(datestamp)
```


## Parameter definitions

The major factors contributing to the computation cost in this simulation example are the number of OTUs retained from the template `GlobalPatterns` dataset, which ultimately is used to dictate the length of the multinomials and their corresponding proportion vectors specified by $\pi$ (R variable `pi` or `pis`); and the number of samples being simulated. 

- `nOTUs` -- the number of most-prevalent OTUs to keep in simulation template, `r nOTUs`.
- `minobs` -- The minimum abundance value, `r minobs`, for which an OTU is "counted" as having been observed in a given sample. This is used for ranking OTUs according to the number of samples in which they appeared.
- `J` -- The number of samples per simulated table. For multitable analysis this needs to be consistent across tables, so included here in the beginning as a global parameter. Value for this simulation is `r J`.
- `mixfacs` -- The fold-weighting of the template multinomial compared with its mixed-in counterpart. Same in both directions (symetric mixing). 
- `ns` -- A vector of the expected values for the sampling depth that are nevertheless subject to the random sampling from the original template totals. The read numbers should not exceed the total number of reads in the template, so this is checked in the "simulation" (subsampling) module, and a ceiling used. The values for this vector of sampling depth means in this particular simulation is: `r ns`.


---
  
  # Define template
  
  Trim to just those samples that are intended as template in this binary effect.

```{r trim-2-sampletypes}
data("GlobalPatterns")
sampsums = sample_sums(GlobalPatterns)
keepsamples = sample_data(GlobalPatterns)$SampleType %in% sampletypes
template = prune_samples(keepsamples, GlobalPatterns)
```

Trim OTUs that do not appear in very many samples in the template. Sort by prevalance (number of samples appeared) and then by abundance value.
```{r trim-by-prevalence}
samobs = apply(otu_table(template), 1, function(x, m) sum(x > m), m=minobs)
otudf = data.frame(prev=samobs, sums=taxa_sums(template))
otudf = otudf[order(-otudf$prev, -otudf$sums), ]
```

Trim all but the first nOTUs
```{r trim-template-nOTUs}
template = prune_taxa(rownames(otudf)[1:nOTUs], template)
template
```


---
  
  # Simulate microbiome census
  In a previous version I used the `simPop` function from [the dirmult package](http://cran.r-project.org/web/packages/dirmult/index.html).

Instead, here I am randomly subsampling from the original templates at different numbers of total samples.

## Define templates

```{r simPop-params}
template1 = subset_samples(template, SampleType==sampletypes[1])
template2 = subset_samples(template, SampleType==sampletypes[2])
```

## Forced mixing.
In the original version 4 result, it appeared that there was very little overlap between the two template types. This could easily be the case, even though they were both body habitats. One biologically intuitive approach would be to choose a different pair of template sample types that are more similar. However, this doesn't afford much additional intuition about the extent that they are similar (without a bunch of extra exploratory analysis), and doesn't give us very much control over the extent that they are similar... we'd be stuck with whatever similarities are "available" among the sample types. Note that it might also help to increase the number of OTUs being included in template and simulation. At the moment it is `r nOTUs`.

Instead, we will intentionally mix two templates that otherwise have very few OTUs in common, now labeled `template1` and `template2`, by adding a known but small proportion of each template to the other.

```{r force-template-mixing}
# Define function for mixing
mix_microbes = function(template1, template2, unmixfac){
require("phyloseq")
# Check that the number of taxa are equal
if( !identical(ntaxa(template1), ntaxa(template2)) ){ stop("number of OTUs b/w template1 and template2 must be equal") }
# Expects templates to be a 1-sample dataset.
if( !identical(nsamples(template1), 1L) ){stop("`template1` should have only 1 sample")}
if( !identical(nsamples(template2), 1L) ){stop("`template2` should have only 1 sample")}
# Enforce taxa_are_rows to FALSE
if(taxa_are_rows(template1)){template1 <- t(template1)}
if(taxa_are_rows(template2)){template2 <- t(template2)}
# Define a vector version of each template for subsampling
x1 = as(otu_table(template1), "numeric")
x2 = as(otu_table(template2), "numeric")
# Create 'dirty' multinomial.
# Defined artificial mixing.
# Create mixed multinomial by adding counts from the other in precise proportion,
# a total of Library Size / Effect Size
addToTemplate1 = round(( sum(x1) * x2 / (sum(x2) * unmixfac) ), 0)
# Add them together to create "dirty" multinomial
# This adds the fractional subsampled abundances to the template
mat1 = matrix((addToTemplate1 + x1), nrow=1)
rownames(mat1) <- sample_names(template1)
colnames(mat1) <- taxa_names(template1)
otu_table(template1) <- otu_table(mat1, taxa_are_rows=FALSE)
return(template1)
}
# merge the template components together into one sample.
template1 = merge_samples(template1, "SampleType")
template2 = merge_samples(template2, "SampleType")
```


---

## Simulation Function

- Input
- `postfix`, `template`, `J`, `n` (number or reads)
- Output
- phyloseq object, incorporating simulation results and inputs

```{r define-microbesim}
microbesim = function(postfix="sim", template, templatex, unmixfac, J, n=1E4){
# Generate `J` simulated microbiomes with `n` total reads each
# (all the same, or n has length equal to the value of `J`),
# with subsamples drawn from `template`.
# `postfix` is a dummy idenitifer added to help distinguish
# simulated samples in downstream code.
require("phyloseq")
# Perform the mixing here, so each replicate is a separate random mix
template = mix_microbes(template, templatex, unmixfac)
# call the proporitions vector `pi`, similar to nomenclature from DMN
pi = taxa_sums(template)
# n must be a scalar (recycled as the number of reads for every simulation)
# or it can be vector of length equal to J, the number of samples being simulated.
if(length(J) != 1){ stop("Length of J should be 1.") }
if(length(n) != 1 & length(n) != J){
stop("n should be length 1, or length J.")
}
# Actually create the simulated abundance table
simat = mapply(function(i, x, sample.size){
if(FALSE){print(i)} # i is a dummy iterator
phyloseq:::rarefaction_subsample(x, sample.size)
}, i=1:J, sample.size=n, MoreArgs=list(x=pi), SIMPLIFY=TRUE)
simat = t(simat)
# Add the OTU names to the OTU (column) indices
colnames(simat) <- names(pi)
# Add new simulated sample_names to the row (sample) indices
rownames(simat) <- paste(i, "::", 1:nrow(simat), postfix, sep="")
# Put simulated abundances together with metadata as a phyloseq object
OTU = otu_table(simat, taxa_are_rows=FALSE)
# Define data.frame that will become sample_data
SDF = data.frame(sample=sample_names(OTU), TableNumber=i, type="simulated", stringsAsFactors=FALSE)
SDF$postfix <- postfix
rownames(SDF) <- sample_names(OTU)
SD  = sample_data(SDF)
# Return a phyloseq object
return(phyloseq(OTU, SD))
}
```

Define a function for sampling from the library sizes as a way of including the "noise" derived from this aspect of the data as well.

```{r simulate-sizes}
# rescale the sum of reads in the original raw(-ish) template data
# to the expected library size being requested here
sumsim = function(n, sumtemplate, J){
# `n` - expected size target
# `sumtemplate` - the template vector of library sizes observed in template
# `J` - The number of sample sizes to return
####
# sumtemplate = sampsums
# n = 2000
# J = 103
scaledSums = round(n*(sumtemplate/median(sumtemplate)))
return(sample(scaledSums, size=J, replace=TRUE))
}
```


---

## Register parallel backend for computing

Register parallel clusters for parallel calculations (e.g. UniFrac, etc.).

```{r register-parallel-cluster, message=FALSE}
# makeCluster
# makeForkCluster()

# The parallel single-node way
cl <- makeCluster(Ncores)

# The multi-node way for PSOCK cluster, from:
# http://www.bioconductor.org/help/bioconductor-cloud-ami/#parallel
# library(parallel)
# lines <- readLines("/usr/local/Rmpi/hostfile.plain")
# hosts <- character()
# for (line in lines){
#     x <- (strsplit(line[[1]], " "))
#     hosts <- c(hosts, rep.int(x[[1]][1], as.integer(x[[1]][2])))
# }
# cl <- makePSOCKcluster(hosts, master=system("hostname -i", intern=TRUE))
# # system.time(clusterCall(cl, Sys.sleep, 1)) # A test

# Register the parallel cluster
registerDoParallel(cl)

# The Rmpi/doMPI way (doMPI package should have been loaded already at the beginning)
# cl <- startMPIcluster(verbose=TRUE)
# registerDoMPI(cl)
```


---

## Simulate microbiome samples from different 'dirty' templates

Repeat the simulation many times, with different values for the number of reads per sample. 

```{r micsimtest}
# Parallelized simulations
simlist <- foreach(i=simparams, .packages=c("phyloseq")) %dopar% {
# i = simparams[4]
# Initialize
n = sim = sim1 = sim2 = n1 = n2 = NULL
#cat(i, "\n")
n = as.numeric(strsplit(i, comdelim)[[1]][1])
mixfac = as.numeric(strsplit(i, comdelim)[[1]][3])
# Rarely a simulation has a weird value and fails.
# Catch these with `try`, and repeat the simulation call
# if error (it will be a new seed)
tryAgain = TRUE; infiniteloopcounter = 1
while(tryAgain & infiniteloopcounter < 5 ){
n1   = sumsim(n, sampsums, J)
n2   = sumsim(n, sampsums, J)
sim1 = microbesim(sampletypes[1], template1, template2, mixfac, J, n1)
sim2 = microbesim(sampletypes[2], template2, template1, mixfac, J, n2)
if( is.null(sim1) | is.null(sim2) | 
is.null(n1) | is.null(n2) | 
inherits(sim1, "try-error") | inherits(sim2, "try-error")){
tryAgain = TRUE
infiniteloopcounter = infiniteloopcounter + 1
} else {
tryAgain = FALSE
}
}
if( infiniteloopcounter >= 5 ){
stop("Consistent error found during simulation. Need to investigate cause.")
}
# Merge the two simulated datasets together into one phyloseq object
# and add back tree.
sim = merge_phyloseq(sim1, sim2)
sim = merge_phyloseq(sim, tax_table(GlobalPatterns), phy_tree(GlobalPatterns))	
return(sim)
}
names(simlist) <- simparams
```


---

# Evaluate clustering after different normalizations, metrics

This is the point of this analysis. Clustering method must be the same (and a good method) for all versions of the simulated dataset.

The following sections define functions that normalize the samples to equal numbers of reads per sample in different ways. 


## Naïve proportions

Need to "preprocess" the simulated microbiome data. Adapted code on standardizing/normalizing/centering(not supported) data. **Rm zeros** Remove the zeros. That is, remove any OTUs with a zero-sum across all samples. Hopefully not many. Also do the same for simulated samples. Both could cause problems later. **Normalize**. Normalize sampling effort of each sample. This won't correct for unobserved OTUs, which will still be zero, but an OTU trimming step helps correct that. **Scale**. Scale by dividing each variable by its standard deviation. **Center**. Center by subtracting the median of each sample from each element of that sample. Skip this for now. The phyloseq package currently forbids negative abundance values. Maybe this should be changed to a warning rather than an error in the package (or get rid of that condition)...

We can take as a given that samples or OTUs that are not reprsented in datasets (did not get any simulated reads), would be cut in any preprocessing/normalization method.

```{r remove-empty}
# Rm any samples or OTUs that don't have any reads.
# Also remove singletons (only 1 count as sample or OTU)
remove_singletons_empties = function(physeq){
  require("phyloseq")
  physeq = prune_taxa(taxa_sums(physeq) > 2.5, physeq)
  physeq = prune_samples(sample_sums(physeq) > 2.5, physeq)
  # Remove OTUs not appearing in at least 3 samples
  if( taxa_are_rows(physeq) ){
    y = otu_table(physeq)
    wh1 = apply(aaply(y, 1, function(x){x >= 1}), MARGIN=1, sum) >= 3
  } else {
    y = as(t(otu_table(physeq)), "matrix")
    wh1 = apply(aaply(y, 1, function(x){x >= 1}), MARGIN=1, sum) >= 3
  }
  physeq = prune_taxa(wh1, physeq)
  return(physeq)
}
# Remove singletons and empties
simnames = names(simlist)
simlist = foreach(physeq=simlist, .packages=c("phyloseq", "plyr")) %dopar% {remove_singletons_empties(physeq)}
names(simlist) <- simnames
```

Define the naïve (simple proportion) normalization function.

```{r rm-zeros-sim}
proportion = function(physeq){
  # Normalize total sequences represented
  normf = function(x, tot=max(sample_sums(physeq))){ tot*x/sum(x) }
  physeq = transform_sample_counts(physeq, normf)
  # Scale by dividing each variable by its standard deviation.
  #physeq = transform_sample_counts(physeq, function(x) x/sd(x))
  # Center by subtracting the median
  #physeq = transform_sample_counts(physeq, function(x) (x-median(x)))
  return(physeq)
}
```

### Random subsampling

Remove some of the lowest-abundance samples (e.g. poorest-represented quartile of samples) prior to random subsampling.

```{r random-subsample-function}
randomsubsample = function(physeq, smalltrim=0.15, replace=TRUE){
  require("phyloseq")
  # Set the minimum value as the smallest library quantile, n`smalltrim` 
  samplemin = sort(sample_sums(physeq))[-(1:floor(smalltrim*nsamples(physeq)))][1]
  physeqr = rarefy_even_depth(physeq, samplemin, rngseed=FALSE,
                              replace=replace, trimOTUs=TRUE)
  return(physeqr)
}
```


### edgeR normalization 

Normalizations provided by the edgeR pacakge. Most notably, "TMM" normalization.

```{r edgeR-norm, fig.keep='none', warning=FALSE}
edgeRnorm = function(physeq, ...){
  require("edgeR")
  require("phyloseq")
  #physeq = simlist[["55000_1e-04"]]
  #z0 = simlisttmm[["55000_1e-04"]]
  #physeq = simlist[["1000_0.2"]]
  #z0 = simlisttmm[["1000_0.2"]]
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # See if adding a single observation, 1, 
  # everywhere (so not zeros) prevents errors
  # without needing to borrow and modify 
  # calcNormFactors (and its dependent functions)
  # It did. This fixed all problems. 
  # Can the 1 be reduced to something smaller and still work?
  x = x + 1
  # Now turn into a DGEList
  y = edgeR::DGEList(counts=x, remove.zeros=TRUE)
  # Perform edgeR-encoded normalization, using the specified method (...)
  z = edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  # Don't need the following additional steps, which are also
  # built-in to some of the downstream distance methods.
  #z1 = estimateCommonDisp(z)
  #z2 = estimateTagwiseDisp(z1)
  return(z)
}
```


### DESeq normalization

From DESeq:
  "The inference in DESeq relies on an estimation of the typical relationship between the data’s variance and their mean, or, equivalently, between the data’s dispersion and their mean. The dispersion can be understood as the square of the coefficient of biological variation. So, if a gene’s expression typically differs from replicate to replicate sample by 20%, this gene’s dispersion is 0.22 = .04. Note that the variance seen between counts is the sum of two components: the sample-to-sample variation just mentioned, and the uncertainty in measuring a concentration by counting reads. The latter, known as shot noise or Poisson noise, is the dominating noise source for lowly expressed genes. The former dominates for highly expressed genes. The sum of both, shot noise and dispersion, is considered in the differential expression inference. Hence, the variance v of count values is modelled as:

v = sμ + αs2μ2

$$
\upsilon = s\mu + \alpha{s}^2{\mu}^2
$$

where $\mu is the expected normalized count value (estimated by the average normalized count value), $s$ is the size factor for the sample under consideration, and $\alpha$ is the dispersion value for the gene under consideration. 
"

To estimate the dispersions, use the `estimateDispersions` function. Note that there are `method`, `sharingMode`, and `fitType` arguments to `estimateDispersions`, the choice of which are presumed to affect power. For the sample clustering simulation here, we have no "replicate" samples from a particular treatment or technical step. 

```{r deseq-variance-stabilization-transform}
deseq_varstab = function(physeq, sampleConditions=rep("A", nsamples(physeq)), ...){
  require("DESeq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # The same tweak as for edgeR to avoid NaN problems
  # that cause the workflow to stall/crash.
  x = x + 1
  # Create annotated data.frame with the taxonomy table, in case it is useful later
  taxADF = as(data.frame(as(tax_table(physeq), "matrix"), stringsAsFactors=FALSE), "AnnotatedDataFrame")
  cds = newCountDataSet(x, sampleConditions, featureData=taxADF)
  # First estimate library size factors
  cds = estimateSizeFactors(cds)
  # Variance estimation, passing along additional options
  cds = estimateDispersions(cds, ...)
  # Determine which column(s) have the dispersion estimates
  dispcol = grep("disp\\_", colnames(fData(cds)))
  # Enforce that there are no infinite values in the dispersion estimates
  if( any(!is.finite(fData(cds)[, dispcol])) ){
    fData(cds)[which(!is.finite(fData(cds)[, dispcol])), dispcol] <- 0.0
  }
  vsmat = exprs(varianceStabilizingTransformation(cds))
  otu_table(physeq) <- otu_table(vsmat, taxa_are_rows=TRUE)
  return(physeq)
}
```

From DESeq vignette:
  
  "As we estimated the dispersions from a small number of replicates, the estimates scatter with quite some sampling variance around their true values. An initial assumption that one could make is that the regression line shown in Figure 1 models the true underlying dispersions, and that the variation of the point estimates around simply reflects sampling variance. This is the assumption that we put forward in the first paper on DESeq. However, subsequent experience with larger data sets indicates that not all of the variability of the points around the regression line seen in Figure 1 is sampling variance: some of it reflects differences of the true, underlying variance between different genes. Hence, the default behaviour of DESeq now uses a more prudent or conservative approach: if a per-gene estimates lies below the regression line, we assume that this might indeed be sampling variance, and shift the estimate upwards to the value predicted by the regression line. If, however, the per-gene estimate lies above the line, we do not shift it downwards to the line, but rather keep it as is.

The option `sharingMode` of the function estimateDispersions can be used to control this behaviour. The value `sharingMode="maximum"` corresponds to the default. If you have many replicates (where many starts at around 7), the choice `sharingMode="gene-est-only"` will typically be more adequate. If you would like to use the behaviour described in [the original DESeq article], this is achieved by specifying `sharingMode="fit-only"`.

Another difference of the current DESeq version to the original method described in the paper is the way how the mean-dispersion relation is fitted. By default, estimateDispersions now performs a parametric fit: Using a gamma-family GLM, two coefficients $\alpha_{0}$, $\alpha_{1}$ are found to parametrize the fit as $\alpha = \alpha_{0} + \alpha_{1}/\mu$. (The values of the two coefficients can be found in the `fitInfo` object, as attribute coefficients to `dispFunc`.) For some data sets, the parametric fit may give bad results [or just fail], in which case one should try a local fit (the method described in the paper), which is available via the option `fitType="local"` to `estimateDispersions`.

In any case, the dispersion values which will be used by the subsequent testing are stored in the feature data slot of cds [use `fData`]."

To adjust the dispersion estimates, fix infinite values, use the `fData` accessor prior to calling the testing function.


### Perform normalizations 

```{r simpletrim}
simpletrim = function(physeq, J){
  if( taxa_are_rows(physeq) ){
    physeq = t(physeq)
  }
  # `prevalence` is the fraction of total samples in which
  # an OTU is observed at least `minobs` times.
  prevalence = apply(as(otu_table(physeq), "matrix"), 2, 
                     function(x, minobs){sum(x > minobs)},
                     minobs)/(2*J)
  # Will only keep OTUs that appear in more than X% of samples
  keepOTUs = prevalence > 0.05
  # Keep only those OTUs with total reads greater
  # than 0.5x the number of samples. 
  keepOTUs = keepOTUs & taxa_sums(physeq) > (0.5*J)
  return(prune_taxa(keepOTUs, physeq))
}
```

Perform initializations, and trimming

```{r perform-initializations-trim}
# Initialize the simulation lists simlist0 and simlist
simlist0 = simlistuntrim =  simlist
# checks
#apply(as(otu_table(simlistuntrim[[1]]), "matrix"), 2, 
#      function(x, tot){sum(x>tot)}, minobs)/(2*J)
#sapply(lapply(simlistuntrim, taxa_sums), function(x, J){sum(x < (0.5*J))}, J)/nOTUs
# Trim simlist0
simlist0 = lapply(simlistuntrim, simpletrim, J)
```

Now perform the normalizations.

```{r perform-normalizations}
# Replace simlist with empty list. Will add different norm types to it.
simlist = vector("list", 0)
simlist$none <- simlist0

# Define the "naive proportion" list
simlist$naiveProp <- foreach(physeq=simlist0, .packages="phyloseq") %dopar% {proportion(physeq)}
names(simlist$naiveProp) <- names(simlist0)

# rarefying.
simlist$rarefy <- foreach(physeq=simlist0, .packages="phyloseq") %dopar% {
  randomsubsample(physeq, smalltrim=0.15, replace=FALSE)
}
names(simlist$rarefy) <- names(simlist0)

# edgeR norms
simlist$edgeRTMM <- foreach(physeq=simlist0, .packages=c("phyloseq", "edgeR")) %dopar% {edgeRnorm(physeq, method="TMM")}
names(simlist$edgeRTMM) <- names(simlist0)
simlist$edgeRRLE <- foreach(physeq=simlist0, .packages=c("phyloseq", "edgeR")) %dopar% {edgeRnorm(physeq, method="RLE")}
names(simlist$edgeRRLE) <- names(simlist0)
simlist$edgeRupqua <- foreach(physeq=simlist0, .packages=c("phyloseq", "edgeR")) %dopar% {edgeRnorm(physeq, method="upperquartile")}
names(simlist$edgeRupqua) <- names(simlist0)

# DESeq variance stabilizing transformations
simlist$DESeqVS <- foreach(physeq=simlist0, .packages=c("phyloseq", "DESeq")) %dopar% {
  deseq_varstab(physeq, method="blind", sharingMode="maximum", fitType="local")
}
names(simlist$DESeqVS) <- names(simlist0)
```


## Save list of normalized simulations

```{r save-simulations-step}
save.image(paste0("simulations-", datestamp, ".RData"))


