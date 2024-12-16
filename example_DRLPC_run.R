# List of required libraries for DRLPC
required_packages <- c(
  "RNetCDF",   # For reading/writing NetCDF files
  "glmnet",    # For elastic net regression
  "car",       # For calculating Variance Inflation Factor (VIF)
  "igraph",    # For graph-based clustering of SNPs
  "FactoMineR",# For Principal Component Analysis (PCA)
  "mgcv"       # For Generalized Additive Models (GAMs)
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load libraries
lapply(required_packages, library, character.only = TRUE)

# Check if input files exist
if (!file.exists("DRLPC/data/SNPinfo.rds") || 
    !file.exists("DRLPC/data/genotype_data.rds") || 
    !file.exists("DRLPC/data/geneSNPinfo.rds")) {
  stop("One or more input files are missing. Please check the 'data' folder.")
}

# Input files
SNPinfo <- readRDS("DRLPC/data/SNPinfo.rds")
genotype_data <- readRDS("DRLPC/data/genotype_data.rds")
geneSNPinfo <- readRDS("DRLPC/data/geneSNPinfo.rds")

# Source the DRLPC functions
source("DRLPC/source/DRLPC_all_functions.R")

# Thresholds and parameters
vifcut <- 20
CLQcut <- 0.5
pccut <- 0.8
mafcut <- 5

# Filter genes with more than one SNP
genelist=unique(geneSNPinfo$gene)
N=length(genelist)

newgenelist=NULL
for(k in 1:N){
  SNPlist=as.character(geneSNPinfo$rsID[which(geneSNPinfo$gene==genelist[k])])
  if(length(SNPlist)<=1) next
  if(length(SNPlist)>1) newgenelist=c(newgenelist,as.character(genelist[k]))
}


N=length(newgenelist)
allSNPlist=colnames(genotype_data)

verbose <- TRUE  # Set to FALSE to suppress progress messages
result_list <- vector("list", N)  # Pre-allocate a list to store results

for (k in 1:N) {
  SNPlist <- as.character(geneSNPinfo$rsID[which(geneSNPinfo$gene == newgenelist[k] & geneSNPinfo$rsID %in% allSNPlist)])
  genename <- newgenelist[k]
  size <- length(SNPlist)
  if (size == 1) next
  
  genodata <- genotype_data[, SNPlist]
  n <- dim(genodata)[1]
  SNPs <- paste("SNP", 1:dim(genodata)[2], sep = "")
  genodata <- data.frame(genodata)
  dimnames(genodata)[[2]] <- SNPs
  PHENOTYPE <- rnorm(n, 0, 10)
  onedata <- cbind(PHENOTYPE, genodata)
  
  oneSNPinfo <- SNPinfo[which(SNPinfo$rsID %in% SNPlist), ]
  DRout <- DRLPC(onedata, oneSNPinfo, vifcut, CLQcut, pccut, Yvar = "PHENOTYPE", itermax = 1000)
  vdata <- DRout$vdata
  numaliasrmv <- length(DRout$aliasremoved)
  numremgrps <- length(DRout$remainedgrps)
  maxgrsize <- ifelse(length(DRout$LPCinfo) > 0, length(unlist(strsplit(DRout$LPCinfo[1], "-"))), 0)
  numPC <- length(DRout$RPCind)
  
  result_list[[k]] <- data.frame(
    datanum = k,
    genename = genename,
    size = size,
    numaliasrmv = numaliasrmv,
    numremgrps = numremgrps,
    maxgrsize = maxgrsize,
    numPC = numPC,
    finalvar = dim(vdata)[2] - 1,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat(paste0("Processing Gene ", k, ": ", genename, 
               " | Initial SNPs = ", size, 
               " | Final Dimension after DRLPC = ", dim(vdata)[2] - 1, "\n"))
  }
}

# Combine results at the end
result <- do.call(rbind, result_list)
print(result)
