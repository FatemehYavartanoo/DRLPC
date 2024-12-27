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

#NOTE: If you have already downloaded the input files, change the directory to your local path.
# Define file paths
local_data_dir <- "DRLPC/data"

# Ensure input files exist
if (!file.exists(file.path(local_data_dir, "SNPinfo.rds")) || 
    !file.exists(file.path(local_data_dir, "genotype_data.rds")) || 
    !file.exists(file.path(local_data_dir, "geneSNPinfo.rds"))) {
  stop("One or more input files are missing. Please check the 'data' folder or update the directory path.")
}

# Load input files
SNPinfo <- readRDS(file.path(local_data_dir, "SNPinfo.rds"))
genotype_data <- readRDS(file.path(local_data_dir, "genotype_data.rds"))
geneSNPinfo <- readRDS(file.path(local_data_dir, "geneSNPinfo.rds"))

# URL of the source code on GitHub
url <- "https://raw.githubusercontent.com/FatemehYavartanoo/DRLPC/main/source/DRLPC_all_functions.R"

# Temporary file to save the script
temp_file <- tempfile(fileext = ".R")

# Download the file
download.file(url, temp_file, mode = "wb")

# Source the downloaded file
source(temp_file)

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
cat("The number of genes with more than one SNP is:", N, "\n")

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
