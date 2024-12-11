# Dimension Reduction using Local Principal Components for Regression-based Multi-SNP Analysis in 1000 Genomes and the Canadian Longitudinal Study on Aging (CLSA)
Dimension Reduction using Local Principal Components (DRLPC) is an approach for gene-level association analysis based on Local PCA. DRLPC enhances the power of multi-SNP test statistics while preserving localized effect interpretability. It clusters highly correlated SNPs, replaces each cluster with a local principal component, and iteratively removes variables with high Variance Inflation Factor (VIF) to address multi-collinearity before regression analysis.

Detailed information about [DRLPC: Dimension Reduction using Local Principal Components] can be found in our [preprint](https://www.biorxiv.org/content/10.1101/2024.05.13.593724v1.abstract).

## **Overviwe of DRLPC algorithm** ##

### **DRLPC algorithm** ###
Step 1: SNPs are clustered using the clique-based algorithm (CLQ).  
Step 2: Clusters are replaced by Local PCs.  
Step 3: Removed variables with high VIF in the updated dataset.  
Step 4: Add additional PCs obtained by the removed SNPs.

<div align="center">
    <img src="DRLPC-algorithm.png" alt="DRLPC Algorithm" width="600">
</div>



### **Usage Instructions** ###
To use the DRLPC function:
1. Download the `DRLPC_main_function.R` script from this repository.
2. Follow the instructions in the [DRLPC_UserGuide.pdf](DRLPC_UserGuide.pdf) to prepare your data and run the analysis.



## **Citation** ##

If you use this repository, please cite:

Dimension Reduction using Local Principal Components for Regression-based Multi-SNP Analysis in 1000 Genomes and the Canadian Longitudinal Study on Aging (CLSA). Yavartanoo, F., Brossard, M., Bull, S. B., Paterson, A. D., & Yoo, Y. J. DOI:
[https://doi.org/10.1101/2024.05.13.593724v1](https://www.biorxiv.org/content/10.1101/2024.05.13.593724v1.abstract)

## **License** ##

This package is released under the [GNU General Public License (GPL) v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).


