# Dimension Reduction using Local Principal Components for Regression-based Multi-SNP Analysis in 1000 Genomes and the Canadian Longitudinal Study on Aging (CLSA)
Dimension Reduction using Local Principal Components (DRLPC) is an approach for gene-level association analysis based on Local PCA. DRLPC clusters correlated SNPs, replaces each cluster with a local principal component, and iteratively removes variables with high Variance Inflation Factor (VIF) to address multicollinearity. 

Detailed information about [DRLPC: Dimension Reduction using Local Principal Components] can be found in our [preprint](https://www.biorxiv.org/content/10.1101/2024.05.13.593724v1.abstract).

#Overviwe of DRLPC algorithm

DRLPC algorithm
	Step1: SNPs are clustered using clique-based algorithm (CLQ).
	Step2: Clusters are replaced by Local Pcs.
	Step3: Removed variables with high VIF in updated dataset.
	Step4: Add additional PCs obtained by removed SNPs. 
![image](https://github.com/user-attachments/assets/9375164a-02ef-4b2f-97a5-56d63b91daa3)

# Citation

If you use this repository, please cite:

Dimension Reduction using Local Principal Components for Regression-based Multi-SNP Analysis in 1000 Genomes and the Canadian Longitudinal Study on Aging (CLSA). Yavartanoo, F., Brossard, M., Bull, S. B., Paterson, A. D., & Yoo, Y. J. DOI:
[https://doi.org/10.1101/2024.05.13.593724v1](https://www.biorxiv.org/content/10.1101/2024.05.13.593724v1.abstract)

# License

This package is released under the [GNU General Public License (GPL) v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).


