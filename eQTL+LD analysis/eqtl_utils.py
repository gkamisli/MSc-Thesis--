# Utility function for Within-LD pairs analysis

def dropWithinPairsLD(df, gene, snps, cisWithinLDSNP_dCN):
    """
    Function to drop variants (SNPs) that are high LD (R^2 > 0.6) within the same gene-SNP pair (cis-eQTL table)
    
    Arguments:
    df -- cis-eQTL dataframe having significant SNPs, genes, statistics, pvalue, FDR, MAF columns 
    gene -- cis-gene to be searched if its variants are high LD or not
    snps -- SNPs of cis-gene 
    
    Return: 
    df -- cis-eQTL dataframe after SNPs with high LD are removed
    
    """
    
    import pandas as pd
    for i in snps:
        for j in snps:
            if i in cisWithinLDSNP_dCN.keys():
                if j in cisWithinLDSNP_dCN[i]:
                    df.drop(df[(df['snps'] == i) & (df['gene'] == gene)].index, inplace=True)
                    df.drop(df[(df['snps'] == j) & (df['gene'] == gene)].index, inplace=True)
    return df