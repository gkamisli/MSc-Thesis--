def tf_Pairs_Finding(transcriptionFactor, sTFs_nP, sTargets_nP, sGeneSymbols_nR):
    """
    Function to search for target in genes (sGeneSymbols_nR) by comparing with known transcription factor (TF) - target pairs in 
    GWAS study and get their indexes 
    
    Argument:
    transcriptionFactor -- transcription factor (TF) to be identified
    sTFs_nP -- known "TFs" from GWAS study to be searched
    sTargets_nP -- known "targets" from GWAS study to be searched
    sGeneSymbols_nR -- all genes obtained from a given data set (i.e., HBTRC)
    
    Returns: 
    tfTargetList -- a list indexes of targets for a given TF so as they are valid pairs
    """
    import numpy as np
    
    index = np.where(transcriptionFactor == sTFs_nP)[0]
    tf_targets = sTargets_nP[index]
    targetCheck = np.array( [tf_targets[i] in sGeneSymbols_nR for i in range(tf_targets.shape[0])] )
    validPairs = index[np.where(targetCheck == True)[0]]
    tfTargetList = validPairs.tolist()
    
    return tfTargetList