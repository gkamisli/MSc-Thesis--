def sort_Gene_Pairs(sGenePairs_lG):
    
    """
    Function that turns sGenePairs_lG into a sorted NumPy array based on TFs
    
    Arguments: 
    sGenePairs_lG -- all gene pairs, list, length of (num_pairs)
    
    Returns:
    sSortedGenePairs_nG -- sorted gene pairs, NumPy array, shape of (2204, 2)
    
    """
    import numpy as np
    sGenePairs_nG = np.empty(shape=[len(sGenePairs_lG), 2], dtype="<U10")
    
    for i in range(len(sGenePairs_nG)):
        sGenePairs_nG[i] = sGenePairs_lG[i][0], sGenePairs_lG[i][1] # n-array of gene pairs 

    sSortedGenePairs_nG = sGenePairs_nG[np.argsort(sGenePairs_nG[:, 0])[::1]]  # sorted n-array of gene pairs, column 1 - gene 0
    
    return sSortedGenePairs_nG

def train_Test_Data(sGenePairs_lG, uLocMax = 200):
    
    """
    Function that divides gene pairs into train and test gene pairs that no gene in the test data
    pairs included in the train gene pairs.
    
    Arguments: 
    sGenePairs_lG -- all gene pairs list of len (num_pairs)
    uLocMax -- maximum number of gene pairs to search before moving to the next TF
    
    Returns:
    sGenePairs_tr_nG -- training data gene pairs, the shape to be decided by the algorithm based on the condition
    that train data contains none of the gene pairs in the test data , NumPy array, shape of (X, 2)
    sGenePairs_tst_nG -- test data gene pairs, NumPy array, shape of (2204 - X, 2)
    sGenePairs_tr_lG -- training data gene pairs, list, length of (X)
    sGenePairs_tst_lG -- test data gene pairs, list, length of (2204 - X)
    """
    
    import numpy as np
    # Find unique Transcription Factors (TFs) in the gene pairs data
    sTFs_lG = [ sGenePairs_lG[ iX ][ 0 ] for iX in range( len( sGenePairs_lG ) ) ] 
    sUniqueTFs_nG = np.unique( sTFs_lG )

    # Sort gene pairs based on TFs
    sSortedGenePairs_nG = sort_Gene_Pairs( sGenePairs_lG )

    # Select the first TFs and find all indexes where TF located
    randomTF = sUniqueTFs_nG[ 0 ]
    idX_TF = np.where( sSortedGenePairs_nG == randomTF )[ 0 ]
    
    # New array to store which genes to be assigned to the training and testing data 
    ## True: training data, False: testing data
    bGeneChecks_nG = np.empty( shape = ( len( sGenePairs_lG ), 1), dtype = 'bool') # np array to check at which positions no gene is in the train data pairs
    bGeneChecks_nG[ : ] = True

    # Start with the first TF and its indexes
    iP = 0
    nLocations_lG = idX_TF.tolist()
    sGenesList_lG = []
    sGenesList_lG.append( randomTF )
    iTF = 0
    
    while True: 
    
        # Find and store targets which selected gene affects and TFs that influence the selected gene 
        sTargets_lG = []
        sTargets_0_lG = [ sSortedGenePairs_nG[ iX ][ 0 ] for iX in idX_TF ] # TFs of selected pairs 
        sTargets_1_lG = [ sSortedGenePairs_nG[ iX ][ 1 ] for iX in idX_TF ] # target genes of selected pairs
    
        sTargets_lG.extend( sTargets_0_lG ) 
        sTargets_lG.extend( sTargets_1_lG )
        sTargets_lG.remove( randomTF )
    
        # Store all genes being found as a list
        sGenesList_lG.extend( sTargets_lG )
    
        # List for temporary targets and their indexes
        nTempIDX_lG = []
        sNewTargets_lG = []
    
        # Find genes and indexes which are influenced by/influence genes in the target list 
        for iT in sTargets_lG:
        
            nNewTargetIDX = np.where( sSortedGenePairs_nG == iT )[ 0 ]
            nTempIDX_lG.extend( nNewTargetIDX.tolist() )
            nLocations_lG.extend( nNewTargetIDX.tolist() )
        
            sTempTargets_0_lG = [ sSortedGenePairs_nG[ iNT ][ 0 ] for iNT in nNewTargetIDX ]
            sTempTargets_1_lG = [ sSortedGenePairs_nG[ iNT ][ 1 ] for iNT in nNewTargetIDX ]
       
            sNewTargets_lG.extend( sTempTargets_0_lG )
            sNewTargets_lG.extend( sTempTargets_1_lG )
            sNewTargets_lG.remove( iT )

        # Unique locations of all genes being stored 
        nUniqueLocations_nG = np.unique( nLocations_lG )

        # Unique new targets 
        sNewUniqueTargets_nG = np.unique( sNewTargets_lG )
        sNewUniqueTargets_nG = sNewUniqueTargets_nG.tolist()
    
        # Boolean to check if new targets are already in the genesList 
        bExists = all( [ iG in sGenesList_lG for iG in sNewUniqueTargets_nG ] )
    
        # If all genes exist, assign them False to put into the testing data 
        if bExists:
    
            bGeneChecks_nG[nUniqueLocations_nG] = False
        
            nPairs = len( nUniqueLocations_nG )
            iP += nPairs
    
            iTF += 1 
        
            if iTF == len(sUniqueTFs_nG):
                break
            
            randomTF = sUniqueTFs_nG[ iTF ]
            idX_TF = np.where( sSortedGenePairs_nG == randomTF )[ 0 ]
            nLocations_lG = idX_TF.tolist()
            sGenesList_lG.clear()
            sGenesList_lG.append( randomTF )
        
        else:
        
            # If the number of unique locations > uLocMax, we need to continue with the next TF to search new pairs
            if len( nUniqueLocations_nG ) > uLocMax:
                
                iTF += 1
                
                if iTF == len( sUniqueTFs_nG ):
                    break
                    
                randomTF = sUniqueTFs_nG[ iTF ]
                idX_TF = np.where( sSortedGenePairs_nG == randomTF )[ 0 ]
                nLocations_lG = idX_TF.tolist()
                sGenesList_lG.clear()
                sGenesList_lG.append( randomTF )
            
            else:
                idX_TF = nTempIDX_lG
                sGenesList_lG.extend( sNewUniqueTargets_nG )
    
    # Retrieve positions of training and testing data
    iPos_tr_X = np.where( bGeneChecks_nG == True )[ 0 ]
    iPos_tst_X = np.where( bGeneChecks_nG == False )[ 0 ]

    sGenePairs_tr_nG = sSortedGenePairs_nG[ iPos_tr_X ]
    sGenePairs_tst_nG = sSortedGenePairs_nG[ iPos_tst_X ]
    
    sGenePairs_tr_lG = sGenePairs_tr_nG.tolist()
    sGenePairs_tst_lG = sGenePairs_tst_nG.tolist()

    return sGenePairs_tr_nG, sGenePairs_tst_nG, sGenePairs_tr_lG, sGenePairs_tst_lG

def print2F( s ) :
    # Prints output "s", but also indicating the function producing this print output
    import inspect
    print( '- ' + inspect.stack()[1][3] + ' = ' + s )