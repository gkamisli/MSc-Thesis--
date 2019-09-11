# Utility functions for the network 

def shuffle_classes(X_A, X_B):
    
    """
    Function that takes input and shuffle samples within each class.
    
    Arguments:
    X_A -- array of the SNPs and RNA levels of gene A, NumPy array, shape of (iXnum, iNnum + 1, iSnum)
    X_B -- array of the SNPs and RNA levels of gene B, NumPy array, shape of (iXnum, iNnum + 1, iSnum)
    
    Returns:
    X_A -- array of shuffled X_A, NumPy array, shape of (iXnum, iNnum + 1, iSnum)
    X_B -- array of shuffled X_B, NumPy array, shape of (iXnum, iNnum + 1, iSnum)
    
    """
    
    import numpy as np
    
    # Retrieve the number of examples per class
    num_examples = int( X_A.shape[ 0 ] / 3 )

    # Split each class of X_A based on iXnum
    first_group_A = X_A[ 0 : num_examples, :, : ]
    second_group_A = X_A[ num_examples : 2 * num_examples, :, : ]
    third_group_A = X_A[ 2 * num_examples : , :, : ]

    # Split each class of X_B based on iXnum
    first_group_B = X_B[ 0 : num_examples, :, : ]
    second_group_B = X_B[ num_examples : 2 * num_examples, :, : ]
    third_group_B = X_B[ 2 * num_examples : , :, : ]

    # Shuffle every class samples within itself based on gene A's samples  
    iXs_shuffled_first_group = np.random.permutation( first_group_A.shape[0] )
    iXs_shuffled_second_group = np.random.permutation( second_group_A.shape[0] )
    iXs_shuffled_third_group = np.random.permutation( third_group_A.shape[0] )

    # Gene A classes after shuffling
    first_group_A = first_group_A[iXs_shuffled_first_group, :, :]
    second_group_A = second_group_A[iXs_shuffled_second_group, :, :]
    third_group_A = third_group_A[iXs_shuffled_second_group, :, :]

    # Gene B classes after shuffling
    first_group_B = first_group_B[iXs_shuffled_first_group, :, :]
    second_group_B = second_group_B[iXs_shuffled_second_group, :, :]
    third_group_B = third_group_B[iXs_shuffled_second_group, :, :]

    # Combine shuffled arrays 
    X_A = np.concatenate( ( first_group_A, second_group_A, third_group_A ), axis = 0 )
    X_B = np.concatenate( ( first_group_B, second_group_B, third_group_B ), axis = 0 )
    
    return X_A, X_B

def input_reshape(rSnp_nXSN, rRna_nXS):
    
    """
    Function that takes two inputs, reshape and convert them to a single input for the RNN model 
    
    Arguments:
    rSnp_nXSN -- array of the SNP encodings of gene G, NumPy array, shape of (iXnum, iSnum, iNnum)
    rRna_nXS -- array of the gene expression levels of gene G, NumPy array, shape of (iXnum, iSnum)
    
    Returns:
    rSnpRna_nXNS -- array of SNP encodings of gene G and gene expression levels of gene G, 
                    NumPy array, shape of (iXnum, iNnum + 1, iSnum)
    
    """
    import numpy as np
    
    # Add one more axis to rRna_nXS
    rRna_nXNS = rRna_nXS[:, np.newaxis, :]
    
    # Swap axes 1 and 2 
    rSnp_nXNS = np.transpose( rSnp_nXSN, (0, 2, 1) )
    

    # Merge rSnp_nXNS and rRnaB_nXNS over SNP axes
    rSnpRna_nXNS = np.hstack( ( rRna_nXNS, rSnp_nXNS ) )
    
    return rSnpRna_nXNS

def extract_batch_size(X, step, batch_size):
    
    """
    Function to fetch a batch_size amount of data from training data. It takes equal (batch_size/3) amount of
    data from each class to feed the network with balanced data
    
    Arguments:
    X -- input array, shape of (iXnum, iNnum + 1, iSnum)
    step -- integer to represent each step/batch
    batch_size -- integer amount to extract data
    
    Returns:
    batch_x -- NumPy array, shape of (batch_size, iNnum + 1, iSnum)
    
    """
    import numpy as np
    
    shape = list(X.shape)
    shape[0] = batch_size
    batch_x = np.empty(shape)
    extra = 0

    if (step * batch_size > X.shape[0]):
        step = 1
        
    for i in range(batch_size):
        
        index = ( ( step - 1 ) * batch_size + i) % len(X)
        
        if i % 3 == 0 and i != 0:
            extra += 1
            
        data_index = int( X.shape[ 0 ] / 3 ) * ( i % 3 ) + extra + ( step - 1 ) * int( batch_size / 3 )
        batch_x[i] = X[ data_index ]
        
    return batch_x

