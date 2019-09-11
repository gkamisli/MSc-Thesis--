# Network creation functions for the network

def dropoutWrapper(n_hidden, dropout):
    
    """
    Function that creates a LSTMCell with DropoutWrapper
    
    Arguments: 
    n_hidden -- the number of nodes in the hidden layer
    dropout -- the probability of dropping nodes
    
    Returns:
    dropout_cell -- the LSTMCell after a DropoutWrapper
    
    """
    import tensorflow as tf
    
    cell = tf.nn.rnn_cell.LSTMCell( n_hidden, state_is_tuple = True )
    dropout_cell = tf.nn.rnn_cell.DropoutWrapper( cell, output_keep_prob = ( 1 - dropout ) )
    
    return dropout_cell

def LSTM_Model(X, init_state, n_layer, n_hidden, dropout):
    
    """
    Function to create a multi-layer LSTM network with static unrolling.  
    - tf.contrib.static_rnn accepts:
      1) Cell = cell of the network, created by Multi_RNN, list of length (n_layer) 
      2) Inputs = inputs with shape of (time_steps x (batch_size, n_input))
      3) Initial_state = initial states as tuple (cell, hidden) with a length of (time_steps)
    Hence, when X and init_state arguments are given, a transformation is applied to both to prepare for the static_rnn
    
    Arguments:
    X -- input for the network, NumPy array, shape of (iXnum, iNnum + 1, iSnum) 
    init_state -- initial state for the network, NumPy array, shape of (n_layer, 2, batch_size, n_hidden)
    n_layer -- the number of layers 
    n_hidden -- the number of hidden neurons in LSTM cells
    dropout -- the probability of dropping nodes    
    
    Returns:
    output -- the output from the network, NumPy array, shape of (batch_size, time_steps, n_hidden)
    state -- the state from the network, tuple (cell, hidden), length of (time_steps)
    
    """
  
    import tensorflow as tf
    from tensorflow.contrib import rnn
    
    time_steps = X.shape[1]
   
    # Encoding layer
    X = tf.transpose( X, perm = [ 1, 0, 2 ] )
    X = tf.reshape( X, [ -1, X.shape[2] ] )
    
    # Split data for the rnn cells, new shape (time_steps x (batch_size, n_input))
    X = tf.split( X, time_steps, axis = 0 )

    # Unpacking the initial state into hidden and cell states
    state_per_layer_list = tf.unstack( init_state, axis = 0 )
    rnn_tuple_state = tuple(
        [ rnn.LSTMStateTuple( state_per_layer_list[idx][0], state_per_layer_list[idx][1] )
        for idx in range( n_layer ) ] )

    # Define LSTM cells of gene
    cell_of_gene = tf.nn.rnn_cell.MultiRNNCell( [ dropoutWrapper( n_hidden, dropout ) 
                                     for _ in range( n_layer ) ], state_is_tuple=True )

    
    output, state = rnn.static_rnn( cell = cell_of_gene, 
                                   inputs = X, 
                                   initial_state = rnn_tuple_state )
    
    return output, state

def length(sequence):
    
    """
    Function to retrieve sequence length of each example 
    
    Arguments:
    sequence (i.e., X) -- input of the network, NumPy array, shape of (iXnum, iNnum + 1, iSnum) 
    
    Returns:
    sequence_length -- the non_padding sequence length for each example, list of length (iXnum)
    
    """
    import tensorflow as tf
    
    # Collapse last dimension into scalars and 
    # Create a binary mask for sequence values: assign 1s to non_padding frames; 0s to padding frames
    binary_mask = tf.sign( tf.reduce_max( tf.abs( sequence ), 2 ) )
    
    # Retrieve the sequence length by summing the number of 1s in binary_mask
    sequence_length = tf.reduce_sum( binary_mask, 1 )
    sequence_length = tf.cast( sequence_length, tf.int32 )
    
    return sequence_length

def dynamicLSTM(X, n_layer, n_hidden, dropout):
    
    """
    Function to create a multi-layer LSTM network with dynamic unrolling which reads sequence until paddings
    - tf.nn.dynamic_rnn accepts
      1) Cell = Cell of the network, created by Multi_RNN, list of length (n_layer) 
      2) Inputs = Inputs with shape of (batch_size, time_steps x n_input))
      3) Sequence_length = Length information of every example, retrieved by length(X) function, list of length (X.shape[0])
    For dynamic_rnn, we do not need to indicate an initial state as it retrieves the initial states with zeros whenever a cell is     created. Additionally, no transformation is not needed for the input.
    
    Arguments:
    X -- input for the network, NumPy array, shape of (iXnum, iNnum + 1, iSnum) 
    n_layer -- the number of layers 
    n_hidden -- the number of hidden neurons in LSTM cells
    dropout -- the dropout rate for DropoutWrapper
    
    Returns:
    output -- the output from the network, NumPy array, shape of (batch_size, time_steps, n_hidden)
    state -- the state from the network, tuple (cell, hidden), length of (time_steps)
    
    """
    import tensorflow as tf
    sequence_length = length( X )
    
    # Create cells for every layer with dropoutWrapper
    cells = tf.nn.rnn_cell.MultiRNNCell( [ dropoutWrapper( n_hidden, dropout ) 
                                     for _ in range( n_layer ) ], state_is_tuple=True )

    output, state = tf.nn.dynamic_rnn( cells,
                                      inputs=X,
                                      dtype=tf.float32,
                                      sequence_length=sequence_length )
    return output, state

def dynamicLSTM_Attention(X, n_layer, n_hidden, dropout):
    
    """
    Function to create a multi-layer LSTM network with dynamic unrolling which reads whole sequence (including paddings)
    - tf.nn.dynamic_rnn accepts
      1) Cell = Cell of the network, created by Multi_RNN, list of length (n_layer) 
      2) Inputs = Inputs with shape of (batch_size, time_steps x n_input))
    For dynamic_rnn, we do not need to indicate an initial state as it retrieves the initial states with zeros whenever a cell is     created. Additionally, no transformation is not needed for the input.
    
    Arguments:
    X -- input for the network, NumPy array, shape of (iXnum, iNnum + 1, iSnum) 
    n_layer -- the number of layers 
    n_hidden -- the number of hidden neurons in LSTM cells
    dropout -- the dropout rate for DropoutWrapper
    
    Returns:
    output -- the output from the network, NumPy array, shape of (batch_size, time_steps, n_hidden)
    state -- the state from the network, tuple (cell, hidden), length of (time_steps)
    
    """
    
    import tensorflow as tf
    cells = tf.nn.rnn_cell.MultiRNNCell( [ dropoutWrapper( n_hidden, dropout ) 
                                     for _ in range( n_layer ) ], state_is_tuple=True )
    
    output, state = tf.nn.dynamic_rnn(cells,
                                      inputs=X,
                                      dtype=tf.float32)
    return output, state

def last_relevant(output, sequence_length):
    
    """
    Function to retrieve the last relevant output created by dynamicLSTM function with variable length sequence
    
    Arguments:
    output -- the output produced by the network, NumPy array, shape of (batch_size, time_steps, n_hidden)
    
    Returns:
    last_relevant_output -- NumPy array, shape of (batch_size, n_hidden)
  
    """
    import tensorflow as tf
    
    # Get necessary information from output
    # 1) Number of examples (batch_size for training data, number of examples for testing data)
    # 2) Maximum length of the input, time_steps
    # 3) Output size, n_hidden
    num_examples = tf.shape( output )[ 0 ]
    time_steps = tf.shape( output )[ 1 ]
    n_hidden = int( output.get_shape()[ 2 ] )
    
    # Create an index to keep track of last index where last non-padding frames are located 
    last_index = tf.range( 0, num_examples ) * time_steps + ( sequence_length - 1 )
    
    # Flatten the output (batch_size x time_steps, n_hidden) to use last_index 
    flat = tf.reshape( output, [ -1, n_hidden ] )
    last_relevant_output = tf.gather( flat, last_index )
    
    return last_relevant_output

def attention(state, output, n_hidden):
    
    """
    Function to create an attention mechanism, introduced by Bahdanau
    
    Arguments:
    state -- the state from the network, tuple (cell, hidden), length of (time_steps)
    output -- the output produced by the network, NumPy array, shape of (batch_size, time_steps, n_hidden)
    n_hidden -- the number of hidden neurons in LSTM cells
    
    Returns:
    context_vector -- encodings from the mechanism, NumPy array, shape of (batch_size, n_hidden)
    attention_weights -- attention weights of the mechanism, NumPy array, shape of (batch_size, time_steps)
    
    """
    import tensorflow as tf
    
    # Add one more axis for time_step
    state_with_time_axis = tf.expand_dims( state, 1 )
    
    # Pass output and state through a tanh layer after a fully-connected (dense) layer
    score = tf.layers.dense(tf.nn.tanh(tf.layers.dense( output, n_hidden ) + 
                               tf.layers.dense( state_with_time_axis, n_hidden ) ), 1 )                  

    # Attention weight and context vector calculations
    attention_weights = tf.nn.softmax( score, axis=1 )                              
    context_vector = attention_weights * output
    context_vector = tf.reduce_sum( context_vector, axis=1 )
    
    return context_vector, attention_weights

def batch_normalised_sum(output_A, output_B):
    
    """
    Function to retrieve last time step's output from gene A and gene B, then normalise the output representations
    and add values to be evaluated by the Static network together for predictions
    
    Arguments:
    output_A -- the output produced by the Static network for gene A, list, length of (time_steps)
    output_B -- the output produced by the Static network for gene B, list, length of (time_steps)
    
    Returns: 
    output -- the normalised and added outputs from the networks of gene A and gene B, NumPy array, 
                                                                                    shape of (batch_size, n_hidden)
    """
    
    import tensorflow as tf
    
    # Normalised outputs from last time step
    last_outputA = tf.layers.batch_normalization( output_A[ -1 ] )
    last_outputB = tf.layers.batch_normalization( output_B[ -1 ] )

    # Add representations from outputs
    output = tf.math.add( last_outputA, last_outputB )
    
    return output