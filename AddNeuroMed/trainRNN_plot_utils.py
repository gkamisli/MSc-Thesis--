# Plot functions for the network

def plot_inputs(f1_score_dict, train_loss_dict, test_loss_dict, n_epoch, hyperparameter):
    
    """
    Function that plot the results of the neural network algorithm; train and test losses, and f1_score for different values
    of a hyperparameter
    
    Arguments:
    f1_score_dict -- F1 scores for test data, dictionary, keys (hyperparameter values)
    train_loss_dict -- train losses, dictionary, keys (hyperparameter values)
    test_loss_dict -- test losses, dictionary, keys (hyperparameter values) 
    n_epoch -- number of epochs
    hyperparameter -- dropout/learning rate
    
    """
    import pandas as pd
    import matplotlib.pyplot as plt
    
    keys = list(f1_score_dict.keys())
    x = range(n_epoch)
    
    width = 20
    height = 5
    plt.figure( figsize = ( width, height ) )
    
    df_f1_score = pd.DataFrame({keys[0]: f1_score_dict[keys[0]],
                                keys[1]: f1_score_dict[keys[1]], 
                                keys[2]: f1_score_dict[keys[2]], 
                                keys[3]: f1_score_dict[keys[3]], 
                                keys[4]: f1_score_dict[keys[4]] } )
    
    df_train_loss = pd.DataFrame({keys[0]: train_loss_dict[keys[0]],
                                  keys[1]: train_loss_dict[keys[1]], 
                                  keys[2]: train_loss_dict[keys[2]], 
                                  keys[3]: train_loss_dict[keys[3]], 
                                  keys[4]: train_loss_dict[keys[4]] } )
    
    df_test_loss = pd.DataFrame({keys[0]: test_loss_dict[keys[0]],
                                 keys[1]: test_loss_dict[keys[1]], 
                                 keys[2]: test_loss_dict[keys[2]], 
                                 keys[3]: test_loss_dict[keys[3]], 
                                 keys[4]: test_loss_dict[keys[4]] } )

    plt.subplot(1, 3, 1)
    plt.plot( x, df_f1_score[keys[0]], marker='', color='skyblue', linewidth=2, label=keys[0])
    plt.plot( x, df_f1_score[keys[1]], marker='', color='green', linewidth=2, label=keys[1])
    plt.plot( x, df_f1_score[keys[2]], marker='', color='orange', linewidth=2, label=keys[2])
    plt.plot( x, df_f1_score[keys[3]], marker='', color='red', linewidth=2, label=keys[3])
    plt.plot( x, df_f1_score[keys[4]], marker='', color='purple', linewidth=2, label=keys[4])
    plt.xlabel( "Epochs" )
    plt.title("F1 score vs. {}".format(hyperparameter))
    plt.legend()
    
    plt.subplot(1, 3, 2)
    plt.plot( x, df_train_loss[keys[0]], marker='', color='skyblue', linewidth=2, label=keys[0])
    plt.plot( x, df_train_loss[keys[1]], marker='', color='green', linewidth=2, label=keys[1])
    plt.plot( x, df_train_loss[keys[2]], marker='', color='orange', linewidth=2, label=keys[2])
    plt.plot( x, df_train_loss[keys[3]], marker='', color='red', linewidth=2, label=keys[3])
    plt.plot( x, df_train_loss[keys[4]], marker='', color='purple', linewidth=2, label=keys[4])
    plt.xlabel( "Epochs" )
    plt.title("Train losses vs. {}".format(hyperparameter))
    plt.legend()
    
    plt.subplot(1, 3, 3)
    plt.plot( x, df_test_loss[keys[0]], marker='', color='skyblue', linewidth=2, label=keys[0])
    plt.plot( x, df_test_loss[keys[1]], marker='', color='green', linewidth=2, label=keys[1])
    plt.plot( x, df_test_loss[keys[2]], marker='', color='orange', linewidth=2, label=keys[2])
    plt.plot( x, df_test_loss[keys[3]], marker='', color='red', linewidth=2, label=keys[3])
    plt.plot( x, df_test_loss[keys[4]], marker='', color='purple', linewidth=2, label=keys[4])
    plt.xlabel( "Epochs" )
    plt.title("Test losses vs. {}".format(hyperparameter))
    plt.legend()
    plt.show()
    
def plot_one_input(f1_score_dict, train_loss_dict, test_loss_dict, n_epoch):
    
    """
    Function that plot the results of the neural network algorithm; train and test losses, and f1_score 
    
    Arguments:
    f1_score_dict -- F1 scores for test data, dictionary, keys (hyperparameter values)
    train_loss_dict -- train losses, dictionary, keys (hyperparameter values)
    test_loss_dict -- test losses, dictionary, keys (hyperparameter values) 
    n_epoch -- number of epochs
    
    """

    import pandas as pd
    import matplotlib.pyplot as plt
    
    keys = list(f1_score_dict.keys())
    x = range(n_epoch)

    width = 20
    height = 5
    plt.figure( figsize = ( width, height ) )

    df_f1_score = pd.DataFrame({keys[0]: f1_score_dict[keys[0]],
                                } )

    df_train_loss = pd.DataFrame({keys[0]: train_loss_dict[keys[0]],
                                  } )

    df_test_loss = pd.DataFrame({keys[0]: test_loss_dict[keys[0]],
                                 } )

    plt.subplot(1, 2, 1)
    plt.plot( x, df_f1_score[keys[0]], marker='', color='purple', linewidth=2, label=keys[0])
    plt.xlabel( "Epochs" )
    plt.title("F1 scores")

    plt.subplot(1, 2, 2)
    plt.plot( x, df_train_loss[keys[0]], marker='', color='red', linewidth=2, label="Train loss")
    plt.plot( x, df_test_loss[keys[0]], marker='', color='blue', linewidth=2, label="Test loss")
    plt.xlabel( "Epochs" )
    plt.title("Train and Test losses ")
    plt.legend()
    plt.show()