**Data Sets**
<br> 1) AddNeuroMed
<br> 2) HBTRC (Harvard Brain Tissue Resource Centre)


**Models**
<br> 1) Static: Read inputs with 0 paddings and create a static network with 101 time steps
<br> 2) Encoder: Read inputs with variable sequence lengths and create a dynamic network
<br> 3) Attention: Read inputs with 0 paddings and add an attention layer for 101 attention size

**File types**
<br> 1) loadData.ipynb: Load raw data and transform into NumPy arrays
<br> 2) preprocessData.ipynb: Process NumPy arrays and generate examples for each class (+1: gene A regulates gene B, -1: gene A is regulated by gene B, 0: no interaction between genes)
<br> 3) trainRNN.ipynb: Create and train the network based on three models with respect to data sets
