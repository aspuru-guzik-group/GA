# How to run the code? : 
Code for the unconstrained optimization experiment can be run uing:  
```python ./core_GA.py
```  

The following settings can be used (found at the end of the file): 
- num_generations: Number of generations to run the GA
- generation_size: Molecular population size encountered in each generation 
- starting_selfies: Initial population of molecules 
- max_molecules_len: Length of the largest molecule string
- disc_epochs_per_generation: Number of epochs of training the discriminator neural network 
- disc_enc_type: Type of molecular encoding shown to the discriminator
- disc_layers : Discriminator architecture
- training_start_gen: generation after which discriminator training begins 
- device: Device the discriminator is trained on 
- properties_calc_ls: Property evaluations to be completed for each molecule of the GA
- num_processors: Number of cpu cores to parallelize calculations over
- beta: Value of parameter beta

# How are the results saved?  : 


