# How to run the code? : 
Code for the unconstrained optimization experiment can be run using:  
```
python ./core_GA.py
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
- beta_preference: Different beta values with which the experiment was run. For each beta, the experiment is repeated 5 times. 


# How are the results saved?  : 
All the results are savents in the 'results' directory. Our results are saved as (Note: 'i' is the run iteration): 
1. images_generation_beta_i:  
   Images of the top 100 molecules of each generation. Below each molecule are the Fitness, logP, SA, ring penalty and discriminator scores
2. results_beta_i:  
   Each sub-directory is named by the generation. The smile strings (ordered by fitness) and corresponding molecular properties are provided as text
   files: 'smiles_ordered.txt', 'logP_ordered.txt', 'sas_ordered.txt', 'ringP_ordered.txt', 'discrP_ordered.txt'. 
   Outside the sub-directories is the information about the best molecules of a generation. 
3. saved_models_beta_i:  
   The trained discriminators after each generation. Please Note: We did not make use of the discriminator predictions in the Fitness for this experiment (beta is set to 0).