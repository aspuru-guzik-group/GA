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

Note: The discriminator scores were not accounted for in the fitness function for this experiment.  
The main fitness function (which incorporates QED scores) is defined in line 157 of file 'generation_props.py'


# How are the results saved?  : 
Note: 'i' is the run iteration 
1. results_0_i:  
   Each sub-directory is named by the generation. The smile strings (ordered by fitness) and corresponding molecular properties are provided as text
   files: 'smiles_ordered.txt', 'logP_ordered.txt', 'sas_ordered.txt', 'ringP_ordered.txt', 'discrP_ordered.txt'. 
   Outside the sub-directories is the information about the best molecules of a generation. 
2. runs:  
   Live plotting tool (for quick visualizations while the experiment was running).  
   To activate, run: tensorboard --logdir ./runs (after cd'ing in the current directory)