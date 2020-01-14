import os
from selfies import decoder 
import time
import multiprocessing
import torch
from tensorboardX import SummaryWriter
from selfies import encoder
      
#### Directory Imports
import discriminator as D
import evolution_functions as evo
import generation_props as gen_func



def initiate_ga(num_generations,            generation_size,    starting_selfies,max_molecules_len,
                disc_epochs_per_generation, disc_enc_type,      disc_layers,     training_start_gen,           
                device,                     properties_calc_ls, num_processors,  beta, starting_smile, desired_delta, save_curve):
    
       
    
    # Obtain starting molecule
    starting_smiles = evo.sanitize_multiple_smiles([decoder(selfie) for selfie in starting_selfies])

    
    # Recording Collective results
    smiles_all         = []    # all SMILES seen in all generations
    selfies_all        = []    # all SELFIES seen in all generation
    smiles_all_counter = {}    # 
    
    
    # Initialize a Discriminator
    discriminator, d_optimizer, d_loss_func = D.obtain_initial_discriminator(disc_enc_type, disc_layers, max_molecules_len, device)
    
    # Read in the Zinc data set 
    molecules_reference = evo.read_dataset_encoding(disc_enc_type)
#    print(molecules_reference[0:2])
#    raise Exception()    
    molecules_reference = dict.fromkeys(molecules_reference, '') # convert the zinc data set into a dictionary

    # Set up Generation Loop 
    total_time = time.time()
    for generation_index in range(1, num_generations+1):
        print("   ###   On generation %i of %i"%(generation_index, num_generations))
        start_time = time.time()

              
        # Obtain molecules from the previous generation 
        smiles_here, selfies_here = gen_func.obtain_previous_gen_mol(starting_smiles,   starting_selfies, generation_size, 
                                                                     generation_index,  selfies_all,      smiles_all)

        # Calculate fitness of previous generation (shape: (generation_size, ))
        fitness_here, order, fitness_ordered, smiles_ordered, selfies_ordered = gen_func.obtain_fitness(disc_enc_type,      smiles_here,   selfies_here, 
                                                                                                        properties_calc_ls, discriminator, generation_index,
                                                                                                        max_molecules_len,  device,        generation_size,  
                                                                                                        num_processors,     writer,        beta,            image_dir, data_dir, starting_smile, desired_delta, save_curve)

        # Obtain molecules that need to be replaced & kept
        to_replace, to_keep = gen_func.apply_generation_cutoff(order, generation_size)
        
        # Obtain new generation of molecules 
        smiles_mutated, selfies_mutated = gen_func.obtain_next_gen_molecules(order,           to_replace,     to_keep, 
                                                                             selfies_ordered, smiles_ordered, max_molecules_len)
        for item in smiles_mutated:
            mol, smi_canon, did_convert = evo.sanitize_smiles(item)
            if did_convert == False:
                raise Exception('Failed with (2): ', item)
                
        # Record in collective list of molecules 
        smiles_all, selfies_all, smiles_all_counter = gen_func.update_gen_res(smiles_all, smiles_mutated, selfies_all, selfies_mutated, smiles_all_counter)

        print('Generation time: ', round((time.time()-start_time), 2), ' seconds')

    print('Total time: ', round((time.time()-total_time)/60, 2), ' mins')
    print('Total number of unique molecules: ', len(smiles_all_counter))
    return smiles_all_counter



if __name__ == '__main__':
    
    file_name      = './datasets/delta40/dearomatization_done.txt'
    with open(file_name) as f:
        starting_smile = f.readlines()
        starting_smile = [x.strip() for x in starting_smile] 

    for smile in starting_smile:
        print('Initiating: ', smile)
        save_curve = []
        beta_preference = [0]
        results_dir = evo.make_clean_results_dir()
        
        exper_time = time.time()
        for i in range(1):
            for beta in beta_preference:
                
                image_dir, saved_models_dir, data_dir = evo.make_clean_directories(beta, results_dir, i) # clear directories 
                
                # Initialize new TensorBoard writers
                torch.cuda.empty_cache()
                writer = SummaryWriter()   
        
                # Initiate the Genetic Algorithm
                smiles_all_counter = initiate_ga(    num_generations            = 20,
                                                     generation_size            = 500,
                                                     starting_selfies           = [encoder(smile)],       
                                                     max_molecules_len          = 81,
                                                     disc_epochs_per_generation = 0,
                                                     disc_enc_type              = 'properties_rdkit',                             # 'selfies' or 'smiles' or 'properties_rdkit'
                                                     disc_layers                = [100, 10],
                                                     training_start_gen         = 200,                                            # generation index to start training discriminator
                                                     device                     = 'cpu', 
                                                     properties_calc_ls         = ['logP', 'SAS', 'RingP', 'SIMILR'],             # None: No properties ; 'logP', 'SAS', 'RingP'
                                                     num_processors             = multiprocessing.cpu_count(),
                                                     beta                       = beta,
                                                     starting_smile             = smile,
                                                     desired_delta              = 0.4 ,
                                                     save_curve                 = save_curve
                                                )   
        print('Total Experiment time: ', (time.time()-exper_time)/60, ' mins')
        
        # SAVE THE AMOUNT OF IMPROVEMENT: 
        f=open("IMPROVEMENT.txt", "a+")
        A = save_curve[1:]
        if max(A) < -100:
            f.write('Failed improvement {} \n'.format(starting_smile))
        else:
            improvement = max(A) - save_curve[0]
            print('IMPROVEMENT: ', improvement)
            f.write('{} \n'.format(improvement))
            f.close()

    

        
        


        
