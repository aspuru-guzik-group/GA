#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: akshat

TODO: 
    - Benzene inside the random string
    - Number after a 'Ring' inside the molecule
    
"""
from selfies import encoder, decoder
import numpy as np
import random
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi
from rdkit.Chem import Descriptors
from SAS_calculator.sascorer import calculateScore


def sanitize_smiles(smi):
    '''Return a canonical smile representation of smi
    
    Parameters:
    smi (string) : smile string to be canonicalized 
    
    Returns:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object                          (None if invalid smile string smi)
    smi_canon (string)          : Canonicalized smile representation of smi (None if invalid smile string smi)
    conversion_successful (bool): True/False to indicate if conversion was  successful 
    '''
    try:
        mol = smi2mol(smi, sanitize=True)
        smi_canon = mol2smi(mol, isomericSmiles=False, canonical=True)
        return (mol, smi_canon, True)
    except:
        return (None, None, False)
         
    
def get_logP(mol):
    '''Calculate logP of a molecule 
    
    Parameters:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object, for which logP is to calculates
    
    Returns:
    float : logP of molecule (mol)
    '''
    return Descriptors.MolLogP(mol)


def get_SA(mol):
    return calculateScore(mol)

def calc_RingP(mol):
    '''Calculate Ring penalty for each molecule in unseen_smile_ls,
       results are recorded in locked dictionary props_collect 
    '''
    cycle_list = mol.GetRingInfo().AtomRings() 
    if len(cycle_list) == 0:
        cycle_length = 0
    else:
        cycle_length = max([ len(j) for j in cycle_list ])
    if cycle_length <= 6:
        cycle_length = 0
    else:
        cycle_length = cycle_length - 6
    return cycle_length


def get_random_selfie_mol():
    '''
    Create random molecules using the SELFIES alphabet. Strings are contrained to 
    a length of 81.
    '''
    valid = False
    alphabet = ['[Branch1_1]', '[Branch1_2]','[Branch1_3]', '[epsilon]', '[Ring1]', '[Ring2]', '[Branch2_1]', '[Branch2_2]', '[Branch2_3]', '[F]', '[O]', '[=O]', '[N]', '[=N]', '[#N]', '[C]', '[=C]', '[#C]', '[S]', '[=S]', '[C][=C][C][=C][C][=C][Ring1][Branch1_1]']

    while valid != True:
        selfie_char_ls = []
        for i in range(81):
            random_char = np.random.choice(alphabet, size=1)[0]
            selfie_char_ls.append(random_char)
        
        selfie_str = ''.join(x for x in selfie_char_ls)
        decoded_smile_str = decoder(selfie_str)
        if decoded_smile_str != -1:
            valid = True
        # check if smile string is recognized by rdkit
        mol, smiles_canon, done = sanitize_smiles(decoded_smile_str)
        if mol == None or smiles_canon == '' or len(smiles_canon) > 81: 
            valid = False
            continue

    return smiles_canon





num_runs = 10
max_scores = []
for i in range(num_runs): 

    A = []
    for j in range(50000):
        A.append(get_random_selfie_mol())

        if len(A[j]) > 81:
            raise Exception('Length fail!')

    logP_scores   = []
    SA_scores    = []
    RingP_scores =  []
    
    # CODE TEST....
    print('Looking at i: ', i)
    for item in A:
        mol, smiles_canon, done = sanitize_smiles(item)
        if mol == None or done == False:
            raise Exception('A molecule is incorrect! Test Failed')

        logP_scores.append(get_logP(mol))
        SA_scores.append(get_SA(mol))
        RingP_scores.append(calc_RingP(mol))
        
        # Save all the smile codes of this data set
        f = open("./results/results_{}.txt".format(i), "a+")
        f.write('smile: {} \n'.format(smiles_canon))
        f.close()
        
    # Calculate J(m)
    logP_norm  = [((x - 2.4729421499641497) / 1.4157879815362406) for x in logP_scores]
    SAS_norm   = [((x - 3.0470797085649894) / 0.830643172314514) for x in SA_scores]
    RingP_norm = [((x - 0.038131530820234766) / 0.2240274735210179) for x in RingP_scores]
    J = []
    for i in range(len(logP_norm)):
        J.append(logP_norm[i] - SAS_norm[i] - RingP_norm[i])
        

        
    print('smile: ', A[J.index(max(J))], max(J))
    #max_scores.append(max(J))
    
    # Save result in text file
    f = open("results.txt", "a+")
    f.write('smile: {}, J: {} \n'.format(A[J.index(max(J))], max(J)))
    f.close()


    
    
    
    
