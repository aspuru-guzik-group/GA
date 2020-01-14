import os
from sklearn.cluster import KMeans
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import colors
import joblib
import time

def smiles_alphabet(disc_enc_type):
   '''Return a list of characters present in the zinc dataset

   Parameters:
   disc_enc_type (string): Indicates whether to return SMILES/SELFiES characters

   Returns:
   alphabet: list of SELFIE/SMILE alphabets in Zinc
   '''
   if disc_enc_type == 'smiles':
       alphabet = ['C', 'c', 'H','O','o', 'N','n', 'S','s', 'F', 'P', 'I',
                   'Cl','Br', '=','#','(',')','[',']','1','2','3','4','5',
                   '6','7','8','9','+','-', '%', '0', 'X'] # SMILES Alphabets in zinc

   elif disc_enc_type == 'selfies':
       alphabet = ['[Ring1]',   '[Branch1_1]', '[Branch1_2]','[Branch1_3]', '[Cl]',
                   '[Ring2]',   '[Branch2_1]', '[Branch2_2]','[Branch2_3]', '[NH3+]',
                   '[N]',       '[=N]',        '[#N]',       '[C]',         '[=C]',
                   '[#C]',      '[S]',         '[=S]',       '[=O]',        '[Br]',
                   '[epsilon]', '[N+]',        '[NH+]',      '[NH2+]',      '[=NH+]',
                   '[=NH2+]',   '[I]',         '[O-]',       '[P]',         '[=P]',
                   '[S-]',      '[=N-]',       '[NH-]',      '[=O+]',       '[CH-]',
                   '[PH+]',     '[=S+]',       '[S+]',       '[CH2-]',      '[P+]',
                   '[O+]',      '[=N+]',       '[N-]' ,       '[=SH+]',     '[=OH+]',
                   '[#N+]',     '[=PH2]',      'X',           '[F]',        '[O]',  '[c]', '[-c]',
                  ] # SELFIES Alphabets in zinc
   else:
       exit('Invalid choice. Only possible choices are: smiles/selfies.')

   return alphabet

def to_onehot(molecule_str, disc_enc_type, max_molecules_len):
   '''Convert given molecule string into a one-hot encoding, with characters
      obtained from function 'smiles_alphabet'.

   One-hot encoding of arbitrary molecules is converted to len
   'max_molecules_len' by padding with character 'X'

   Parameters:
   molecule_str      (string): SMILE/SELFIE string of molecule
   disc_enc_type     (string): Indicating weather molecule string is either
                               SMILE or SELFIE
   max_molecules_len (string): Length of the one-hot encoding


   Returns:
   one_hots   (list of lists): One-Hot encoding of molecule string, padding
                               till length max_molecules_len (dim: len(alphabet) * max_molecules_len)
   '''
   alphabet = smiles_alphabet(disc_enc_type)
   alphabet_length = len(alphabet)
   one_hots=np.zeros(shape=(len(molecule_str), alphabet_length * max_molecules_len)).astype(np.int32)

   if disc_enc_type == 'smiles':
       alphabet.remove('Cl')      # Replace 'Cl' & 'Br' with 'Y' & 'Z' for convenience
       alphabet.remove('Br')      # (Searching for single characters is easier)
       alphabet.append('Y')
       alphabet.append('Z')

   for smi_index, smi in enumerate(molecule_str):
       # Relace 'Cl' and 'Br' with 'Y', 'Z' from smi (for conveninece)
       if disc_enc_type == 'smiles':
           smi = smi.replace('Cl', 'Y')
           smi = smi.replace('Br', 'Z')

       if disc_enc_type == 'selfies':
           smi = get_selfie_chars(smi)
       if len(smi) > max_molecules_len:
           raise Exception('smile is too large :): ', smi)
       for char_index, char in enumerate(smi):
           if char not in alphabet:
               print("smiles character %s not in alphabet MOLECULE: %s"%(char, smi))

           one_hots[smi_index][(char_index * alphabet_length) + alphabet.index(char)] = 1

   return (one_hots)



def collect_hist_data():
    time0=time.time()
    num_generations=1000
    smiles_all=[]
    smiles_per_generation=[]
    for g in range(1, num_generations+1):
        smiles = []
        for lineidx,line in enumerate(open("results_0_2/%i/smiles_ordered.txt"%(g), "r")):
            smiles.append(line.split()[0])
            if lineidx==50:
                break
        for smi in smiles:
            smiles_all.append(smi)


        smiles_per_generation.append(smiles)
    time1=time.time()
    print("read %i smiles codes after %.2f seconds"%(len(smiles_all), time1-time0))
    time0=time.time()
    #smiles_all = list(set(smiles_all))
    time1=time.time()
    print("found %i unique ones after %.2f seconds"%(len(smiles_all), time1-time0))
    return(smiles_all, smiles_per_generation)



smiles_all, smiles_per_generation = collect_hist_data()
print("found %i unique smiles"%(len(smiles_all)))
#smiles_all=smiles_all[:1000]


smiles_training = np.random.choice(smiles_all, size=2000)

print("selected %i unique smiles for training"%(len(smiles_training)))

n_clusters=20
rdkit_l=7
rdkit_s=4000
if not os.path.exists("kmeans.joblib"):
    time0=time.time()
    print("start with kmeans")
    print("convert to mol objects")
    mols=[Chem.MolFromSmiles(smi) for smi in smiles_training]
    print("convert to fingerprints")
    x_unscaled_fp = np.array([Chem.RDKFingerprint(mol, maxPath=rdkit_l, fpSize=rdkit_s) for mol in mols]).astype(float)
    x_unscaled_fp=np.array(x_unscaled_fp)
    print("converted to fingerprints:")
    print(x_unscaled_fp.shape)

    print("start with kmeans")
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    kmeans.fit(x_unscaled_fp)
    print(kmeans.labels_.shape)
    print(kmeans.cluster_centers_.shape)
    joblib.dump(kmeans, "kmeans.joblib")
    time1=time.time()
    print("kmeans done in %.1f seconds"%(time1-time0))
else:
    kmeans = joblib.load("kmeans.joblib")

repeats=5
if not os.path.exists("labels_sorted.dat"):
    time0=time.time()
    print("start labeling all molecules")
    #labels=[]
    labels_sorted=[]
    for g in range(len(smiles_per_generation)):
        print("get labels for generation %i if %i"%(g+1, len(smiles_per_generation)))
        fps = np.array([Chem.RDKFingerprint(Chem.MolFromSmiles(smi), maxPath=rdkit_l, fpSize=rdkit_s) for smi in smiles_per_generation[g]]).astype(float)
        labels_here = kmeans.predict(fps)

        labels_here=sorted(labels_here)
        for i in range(repeats):
            labels_sorted.append(labels_here)
        #if g==40:
        #    break
    labels_sorted=np.array(labels_sorted)
    np.savetxt("labels_sorted.dat", labels_sorted)
    time1=time.time()
    print("labels done in %.1f seconds"%(time1-time0))
else:
    labels_sorted=np.loadtxt("labels_sorted.dat")


print("start plotting")
print("labels:")
print(labels_sorted.shape)



plt.figure(figsize=(25,5))
img = plt.imshow(labels_sorted.T, interpolation='nearest', origin='lower', cmap="nipy_spectral", aspect="auto", vmin=-0.5, vmax = float(n_clusters)+0.2)
plt.xlabel("Generation (times %i)"%(repeats))
plt.ylabel("Molecule")
plt.savefig("families.png", dpi=600)
plt.close()













