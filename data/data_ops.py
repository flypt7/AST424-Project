"""
AST424 PROJECT - SAVE DATA

Written by Luna.
No artificial intelligence was used in the implementation of this code.

This file loads the SNe-Ia and SNe-Iax rows from MMU_PLAsTiCC and saves them as binary 
files using Pickle. Data files are saved in the <type>_<sample_size>.pkl format in the
same directory as this file.
"""
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
ROOT = os.path.dirname(currentdir)
sys.path.insert(0, ROOT)

import lsdb
import numpy as np
import pickle

from misc.filters import *


def save_data(data, name, seed):
    """
    Save the data into a file with filename <name>_<N>.pkl.
    """
    try:
        with open(f"{ROOT}/data/{name}_z{Z_MAX}_s{seed}_{N}.pkl", "wb") as f:
            pickle.dump(data, f)
    
    except:
        print(f"I/O calls failed for {name}.")
        
def save_data_full(data, name):
    """
    Save the data into a file with filename <name>_<details>_FULL.pkl
    """
    try:
        with open(f"{name}_z{Z_MAX}_FULL.pkl", "wb") as f:
            pickle.dump(data, f)
    
    except:
        print(f"I/O calls failed for {name} full saving.")

def load_data(path):
    """
    Return the object stored in the file with the corresponding path.
    If the file does not exist, return None.
    """
    try:
        with open(path, "rb") as f:
            return pickle.load(f)
        
    except IOError:
        print(f"Error: No such file or directory ({path}). No contents were loaded.")
     

def main(spl, seed):
    """
    Main function. Write the full MMU_PLAsTiCC SNe-Ia and SNe-Iax queries (with redshift requirements 
    set in misc/filters.py) if spl is None, or the specified sample size (spl) randomly selected from 
    the full queries, to separate .pkl files in the same directory as this file.
    The seed parameter is used for the random sample selection in the case where spl is not None. The
    default seed is 0.
    The full SNe-Ia file contains ~78000 objects (~1.1GB) and the SNe-Iax ~9000 (~110MB). The queries 
    require approximately 4.2GB of memory, and this function takes roughly 10-15 minutes to return.
    Return 0 on success.
    """
    np.random.seed(seed)

    print(f"Seed: {seed}\t Sample size: {spl}")
    
    print("Loading PLAsTiCC...")
    plasticc = lsdb.open_catalog("hf://datasets/UniverseTBD/mmu_plasticc")
    print("PLAsTiCC Loaded.")

    # Filter by type and redshift
    print("Querying SNe-Ia...")
    sne_ia = plasticc[(plasticc["obj_type"] == "SNIa") &
                    (plasticc["redshift"] >= 0) &
                    (plasticc["redshift"] < Z_MAX)].compute()
    print("SNe-Ia Selected.")

    print("Querying SNe-Iax...")
    sne_iax = plasticc[(plasticc["obj_type"] == "SNIax") &
                    (plasticc["redshift"] >= 0) &
                    (plasticc["redshift"] < Z_MAX)].compute()
    print("SNe-Iax Selected.")
    
    # MMU_PLAsTiCC contains each object twice when loaded with LSDB, 
    # so we drop duplicates.
    ia_full = sne_ia.drop_duplicates(subset="object_id")
    iax_full = sne_iax.drop_duplicates(subset="object_id")

    print("Duplicates removed from both selections.")

    # saving candidates
    if spl is None:
        print("Saving SNe-Ia data to file...")
        save_data_full(ia_full, "Ia")
        print(f"SNe-Ia saved successfully to Ia_z{Z_MAX}_FULL.pkl. Any data in the file was overwritten.")

        print("Saving SNe-Iax data to file...")
        save_data_full(iax_full, "Iax")
        print(f"SNe-Ia saved successfully to Iax_z{Z_MAX}_FULL.pkl. Any data in the file was overwritten.")

    else:

        ia_rand = np.random.choice(np.arange(0, len(ia_full), 1), spl, replace=False)
        iax_rand = np.random.choice(np.arange(0, len(iax_full), 1), spl, replace=False)

        print("Random sampling...")
        ia_sub = ia_full.iloc[ia_rand.tolist(), :]
        iax_sub = iax_full.iloc[iax_rand.tolist(), :]

        print("Saving SNe-Ia data to file...")
        save_data(ia_sub, "Ia", seed)
        print(f"SNe-Ia saved successfully to Ia_z{Z_MAX}_s{seed}_{N}.pkl. Any data in the file was overwritten.")

        print("Saving SNe-Iax data to file...")
        save_data(iax_sub, "Iax", seed)
        print(f"SNe-Iax saved successfully to Iax_z{Z_MAX}_s{seed}_{N}.pkl. Any data in the file was overwritten.")

    return 0

if __name__ == "__main__":
    seed = 0    # default seed for random sample selection
    
    spl = input(f"Number of samples to save (int) [no input = all, 'cst' = {N}]: ")
    try:
        if spl == '':
            spl = None
        elif spl == "cst":
            spl = N
            int(spl)

    except:
        print(f"Incorrect value {spl}. Input must be positive non-zero integer.")

    seed = input(f"Seed for random selection (int) [no input = 0]: ")
    if seed == '':
        seed = 0
    else:    
        int(seed)

    main(spl, seed)
