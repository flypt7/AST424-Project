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
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import lsdb
import numpy as np
import pandas as pd
import pickle

from misc.filters import *


def save_data(data, name):
    """
    Save the data into a file with filename <name>_<N>.pkl
    """
    try:
        with open(f"{name}_{N}.pkl", "wb") as f:
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
     

def main():
    """
    Main function. Returns 0 on success, 1 on failure to save.
    """
    try:
        plasticc = lsdb.open_catalog("hf://datasets/UniverseTBD/mmu_plasticc")
        print("PLAsTiCC Loaded.")

        # Filter by type and redshift
        sne_ia = plasticc[(plasticc["obj_type"] == "SNIa") &
                        (plasticc["redshift"] >= 0) &
                        (plasticc["redshift"] < Z_MAX)].compute()
        print("SNe-Ia Selected.")
        
        sne_iax = plasticc[(plasticc["obj_type"] == "SNIax") &
                        (plasticc["redshift"] >= 0) &
                        (plasticc["redshift"] < Z_MAX)].compute()
        print("SNe-Iax Selected.")

        # MMU_PLAsTiCC contains each object twice when loaded with LSDB, 
        # so we drop duplicates.
        ia_full = sne_ia.drop_duplicates(subset="object_id")
        iax_full = sne_iax.drop_duplicates(subset="object_id")

        print("Duplicates removed.")

        # saving candidates
        save_data_full(ia_full, "Ia")
        print(f"SNe-Ia saved successfully to Ia_z{Z_MAX}_FULL.pkl. Any data in the file was overwritten.")

        save_data_full(iax_full, "Iax")
        print(f"SNe-Ia saved successfully to Iax_z{Z_MAX}_FULL.pkl. Any data in the file was overwritten.")

        return 0

    except:
        print("Something went wrong...")
        return 1

if __name__ == "__main__":
    main()