"""
AST424 PROJECT - DATA SAMPLING MACRO SCRIPT

Written by Luna.
No artificial intelligence was used in the implementation of this code.

This script fulfils the same purpose as data_ops.py, but generates a package of .pkl files
with varying sample sizes (1000, 2000, 4000, 8000) and a user-selected seed (default 0).

The files cumulatively take up ~410MB, ensure to have enough storage space available
before running this script!
"""

import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
ROOT = os.path.dirname(currentdir)
sys.path.insert(0, ROOT)

import lsdb
import numpy as np
import pickle as pkl

from data_ops import save_data
from misc.filters import *


SPL_SIZES = [1000, 2000, 4000, 8000]


def main(seed):
    """
    Main function.
    """
    if seed != '':
        np.random.seed(int(seed))
    else:
        np.random.seed(0)
        
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

    ia_rands = [np.random.choice(np.arange(0, len(ia_full), 1), spl_size, replace=False) for spl_size in SPL_SIZES]
    iax_rands = [np.random.choice(np.arange(0, len(iax_full), 1), spl_size, replace=False) for spl_size in SPL_SIZES]

    print("Random sampling...")

    for i in range(len(ia_rands)):
        print(f"Saving {SPL_SIZES[i]}-sample data to files...")

        ia_sub = ia_full.iloc[ia_rands[i].tolist(), :]
        save_data(ia_sub, "Ia", seed, SPL_SIZES[i])

        iax_sub = iax_full.iloc[iax_rands[i].tolist(), :]
        save_data(iax_sub, "Iax", seed, SPL_SIZES[i])

        print(f"{SPL_SIZES[i]}-sample data files saved successfully.")

    return 0

if __name__ == "__main__":
    seed = input("Input random selection seed [default = 0]: ")
    main(seed)