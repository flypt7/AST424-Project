"""
AST424 PROJECT - PLAsTiCC TO CLEAN LIGHTCURVE DICT

Written by Luna.
No artificial intelligence was used in the implementation of this code.

This file uses the NestedFrame data of SNe-Ia and SNe-Iax lightcurves and
processes the data to return two dictionaries of object_id: (lightcurve, redshift) 
pairs in the BAND specified in misc/filters.py (NONE => all bands are kept).
"""

import os
import sys
import inspect

ROOT = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0, ROOT)

import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit

from data.data_ops import load_data
from misc.filters import *

def clean_lightcurve(lc, band=None):
    """
    Cleaned up lightcurve data from lc in the specified band in-place.
    Cleaning up involves the removal of any non-existent and infinite values,
    and any non-positive errors.
    """
    if band is None:
        lc = lc[(lc["time"] > 0) & np.isfinite(lc["time"]) &
                np.isfinite(lc["flux"]) & np.isfinite(lc["flux_err"]) &
                (lc["flux_err"] > 0) & (np.abs(lc["flux"] / lc["flux_err"]) > 5)]
    else:
        lc = lc[(lc["time"] > 0) & np.isfinite(lc["time"]) &
                np.isfinite(lc["flux"]) & np.isfinite(lc["flux_err"]) &
                (lc["flux_err"] > 0) & (lc["band"] == band) & 
                (np.abs(lc["flux"] / lc["flux_err"]) > 5)]

    # remove lightcurves with no noticeable peak + lightcurves with too few points
    if np.max(lc["flux"]) > 200 and len(lc) > 5:
        return lc.sort_values("time")


def centroid(peak, around, mags):
    """
    Return the centroid location of the peak.
    peak is the location of the peak in the lightcurve data
    around is the number of lightcurve values to take around peak for the 
        centroid calculation
    mags is the lightcurve magnitude data
    """
    val = 0
    value_sum = 0
    for i in range((peak - around), (peak + around + 1)):
        val += mags[i] * i
        value_sum += mags[i]
    
    return val/value_sum


def main(spl_ct):
    """
    Main function. Uses all other functions to turn NestedFrame lightcurve data
    into two dictionaries of object_id: (lightcurve, redshift) pairs in a tuple.
    spl_ct must be a number of samples specified in the filename of at least one
    .pkl file for both SNe-Ia and SNe-Iax.
    Returns the dictionaries on success, None with an error message on failure.
    """
    try:
        # load lightcurve data
        ia = load_data(f"{ROOT}/data/Ia_z0.2_s0_{spl_ct}.pkl")
        print(f"SNe-Ia dataset loaded with {spl_ct} samples.")

        iax = load_data(f"{ROOT}/data/Iax_z0.2_s0_{spl_ct}.pkl")
        print(f"SNe-Iax dataset loaded with {spl_ct} samples.")

        # initialize dicts for SNe-Ia and SNe-Iax (obj_id: lightcurve)
        clean_ia, clean_iax = dict(), dict()

        # remove garbage data, order lightcurve in time, redshift-correct time values
        print("Cleaning up data...")

        for row in ia.iterrows():
            clean = clean_lightcurve(row[1]["lightcurve"], band=None)
            if clean is not None:
                clean_ia[row[0]] = (clean, row[1]["redshift"])
        print(f"SNe-Ia lightcurves cleaned up. {len(clean_ia)} SNe-Ia returned data.")

        for row in iax.iterrows():
            clean = clean_lightcurve(row[1]["lightcurve"], band=None)
            if clean is not None:
                clean_iax[row[0]] = (clean, row[1]["redshift"])
        print(f"SNe-Iax lightcurves cleaned up. {len(clean_iax)} SNe-Iax returned data.")

        return (clean_ia, clean_iax)
    
    except:
        print("An error occurred. The data was likely not loaded properly.")
        return None

# TODO:
#   1. redshift correction                          done
#   2. filter data                                  min number of points before peak?
#   3. peak centroid location                       implemented, to apply
#   4. zero location (Gaussian/polyfit)             curve_fit(norm), value threshold
#                                                           -> do different fits work better for SNIa vs SNIax? check covariance matrix
#   5. time difference between the two              simple subtraction
#   6. boom, send back the full thing               some data storage, done in main function

if __name__ == "__main__":
    main(8000)
