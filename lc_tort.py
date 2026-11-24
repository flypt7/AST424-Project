"""
AST424 PROJECT - LIGHTCURVE TO RISE TIME ALGORITHM

Written by Luna.
No artificial intelligence was used in the implementation of this code.

This file uses the NestedFrame data of SNe-Ia and SNe-Iax lightcurves and
processes the data to return a list of (obj_id, object_type, rt) 3-ples,
where rt is a dictionary of rise times per band. Specific band data can then
be selected by dictionary lookup in the last item of the 3-ple.
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


def clean_lightcurve(lc, z, band=None):
    """
    Cleaned up lightcurve data from lc in the specified band in-place.
    Cleaning up involves the removal of any non-existent and infinite values,
    and any non-positive errors. Also accounts for redshift z.
    """
    if band is None:
        lc = lc[(lc["time"] > 0) & np.isfinite(lc["time"]) &
                np.isfinite(lc["flux"]) & np.isfinite(lc["flux_err"]) &
                (lc["flux_err"] > 0)]
    else:
        lc = lc[(lc["time"] > 0) & np.isfinite(lc["time"]) &
                np.isfinite(lc["flux"]) & np.isfinite(lc["flux_err"]) &
                (lc["flux_err"] > 0) & (lc["band"] == band)]

    # account for redshift
    time_ratio = 1 / (1 + z)
    lc["time"] = lc["time"].apply(lambda x: x * time_ratio)

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


def main():
    """
    Main function. Uses all other functions to turn NestedFrame lightcurve data
    into the (obj_id, object_type, rt) 3-ple.
    Returns the 3-ple on success, None with an error message on failure.
    """
    # load lightcurve data
    ia = load_data(f"{ROOT}/data/Ia_z0.2_s0_8000.pkl")
    print("SNe-Ia dataset loaded with 8000 samples.")

    iax = load_data(f"{ROOT}/data/Iax_z0.2_s0_8000.pkl")
    print("SNe-Iax dataset loaded with 8000 samples.")

    # initialize dicts for SNe-Ia and SNe-Iax (obj_id: lightcurve)
    clean_ia, clean_iax = dict(), dict()

    # remove garbage data, order lightcurve in time, redshift-correct time values
    print("Cleaning up data...")

    for row in ia.iterrows():
        clean_ia[row[0]] = clean_lightcurve(row[1]["lightcurve"], row[1]["redshift"], band=BAND)
    print(f"SNe-Ia lightcurves cleaned up. {len(clean_ia)} SNe-Ia returned data.")

    for row in iax.iterrows():
        clean_iax[row[0]] = clean_lightcurve(row[1]["lightcurve"], row[1]["redshift"], band=BAND)
    print(f"SNe-Iax lightcurves cleaned up. {len(clean_iax)} SNe-Iax returned data.")


# TODO:
#   1. redshift correction                          done
#   2. filter data                                  min number of points before peak?
#   3. peak centroid location                       implemented, to apply
#   4. zero location (Gaussian/polyfit)             curve_fit(norm), value threshold
#                                                           -> do different fits work better for SNIa vs SNIax? check covariance matrix
#   5. time difference between the two              simple subtraction
#   6. boom, send back the full thing               some data storage, done in main function

if __name__ == "__main__":
    main()