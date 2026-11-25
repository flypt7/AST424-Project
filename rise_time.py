"""
AST424 PROJECT - LIGHTCURVE TO RISE TIME ALGORITHM

Written by Luna.
No artificial intelligence was used in the implementation of this code.

This file uses the dictionaries from PLAsTiCC data created in clean_lc.py
to estimate the rise time of each object.
The "rise time" of a supernova lightcurve is defined as the time between 
the supernova's first changes in flux to its peak flux. Both are determined
using SALT2 fitting through the sncosmo library.

This file also plots a histogram of rise times of SNe-Ia and SNe-Iax rise time
for those with less than 10% fractional uncertainty (this is subject to change
as the code gets developed).
"""

from astropy.table import Table
import numpy as np
import sncosmo

from clean_lc import main as clean_lc 
from misc.filters import *
import matplotlib.pyplot as plt


def fill_table(table, data):
    """
    TODO: write this
    """
    table["time"] = data["time"]
    table["band"] = f"lsst" + data["band"]
    table["flux"] = data["flux"]
    table["fluxerr"] = data["flux_err"]
    table["zp"] = ZERO_POINT # see https://doi.org/10.1051/0004-6361/202245189 
    table["zpsys"] = "ab"    # idem.


def main():
    """
    TODO: write this
    """
    cleans = clean_lc()

    rises = (ia_rise, iax_rise) = (dict(), dict())

    i = 0

    for clean in cleans:
        for sn in clean:
            if not clean[sn][0].empty:
                # prep model
                model = sncosmo.Model(source='salt2')
                model.set(z=clean[sn][1])

                # reshape data to appropriate format for sncosmo
                ia_table = Table()
                fill_table(ia_table, clean[sn][0])

                # run fit
                result, fitted_model = sncosmo.fit_lc(ia_table, model,
                                            ['t0', 'x0', 'x1', 'c'], guess_z=False, minsnr=3)

                peak_phases = []

                for band in ["lsstu", "lsstg", "lsstr", "lssti"]:  # TODO: lsstz throws an error when called here. why?
                    peak_phases.append(fitted_model.source.peakphase(band))

                flat_result = sncosmo.flatten_result(result)

                t0_t0_cov = flat_result["t0_t0_cov"]

                # TODO: might want to split rise times per band instead of aggregating, especially since not all
                #       SNe have the same bands given their data
                peak_phases = np.array(peak_phases)
                peak_avg = np.average(peak_phases)
                peak_err = np.std(peak_phases)

                rise_time = peak_avg - fitted_model.source.minphase()
                rise_time_err = (peak_err ** 2 + t0_t0_cov) ** (1/2)  # simple quadrature uncertainty propagation

                if rise_time / rise_time_err > 10: # arbitrary threshold; 10% fractional uncertainty maximum
                    rises[i][sn] = [rise_time, rise_time_err]

        i += 1

    fig, ax = plt.subplots(1, 1)
    ax.hist([ia_rise[sn][0] for sn in ia_rise], bins=np.arange(18.5, 20, 0.05), color="#1F1E33", label="SNe-Ia", alpha=0.3)
    ax.hist([iax_rise[sn][0] for sn in iax_rise], bins=np.arange(18.5, 20, 0.05), color="#1AEC55", label="SNe-Iax", alpha=0.3)
    ax.legend(loc="upper right", fontsize=16)

    plt.show()


if __name__ == "__main__":
    main()