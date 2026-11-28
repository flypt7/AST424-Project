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

    chi2 = (ia_chi2, iax_chi2) = (dict(), dict())

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

                flat_result = sncosmo.flatten_result(result)

                if 0 < flat_result["chisq"] < 400:
                    chi2[i][sn] = [flat_result["chisq"], flat_result["c"], flat_result["x0"], flat_result["x1"]]

        i += 1

    params = ["chisq", "c", "x0", "x1"]
    param_ticks = [5, 0.05, 0.001, 0.2]

    fig, axs = plt.subplots(len(params), 1)
    plt.subplots_adjust(hspace=0.4)

    plt.title(f"SALT2 Properties [\chi^2 < 50]")

    for i in range(len(params)):
        lo = np.round(np.min([ia_chi2[sn][i] for sn in ia_chi2]), 1) - 0.5
        hi = np.round(np.max([ia_chi2[sn][i] for sn in ia_chi2]), 1) + 0.5

        axs[i].hist([ia_chi2[sn][i] for sn in ia_chi2], bins=np.arange(lo, hi, param_ticks[i]), color="#785EF0", label=f"SNe-Ia ({len(ia_chi2)} objects)", alpha=0.5)
        axs[i].hist([iax_chi2[sn][i] for sn in iax_chi2], bins=np.arange(lo, hi, param_ticks[i]), color="#FE6100", label=f"SNe-Iax ({len(iax_chi2)} objects)", alpha=0.5)
        
        if i == 0:
            axs[i].legend(loc="upper right", fontsize=16)
        
        axs[i].set_title(params[i].capitalize(), fontsize=22)
        axs[i].set_ylabel("Count", fontsize=18)
        axs[i].set_xlabel(params[i].capitalize(), fontsize=18)
        axs[i].tick_params(labelsize=16)

    plt.show()


if __name__ == "__main__":
    main()