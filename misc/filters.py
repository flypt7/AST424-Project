"""
AST424 PROJECT - FILTERS FILE

This file contains all of the fixed requirements for supernova selection. 
These are constant between detection methods, meaning that this file is imported
by all feature detection notebooks.
"""


##########
# GENERAL FILTERS
##########

N = 8000        # Sample count
Z_MAX = 0.2
BAND = "z"

##########
# SECOND PEAK DETECTION
##########

MIN_PRIMARY_SNR = 5
MIN_SECONDARY_SNR = 3
MIN_RISE_SIGMA = 1
MIN_FALL_SIGMA = 1

MAX_GAP = 10
MIN_DAYS = 10
MAX_DAYS = 40

# How far back to search for a preceding minimum
MINIMUM_WINDOW = 15