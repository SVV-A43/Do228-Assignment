import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

bending_region1 = pd.read_fwf("B737.rpt", header=None, skiprows=20, nrows=5778, error_bad_lines=False)
bending_region2 = pd.read_fwf("B737.rpt", header=None, skiprows=5816, nrows=856, error_bad_lines=False)

jam_bent_region1 = pd.read_fwf("B737.rpt", header=None, skiprows=6705, nrows=5778, error_bad_lines=False)
jam_bent_region2 = pd.read_fwf("B737.rpt", header=None, skiprows=12501, nrows=856, error_bad_lines=False)

jam_straight_region1 = pd.read_fwf("B737.rpt", header=None, skiprows=13390, nrows=5778, error_bad_lines=False)
jam_straight_region2 = pd.read_fwf("B737.rpt", header=None, skiprows=19186, nrows=856, error_bad_lines=False)

