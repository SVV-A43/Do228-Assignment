import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_fwf("B737.rpt", skiprows=20, error_bad_lines=False)

print(df)