import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from math import *

# import fO2bufferplotter, which is in another directory
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from fO2bufferplotter import buffers

# Parameters the user might wish to change
pressure = 1 # Pressure in bar
temperature_min = 600.0
temperature_max = 1200.0
temperature_step = 1 # Default = 1

# Figure quality
dpi = 250 # Default = 250

# Choose which buffers to plot
# comment out any buffers you don't wan to plot
buffer_list = [
	'HM',
    'ReReO',
    'NNO',	
    'QFM', 
    'CoCoO',	
    'IW',  	
    'Graphite',
    # 'QIF',
    # 'SiSiO2',	
    # 'CrCr2O3',	
    # 'MoMoO2',
    # 'CaCaO',
    # 'AlAl2O3',
    # 'KK2O',
    # 'MgMgO',
    # 'MnMnO',
    # 'NaNa2O',
    # 'TiTiO2',
    ]

# The python call
buffers.plot_log_fO2(
	pressure,
	temperature_min,
	temperature_max,
	temperature_step=temperature_step,
	buffers=buffer_list,
	dpi=dpi,
	)
