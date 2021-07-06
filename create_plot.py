import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from math import *
from fO2bufferplotter import buffers

# Parameters the user might wish to change
pressure = 1 # Pressure in bar
temperature_min = 600.0
temperature_max = 1200.0
temperature_step = 1 # Default

# Choose which buffers to plot
# comment out any buffers you don't wan to plot
buffer_list = ['NNO',	
		       'QFM', 	
		       'IW',  	
		       'HM', 	
		       'CoCoO', 	
		       'ReReO',
		       'Graphite',
		       'QIF',
		       'SiSiO2'	
		       'CrCr2O3',	
		       'MoMoO2',
		       'CaCaO',
		       'AlAl2O3',
		       'KK2O',
		       'MgMgO',
		       'MnMnO',
		       'NaNa2O',
		       'TiTiO2'
		       ]

# The python call
buffers.plot_log_fO2(pressure, temperature_min, temperature_max, temperature_step=temperature_step, buffers=buffer_list)