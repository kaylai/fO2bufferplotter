import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from math import *
from Plot_fO2_buffers import *

#Parameters the user might wish to change
pressure = 1 #Pressure in bar
temperature_min = 600.0
temperature_max = 1200.0
temperature_step = 1 #Default

buffers = ['NNO',
		   'IW',
		   'QFM',
		   'QIF',
		   'HM',
		   'CoCoO',
		   'ReReO',
		   'Graphite',
		   'SiSiO2']

#The python call
plotfO2(pressure, temperature_min, temperature_max, temperature_step=temperature_step, buffers=buffers)