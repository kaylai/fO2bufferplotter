import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from math import *


#Parameters the user might wish to change
P = 1 #Pressure in bar

#Define X values (Temp in degrees C)
x = np.arange(50.0, 1400.0, 1)
x_K = x + 273.15 #Temp in K

#Define temp ranges over which buffers can be plotted
def get_usable_temp_range(minT, maxT):
	extrap_lower_min = None
	extrap_upper_min = None
	if x[0] >= minT:
		min_range = x[0]
	else:
		min_range = minT
		extrap_lower_min = x[0]
		extrap_lower_max = minT

	if x[-1]+1 <= maxT:
		max_range = x[-1]+1
	else:
		max_range = maxT
		extrap_upper_min = maxT
		extrap_upper_max = x[-1] + 1

	usable_temp_range_C = np.arange(min_range, max_range, 1)
	usable_temp_range_K = usable_temp_range_C + 273.15

	if extrap_lower_min is not None:
		extrap_lower_range = np.arange(extrap_lower_min, extrap_lower_max, 1)
	else:
		extrap_lower_range = None

	if extrap_upper_min is not None:
		extrap_upper_range = np.arange(extrap_upper_min, extrap_upper_max, 1)
	else:
		extrap_upper_rane = None

	return usable_temp_range_C, usable_temp_range_K, extrap_lower_range, extrap_upper_range

#Buffer values in terms of logfO2
#Buffer equations from B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25 unless noted.

#Define NNO buffer value at P
#Equation from O'Neill and Pownceby (1993) Thermodynamic data from redox reactions at high temperatures. I.
temp_range_NNO, temp_range_NNO_K, NNO_extrap_lower, NNO_extrap_upper = get_usable_temp_range(700.0-273.15, 1700.0-273.15)
NNO_y = -478967 + 248.514*temp_range_NNO_K - 12.9053*temp_range_NNO_K*np.log(temp_range_NNO_K)
if NNO_extrap_lower is not None:
	NNO_extrap_lower_y = -478967 + 248.514*(NNO_extrap_lower+273.15) - 12.9053*(NNO_extrap_lower+273.15)*np.log(NNO_extrap_lower+273.15)
if NNO_extrap_upper is not None:
	NNO_extrap_lower_y = -478967 + 248.514*(NNO_extrap_upper+273.15) - 12.9053*(NNO_extrap_upper+273.15)*np.log(NNO_extrap_upper+273.15)

#Define QFM buffer value at P
QFM_y = (-25096.3/x_K) + 8.735 + 0.11 * (P-1)/x_K

#Define IW buffer value at P
IW_y = (-27489/x_K) + 6.702 + 0.055 * (P-1)/x_K

#Define MH buffer value at P
MH_y = (-25700.6/x_K) + 14.558 + 0.019 * (P-1)/x_K

#Define CoCoO buffer value at P
CCO_y = (-24332.6/x_K) + 7.295 + 0.052 * (P-1)/x_K

#Define ReReO buffer value at P
#Equation from Pownceby and O'Neill (1994) Thermodynamic data from redox reactions at high temperatures. IV. Calibration of the Re-ReO2 oxygen buffer from EMF and NiO+Ni-Pd redox sensor measurements
RRO_y = (-451020 + 297.595 * x_K - 14.6585 * x_K * np.log(x_K))/(8.31441 * x_K * log(10))

#Define graphite-gas (graphite-CO-CO2) buffer value at P
#Equation from French and Eugster (1965) Journal of Geophysical Research
Graphite_y = (-20586/x_K) - 0.044 + np.log10(P) - 0.028 * (P-1)/x_K

#Define QIF low temperature buffer value at P
temp_range_QIF_lowT, temp_range_QIF_lowT_K, QIF_lowT_extrap_lower, QIF_lowT_extrap_upper = get_usable_temp_range(150.0, 573.0)

QIF_lowT_y = (-29435.7/temp_range_QIF_lowT_K) + 7.391 + 0.44 * (P-1)/temp_range_QIF_lowT_K
if QIF_lowT_extrap_lower is not None:
	QIF_lowT_extrap_lower_y = (-29435.7/(QIF_lowT_extrap_lower+273.15)) + 7.391 + 0.44 * (P-1)/(QIF_lowT_extrap_lower+273.15)

#Define QIF high temperature buffer value at P
temp_range_QIF_highT, temp_range_QIF_highT_K, QIF_highT_extrap_lower, QIF_highT_extrap_upper = get_usable_temp_range(574.0, 1200.0)
if QIF_highT_extrap_upper is not None:
	QIF_highT_extrap_upper_y = (-29520.8/(QIF_highT_extrap_upper+273.15)) + 7.492 + 0.050 * (P-1)/(QIF_highT_extrap_upper+273.15)

QIF_highT_y = (-29520.8/temp_range_QIF_highT_K) + 7.492 + 0.050 * (P-1)/temp_range_QIF_highT_K




#DRAW THE FIGURE
fig, ax1 = plt.subplots()

figure_labels = []

#The comma here must be here for the legend to work. I HAVE NO IDEA WHY.
#NNO
NNO_plt, = ax1.plot(temp_range_NNO, NNO_y, color='green', linewidth=1.5)
figure_labels.append('NNO Buffer')
if NNO_extrap_lower is not None:
	NNO_extrap_lower_plt, = ax1.plot(NNO_extrap_lower, NNO_extrap_lower_y, color='green', linestyle='dashed', linewidth=1.5)
	figure_labels.append('NNO Extrapolated Lower')
if NNO_extrap_upper is not None:
	NNO_extrap_upper_plt, = ax1.plot(NNO_extrap_upper, NNO_extrap_upper_y, color='green', linestyle='dashed', linewidth=1.5, label='NNO Extrapolated')
	figure_labels.append('NNO Extrapolated Upper')

#QFM
QFM_plt, = ax1.plot(x, QFM_y, color='purple',linewidth=1.5)
figure_labels.append('QFM Buffer')

#IW
IW_plt, = ax1.plot(x, IW_y, color='blue', linewidth=1.5)
figure_labels.append('IW Buffer')

#MH
MH_plt, = ax1.plot(x, MH_y, color='red', linewidth=1.5, label='MH Buffer')
figure_labels.append('MH Buffer')

#CCO
CCO_plt, = ax1.plot(x, CCO_y, color='magenta', linewidth=1.5, label='CCO Buffer')
figure_labels.append('CCO Buffer')

#RRO
RRO_plt, = ax1.plot(x, RRO_y, color='orange', linewidth=1.5, label='RRO Buffer')
figure_labels.append('RRO Buffer')

#Graphite
Graphite_plt, = ax1.plot(x, Graphite_y, color='gray', linewidth=1.5, label='Graphite')
figure_labels.append('Graphite Buffer')

#QIF
QIF_lowT_plt, = ax1.plot(temp_range_QIF_lowT, QIF_lowT_y, color='mediumaquamarine', linewidth=1.5)
figure_labels.append('QIF Low T')
if QIF_lowT_extrap_lower is not None:
	QIF_lowT_extrap_lower_plt, = ax1.plot(QIF_lowT_extrap_lower, QIF_lowT_extrap_lower_y, color='mediumaquamarine', linestyle='dashed', linewidth=1.5)
	figure_labels.append('QIF Low T Extrapolated')
QIF_highT_plt, = ax1.plot(temp_range_QIF_highT, QIF_highT_y, color='brown', linewidth=1.5)
figure_labels.append('QIF High T')
if QIF_highT_extrap_upper is not None:
	QIF_highT_extrap_upper_plt, = ax1.plot(QIF_highT_extrap_upper, QIF_highT_extrap_upper_y, color='brown', linestyle='dashed', linewidth=1.5)
	figure_labels.append('QIF High T Extrapolated')


#Define legend styles
labels = [i for i in figure_labels]

#Label axes
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('log $f$O$_2$')

ax1.legend(labels, loc='lower right')
plt.title('Redox buffers at ' + str(P) + ' bar')
plt.show()