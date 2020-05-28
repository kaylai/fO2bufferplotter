import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from math import *

#Parameters the user might wish to change
P = 1 #Pressure in bar
P_GPa = P/10000

#Define X values (Temp in degrees C)
x = np.arange(950.0, 1400.0, 1)
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
		extrap_upper_range = None

	return usable_temp_range_C, usable_temp_range_K, extrap_lower_range, extrap_upper_range

"""
-------------Buffer values in terms of logfO2-------------
Buffer equations from B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25 unless noted.
"""

""" 
-----------NNO-----------
Define NNO buffer value at P
Equation from Campbell et al. (2009) High-pressure effects on the iron-iron oxide and nickel-nickel oxide oxygen fugacity buffers

T in K, P in GPa
Polynomial coefficients to:  log fO2  =  (a0+a1*P+a2*P^2+a3*P^3+a4*P^4) + (b0+b1*P+b2*P^2+b3*P^3)/T
a0	a1	a2	a3	a4	b0	b1	b2	b3
8.699	0.01642	-0.0002755	0.000002683	-1.015E-08	-24205	444.73	-0.59288	0.0015292
"""
def calc_NNO(P, T):
	fO2 = (8.699 + 0.01642*P -0.0003*P**2 + (2.7*10**(-6))*P**3 - (10**(-8))*P**4) + (-24205 + 444.73*P - 0.5929*P**2 + 0.00153*P**3)/T
	return fO2

temp_range_NNO, temp_range_NNO_K, NNO_extrap_lower, NNO_extrap_upper = get_usable_temp_range((700.0-273.15), (1700.0-273.15))

NNO_y = calc_NNO(P_GPa, temp_range_NNO_K)
if NNO_extrap_lower is not None:
	NNO_extrap_lower_y = calc_NNO(P_GPa, NNO_extrap_lower)
if NNO_extrap_upper is not None:
	NNO_extrap_upper_y = calc_NNO(P_GPa, NNO_extrap_upper)

#Define QFM buffer value at P
QFM_y = (-25096.3/x_K) + 8.735 + 0.11 * (P-1)/x_K


def calc_IW(P, T):
	"""
	Fe-FeO (Iron-Wustite)
	=====================
	Define IW buffer value at P
	
	References
	----------
	Campbell et al. (2009) High-pressure effects on the iron-iron oxide and nickel-nickel oxide oxygen fugacity buffers

	Parameters
	----------
	P float
		Pressure in GPa

	T flaot
		Temperature in degrees K

	Returns
	-------
	float
		fO2

	Polynomial coefficients
	-----------------------
	log fO2  =  (a0+a1*P) + (b0+b1*P+b2*P^2+b3*P^3)/T
	a0: 6.54106	
	a1: 0.0012324
	b0:	-28163.6
	b1:	546.32
	b2:	-1.13412
	b3: 0.0019274
						
	"""
	fO2 = (6.54106+0.0012324*P) + (-28163.6+546.32*P-1.13412*P**2+0.0019274*P**3)/T
	return fO2

temp_range_IW, temp_range_IW_K, IW_extrap_lower, IW_extrap_upper = get_usable_temp_range((1000.0-273.15), (2600.0-273.15))

IW_y = calc_IW(P_GPa, temp_range_IW_K)
if IW_extrap_lower is not None:
	IW_extrap_lower_y = calc_NNO(P_GPa, IW_extrap_lower)
if IW_extrap_upper is not None:
	IW_extrap_upper_y = calc_NNO(P_GPa, IW_extrap_upper)

#Define MH buffer value at P
MH_y = (-25700.6/x_K) + 14.558 + 0.019 * (P-1)/x_K

#Define CoCoO buffer value at P
CoCoO_y = (-24332.6/x_K) + 7.295 + 0.052 * (P-1)/x_K

#Define ReReO buffer value at P
#Equation from Pownceby and O'Neill (1994) Thermodynamic data from redox reactions at high temperatures. IV. Calibration of the Re-ReO2 oxygen buffer from EMF and NiO+Ni-Pd redox sensor measurements
RRO_y = (-451020 + 297.595 * x_K - 14.6585 * x_K * np.log(x_K))/(8.31441 * x_K * log(10))

#Define graphite-gas (graphite-CO-CO2 aka CCO) buffer value at P
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

def calc_SiSiO2(P, T):
	"""
	Si-SiO2
	=======
	Define the silicon-silicon dioxide buffer value at P
	Equation from Kathleen Vander Kaaden (pers. comm), with thermodynamic data taken from the JANAF tables in Chase (1998)

	Parameters
	----------
	P float
		Pressure in GPa

	T flaot
		Temperature in degrees K

	Returns
	-------
	float
		fO2

	References
	----------
	Chase, M. W. (1998). NIST-JANAF thermochemical tables. In Journal of Physical and Chemical Reference Data, 9 (1961 pp.). Gaithersburg, MD:
	National Institute of Standards and Technology.

	Barin 1993, Phases of Silicon at High Pressure

	Hu 1984, Thermochemical Data of Pure Substances (v.1 and V.2), CSEL QD 511.8 B369 1993

	Fried, Howard, and Souers. EXP6: A new EOS library for HP Thermochemistry

	Murnaghan parameters
	--------------------
	Vo (ML/mol) - Fe:7.11, FeO:12.6, Si:12.06, SiO2:22.68
	Bo (GPa) - Fe:139.0, FeO:142.5, Si:97.9, SiO2:27.0
	Bo' - Fe:4.7, FeO:5.0, Si:4.16, SiO2:3.8
	"""




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
	NNO_extrap_upper_plt, = ax1.plot(NNO_extrap_upper, NNO_extrap_upper_y, color='green', linestyle='dashed', linewidth=1.5)
	figure_labels.append('NNO Extrapolated Upper')

#QFM
QFM_plt, = ax1.plot(x, QFM_y, color='purple',linewidth=1.5)
figure_labels.append('QFM Buffer')

#IW
IW_plt, = ax1.plot(x, IW_y, color='blue', linewidth=1.5)
figure_labels.append('IW Buffer')
if IW_extrap_lower is not None:
	IW_extrap_lower_plt, = ax1.plot(IW_extrap_lower, IW_extrap_lower_y, color='blue', linestyle='dashed', linewidth=1.5)
	figure_labels.append('IW Extrapolated Lower')
if IW_extrap_upper is not None:
	IW_extrap_upper_plt, = ax1.plot(IW_extrap_upper, IW_extrap_upper_y, color='blue', linestyle='dashed', linewidth=1.5)
	figure_labels.append('IW Extrapolated Upper')

#MH
# MH_plt, = ax1.plot(x, MH_y, color='red', linewidth=1.5, label='MH Buffer')
# figure_labels.append('MH Buffer')

#CoCoO
# CoCoO_plt, = ax1.plot(x, CoCoO_y, color='magenta', linewidth=1.5, label='CoCoO Buffer')
# figure_labels.append('CCO Buffer')

#RRO
# RRO_plt, = ax1.plot(x, RRO_y, color='orange', linewidth=1.5, label='RRO Buffer')
# figure_labels.append('RRO Buffer')

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