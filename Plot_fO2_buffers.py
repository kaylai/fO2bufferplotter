import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from math import *

#-------------DEFINE SOME UNIVERSAL FUNCTIONS------------#
def set_calibration_temp_range(temperature_min, temperature_max, buffer_min_T, buffer_max_T):
	"""
	Set calibrated temperature range for plotting
	=============================================
	Returns variables representing the temperature range over which a particular buffer is calibrated and
	the user-defined temperatures that fall outside of this range. Outside of the calibrated temperature
	range, buffer curves are plotted as dashed lines (inside, as solid lines).

	Parameters
	----------
	buffer_min_T: float
		Minimum calibrated temperature of the buffer in degrees C.

	buffer_max_T: float
		Maximum calibrated temperature of the buffer in degrees C.

	Returns
	-------
	Four numpy arrays
		1. Usable temperature range in degrees C, with steps of 1 degree
		2. Usable temperature range in degrees K, with steps of 1 degree
		3. Any user input temperature values that fall below the calibrated temperature range for this buffer, in K
		4. Any user input temperature values that fall above the calibrated temperature range for this buffer, in K
	"""
	temp_range = np.arange(temperature_min, temperature_max, 1)
	extrap_lower_min = None
	extrap_upper_min = None
	if temp_range[0] >= buffer_min_T:
		min_range = temp_range[0]
	else:
		min_range = buffer_min_T
		extrap_lower_min = temp_range[0]
		extrap_lower_max = buffer_min_T

	if temp_range[-1]+1 <= buffer_max_T:
		max_range = temp_range[-1]+1
	else:
		max_range = buffer_max_T
		extrap_upper_min = buffer_max_T
		extrap_upper_max = temp_range[-1] + 1

	usable_temp_range_C = np.arange(min_range, max_range, 1)
	usable_temp_range_K = usable_temp_range_C + 273.15

	if extrap_lower_min is not None:
		extrap_lower_range = np.arange(extrap_lower_min, extrap_lower_max, 1) + 273.15
	else:
		extrap_lower_range = None

	if extrap_upper_min is not None:
		extrap_upper_range = np.arange(extrap_upper_min, extrap_upper_max, 1) + 273.15
	else:
		extrap_upper_range = None

	return usable_temp_range_C, usable_temp_range_K, extrap_lower_range, extrap_upper_range

#-----------------DEFINE BUFFER EQUATIONS-------------#
def calc_NNO(P, T):
	""" 
	Ni-NiO (NNO)
	============
	Define NNO buffer value at P

	References
	----------
	Campbell et al. (2009) High-pressure effects on the iron-iron oxide and nickel-nickel oxide oxygen fugacity buffers

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	Polynomial coefficients
	-----------------------  
	log fO2  =  (a0+a1*P+a2*P^2+a3*P^3+a4*P^4) + (b0+b1*P+b2*P^2+b3*P^3)/T
	a0: 8.699
	a1: 0.01642
	a2: -0.0002755
	a3: 0.000002683
	a4: -1.015E-08
	b0: -24205
	b1: 444.73
	b2: -0.59288
	b3:0.0015292                            
	"""
	fO2 = (8.699 + 0.01642*P -0.0003*P**2 + (2.7*10**(-6))*P**3 - (10**(-8))*P**4) + (-24205 + 444.73*P - 0.5929*P**2 + 0.00153*P**3)/T
	return fO2

def calc_QFM(P, T):
	"""
	Quartz-Fayalite-Magnetite (QFM)
	===============================
	Define QFM buffer value at P

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	if isinstance(T, float) or isinstance(T, int):
		if T<573:
			fO2 = (-26445.3/T) + 10.344 + 0.092 * (P-1)/T
		if T>=573:
			fO2 = (-25096.3/T) + 8.735 + 0.11 * (P-1)/T

	if isinstance(T, np.ndarray):
		fO2_list = []
		for temp in T:
			if temp<573:
				fO2_list.append((-26445.3/temp) + 10.344 + 0.092 * (P-1)/temp)
			if temp>=573:
				fO2_list.append((-25096.3/temp) + 8.735 + 0.11 * (P-1)/temp)
		fO2 = np.array(fO2_list)

	return fO2

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
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	Polynomial coefficients
	-----------------------
	log fO2  =  (a0+a1*P) + (b0+b1*P+b2*P^2+b3*P^3)/T
	a0: 6.54106 
	a1: 0.0012324
	b0: -28163.6
	b1: 546.32
	b2: -1.13412
	b3: 0.0019274               
	"""
	fO2 = (6.54106+0.0012324*P) + (-28163.6+546.32*P-1.13412*P**2+0.0019274*P**3)/T
	return fO2

def calc_SiSiO2(P, T):
	"""
	Si-SiO2
	=======
	Define the silicon-silicon dioxide buffer value at P
	Equation from Kathleen Vander Kaaden (pers. comm), with thermodynamic data taken from the JANAF tables in Chase (1998)

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
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
	if isinstance(T, int) or isinstance(T, float):
		fO2_1bar = log10(exp((-911336 + 212.2423*T - 4.512*T*log(T))/(8.314*T)))
		fO2 = (1/T*0.4343*120*((0.00251*P**3 - 0.2044*P**2 + 10.32257*P)-(0.00251*0.0001**3 - 0.2004*0.0001**2 + 10.32257*0.001)))+fO2_1bar
	if isinstance(T, np.ndarray):
		fO2_list = []
		for temp in T:
			fO2_1bar = log10(exp((-911336 + 212.2423*temp - 4.512*temp*log(temp))/(8.314*temp)))
			fO2_list.append((1/temp*0.4343*120*((0.00251*P**3 - 0.2044*P**2 + 10.32257*P)-(0.00251*0.0001**3 - 0.2004*0.0001**2 + 10.32257*0.001)))+fO2_1bar)
		fO2 = np.array(fO2_list)
	return fO2

def calc_HM(P, T):
	"""
	Hematite-Magnetite (HM)
	=======================
	Define HM buffer value at P

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	fO2 = (-25700.6/T) + 14.558 + 0.019 * (P-1)/T
	return fO2

def calc_CoCoO(P, T):
	"""
	Cobalt-Cobalt Oxide (CoCoO)
	===========================
	Define CoCoO buffer value at P

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	fO2 = (-24332.6/T) + 7.295 + 0.052 * (P-1)/T
	return fO2

def calc_ReReO(P, T):
	"""
	Rhenium-Rhenium Oxide (ReReO)
	=============================
	Define ReReO buffer value at P

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	References
	----------
	Pownceby and O'Neill (1994) Thermodynamic data from redox reactions at high temperatures. 
	IV. Calibration of the Re-ReO2 oxygen buffer from EMF and NiO+Ni-Pd redox sensor measurements
	"""
	fO2 = (-451020 + 297.595 * T - 14.6585 * T * np.log(T))/(8.31441 * T * log(10))
	return fO2

def calc_Graphite(P, T):
	"""
	Graphite-CO-CO2 gas aka CCO
	===========================
	Define Graphite buffer value at P

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	References
	----------
	French and Eugster (1965) Journal of Geophysical Research
	"""
	fO2 = (-20586/T) - 0.044 + np.log10(P) - 0.028 * (P-1)/T
	return fO2

def calc_QIF(P, T):
	"""
	Quartz-Iron-Fayalite (QIF)
	==========================
	Define QIF buffer value at P

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		fO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	if isinstance(T, float) or isinstance(T, int):
		if T<=573.0:
			fO2 = (-29435.7/T) + 7.391 + 0.44 * (P-1)/T
		if T>573.0:
			fO2 = (-29520.8/T) + 7.492 + 0.050 * (P-1)/T
	if isinstance(T, np.ndarray):
		fO2_list = []
		for temp in T:
			if temp<573:
				fO2_list.append((-29435.7/temp) + 7.391 + 0.44 * (P-1)/temp)
			if temp>=573:
				fO2_list.append((-29520.8/temp) + 7.492 + 0.050 * (P-1)/temp)
		fO2 = np.array(fO2_list)
	return fO2

#--------------PLOTTING-------------#
def plotfO2(pressure, temperature_min, temperature_max, temperature_step=1, buffers=['NNO', 'QFM']):
	"""
	Returns a matplotlib plot of buffer curves at specified P and range of T's.

	Parameters
	----------
	pressure: float
		pressure in bars

	temperature_min: float
		minimum temperature in degrees C

	temperature_max: float
		maximum temperature in degrees C

	temperature_step: int
		OPTIONAL. Default is 1. Step size between temperature values.

	buffers: list
		OPTIONAL. Default is ['NNO', 'QFM']. List with strings of all buffers you wish to plot.
		Possible buffers are: NNO, QFM, IW, SiSiO2, HM, CoCoO, ReReO, Graphite, QIF.
	"""
	PGPa = pressure/10000
	temp_range = np.arange(temperature_min, temperature_max, temperature_step)
	temp_range_K = temp_range + 273.15

	#DRAW THE FIGURE
	fig, ax1 = plt.subplots()
	figure_labels = []

	#for temp in temp_range:
	#NNO
	if 'NNO' in buffers:
		temp_range_NNO, temp_range_NNO_K, NNO_extrap_lower, NNO_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 600, 1200)
		if NNO_extrap_lower is not None:
			NNO_extrap_lower_C = NNO_extrap_lower - 273.15
		if NNO_extrap_upper is not None:
			NNO_extrap_upper_C = NNO_extrap_upper - 273.15

		NNO_y = calc_NNO(PGPa, temp_range_NNO_K)
		if NNO_extrap_lower is not None:
			NNO_extrap_lower_y = calc_NNO(PGPa, NNO_extrap_lower)
		if NNO_extrap_upper is not None:
			NNO_extrap_upper_y = calc_NNO(PGPa, NNO_extrap_upper)

		NNO_plt, = ax1.plot(temp_range_NNO, NNO_y, color='green', linewidth=1.5)
		figure_labels.append('NNO Buffer')
		if NNO_extrap_lower is not None:
			NNO_extrap_lower_plt, = ax1.plot(NNO_extrap_lower_C, NNO_extrap_lower_y, color='green', linestyle='dashed', linewidth=1.5)
			figure_labels.append('NNO Extrapolated Lower')
		if NNO_extrap_upper is not None:
			NNO_extrap_upper_plt, = ax1.plot(NNO_extrap_upper_C, NNO_extrap_upper_y, color='green', linestyle='dashed', linewidth=1.5)
			figure_labels.append('NNO Extrapolated Upper')

	#QFM
	if 'QFM' in buffers:
		temp_range_QFM, temp_range_QFM_K, QFM_extrap_lower, QFM_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 400, 1200)
		if QFM_extrap_lower is not None:
			QFM_extrap_lower_C = QFM_extrap_lower - 273.15
		if QFM_extrap_upper is not None:
			QFM_extrap_upper_C = QFM_extrap_upper - 273.15

		QFM_y = calc_QFM(PGPa, temp_range_QFM_K)
		if QFM_extrap_lower is not None:
			QFM_extrap_lower_y = calc_QFM(PGPa, QFM_extrap_lower)
		if QFM_extrap_upper is not None:
			QFM_extrap_upper_y = calc_QFM(PGPa, QFM_extrap_upper)

		QFM_plt, = ax1.plot(temp_range_QFM, QFM_y, color='darkorange', linewidth=1.5)
		figure_labels.append('QFM Buffer')
		if QFM_extrap_lower is not None:
			QFM_extrap_lower_plt, = ax1.plot(QFM_extrap_lower_C, QFM_extrap_lower_y, color='darkorange', linestyle='dashed', linewidth=1.5)
			figure_labels.append('QFM Extrapolated Lower')
		if QFM_extrap_upper is not None:
			QFM_extrap_upper_plt, = ax1.plot(QFM_extrap_upper_C, QFM_extrap_upper_y, color='darkorange', linestyle='dashed', linewidth=1.5)
			figure_labels.append('QFM Extrapolated Upper')

	#IW
	if 'IW' in buffers:
		temp_range_IW, temp_range_IW_K, IW_extrap_lower, IW_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 565, 1200)
		if IW_extrap_lower is not None:
			IW_extrap_lower_C = IW_extrap_lower - 273.15
		if IW_extrap_upper is not None:
			IW_extrap_upper_C = IW_extrap_upper - 273.15

		IW_y = calc_IW(PGPa, temp_range_IW_K)
		if IW_extrap_lower is not None:
			IW_extrap_lower_y = calc_IW(PGPa, IW_extrap_lower)
		if IW_extrap_upper is not None:
			IW_extrap_upper_y = calc_IW(PGPa, IW_extrap_upper)

		IW_plt, = ax1.plot(temp_range_IW, IW_y, color='black', linewidth=1.5)
		figure_labels.append('IW Buffer')
		if IW_extrap_lower is not None:
			IW_extrap_lower_plt, = ax1.plot(IW_extrap_lower_C, IW_extrap_lower_y, color='black', linestyle='dashed', linewidth=1.5)
			figure_labels.append('IW Extrapolated Lower')
		if IW_extrap_upper is not None:
			IW_extrap_upper_plt, = ax1.plot(IW_extrap_upper_C, IW_extrap_upper_y, color='black', linestyle='dashed', linewidth=1.5)
			figure_labels.append('IW Extrapolated Upper')

	#HM
	if 'HM' in buffers:
		temp_range_HM, temp_range_HM_K, HM_extrap_lower, HM_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 565, 1200)
		if HM_extrap_lower is not None:
			HM_extrap_lower_C = HM_extrap_lower - 273.15
		if HM_extrap_upper is not None:
			HM_extrap_upper_C = HM_extrap_upper - 273.15

		HM_y = calc_HM(PGPa, temp_range_HM_K)
		if HM_extrap_lower is not None:
			HM_extrap_lower_y = calc_HM(PGPa, HM_extrap_lower)
		if HM_extrap_upper is not None:
			HM_extrap_upper_y = calc_HM(PGPa, HM_extrap_upper)

		HM_plt, = ax1.plot(temp_range_HM, HM_y, color='red', linewidth=1.5)
		figure_labels.append('HM Buffer')
		if HM_extrap_lower is not None:
			HM_extrap_lower_plt, = ax1.plot(HM_extrap_lower_C, HM_extrap_lower_y, color='red', linestyle='dashed', linewidth=1.5)
			figure_labels.append('HM Extrapolated Lower')
		if HM_extrap_upper is not None:
			HM_extrap_upper_plt, = ax1.plot(HM_extrap_upper_C, HM_extrap_upper_y, color='red', linestyle='dashed', linewidth=1.5)
			figure_labels.append('HM Extrapolated Upper')

	#CoCoO
	if 'CoCoO' in buffers:
		temp_range_CoCoO, temp_range_CoCoO_K, CoCoO_extrap_lower, CoCoO_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 565, 1200)
		if CoCoO_extrap_lower is not None:
			CoCoO_extrap_lower_C = CoCoO_extrap_lower - 273.15
		if CoCoO_extrap_upper is not None:
			CoCoO_extrap_upper_C = CoCoO_extrap_upper - 273.15

		CoCoO_y = calc_CoCoO(PGPa, temp_range_CoCoO_K)
		if CoCoO_extrap_lower is not None:
			CoCoO_extrap_lower_y = calc_CoCoO(PGPa, CoCoO_extrap_lower)
		if CoCoO_extrap_upper is not None:
			CoCoO_extrap_upper_y = calc_CoCoO(PGPa, CoCoO_extrap_upper)

		CoCoO_plt, = ax1.plot(temp_range_CoCoO, CoCoO_y, color='blue', linewidth=1.5)
		figure_labels.append('CoCoO Buffer')
		if CoCoO_extrap_lower is not None:
			CoCoO_extrap_lower_plt, = ax1.plot(CoCoO_extrap_lower_C, CoCoO_extrap_lower_y, color='blue', linestyle='dashed', linewidth=1.5)
			figure_labels.append('CoCoO Extrapolated Lower')
		if CoCoO_extrap_upper is not None:
			CoCoO_extrap_upper_plt, = ax1.plot(CoCoO_extrap_upper_C, CoCoO_extrap_upper_y, color='blue', linestyle='dashed', linewidth=1.5)
			figure_labels.append('CoCoO Extrapolated Upper')

	#ReReO
	if 'ReReO' in buffers:
		temp_range_ReReO, temp_range_ReReO_K, ReReO_extrap_lower, ReReO_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 565, 1200)
		if ReReO_extrap_lower is not None:
			ReReO_extrap_lower_C = ReReO_extrap_lower - 273.15
		if ReReO_extrap_upper is not None:
			ReReO_extrap_upper_C = ReReO_extrap_upper - 273.15

		ReReO_y = calc_ReReO(PGPa, temp_range_ReReO_K)
		if ReReO_extrap_lower is not None:
			ReReO_extrap_lower_y = calc_ReReO(PGPa, ReReO_extrap_lower)
		if ReReO_extrap_upper is not None:
			ReReO_extrap_upper_y = calc_ReReO(PGPa, ReReO_extrap_upper)

		ReReO_plt, = ax1.plot(temp_range_ReReO, ReReO_y, color='magenta', linewidth=1.5)
		figure_labels.append('ReReO Buffer')
		if ReReO_extrap_lower is not None:
			ReReO_extrap_lower_plt, = ax1.plot(ReReO_extrap_lower_C, ReReO_extrap_lower_y, color='magenta', linestyle='dashed', linewidth=1.5)
			figure_labels.append('ReReO Extrapolated Lower')
		if ReReO_extrap_upper is not None:
			ReReO_extrap_upper_plt, = ax1.plot(ReReO_extrap_upper_C, ReReO_extrap_upper_y, color='magenta', linestyle='dashed', linewidth=1.5)
			figure_labels.append('ReReO Extrapolated Upper')

	#Graphite
	if 'Graphite' in buffers:
		temp_range_Graphite, temp_range_Graphite_K, Graphite_extrap_lower, Graphite_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 565, 1200)
		if Graphite_extrap_lower is not None:
			Graphite_extrap_lower_C = Graphite_extrap_lower - 273.15
		if Graphite_extrap_upper is not None:
			Graphite_extrap_upper_C = Graphite_extrap_upper - 273.15

		Graphite_y = calc_Graphite(PGPa, temp_range_Graphite_K)
		if Graphite_extrap_lower is not None:
			Graphite_extrap_lower_y = calc_Graphite(PGPa, Graphite_extrap_lower)
		if Graphite_extrap_upper is not None:
			Graphite_extrap_upper_y = calc_Graphite(PGPa, Graphite_extrap_upper)

		Graphite_plt, = ax1.plot(temp_range_Graphite, Graphite_y, color='gray', linewidth=1.5)
		figure_labels.append('Graphite Buffer')
		if Graphite_extrap_lower is not None:
			Graphite_extrap_lower_plt, = ax1.plot(Graphite_extrap_lower_C, Graphite_extrap_lower_y, color='gray', linestyle='dashed', linewidth=1.5)
			figure_labels.append('Graphite Extrapolated Lower')
		if Graphite_extrap_upper is not None:
			Graphite_extrap_upper_plt, = ax1.plot(Graphite_extrap_upper_C, Graphite_extrap_upper_y, color='gray', linestyle='dashed', linewidth=1.5)
			figure_labels.append('Graphite Extrapolated Upper')

	#QIF
	if 'QIF' in buffers:
		#Define QIF low temperature buffer value at P
		temp_range_QIF, temp_range_QIF_K, QIF_extrap_lower, QIF_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 150, 1200)
		if QIF_extrap_lower is not None:
			QIF_extrap_lower_C = QIF_extrap_lower - 273.15
		if QIF_extrap_upper is not None:
			QIF_extrap_upper_C = QIF_extrap_upper - 273.15

		QIF_y = calc_QIF(PGPa, temp_range_QIF_K)
		if QIF_extrap_lower is not None:
			QIF_extrap_lower_y = calc_QIF(PGPa, QIF_extrap_lower)
		if QIF_extrap_upper is not None:
			QIF_extrap_upper_y = calc_QIF(PGPa, QIF_extrap_upper)

		QIF_plt, = ax1.plot(temp_range_QIF, QIF_y, color='mediumaquamarine', linewidth=1.5)
		figure_labels.append('QIF')
		if QIF_extrap_lower is not None:
			QIF_extrap_lower_plt, = ax1.plot(QIF_extrap_lower_C, QIF_extrap_lower_y, color='mediumaquamarine', linestyle='dashed', linewidth=1.5)
			figure_labels.append('QIF Extrapolated')
		if QIF_extrap_upper is not None:
			QIF_extrap_upper_plt, = ax1.plot(QIF_extrap_upper_C, QIF_extrap_upper_y, color='mediumaquamarine', linestyle='dashed', linewidth=1.5)
			figure_labels.append('QIF Extrapolated Upper')

	#SiSiO2
	if 'SiSiO2' in buffers:
		temp_range_SiSiO2, temp_range_SiSiO2_K, SiSiO2_extrap_lower, SiSiO2_extrap_upper = set_calibration_temp_range(temperature_min, temperature_max, 565, 1200)
		if SiSiO2_extrap_lower is not None:
			SiSiO2_extrap_lower_C = SiSiO2_extrap_lower - 273.15
		if SiSiO2_extrap_upper is not None:
			SiSiO2_extrap_upper_C = SiSiO2_extrap_upper - 273.15

		SiSiO2_y = calc_SiSiO2(PGPa, temp_range_SiSiO2_K)
		if SiSiO2_extrap_lower is not None:
			SiSiO2_extrap_lower_y = calc_SiSiO2(PGPa, SiSiO2_extrap_lower)
		if SiSiO2_extrap_upper is not None:
			SiSiO2_extrap_upper_y = calc_SiSiO2(PGPa, SiSiO2_extrap_upper)

		SiSiO2_plt, = ax1.plot(temp_range_SiSiO2, SiSiO2_y, color='purple', linewidth=1.5)
		figure_labels.append('SiSiO2 Buffer')
		if SiSiO2_extrap_lower is not None:
			SiSiO2_extrap_lower_plt, = ax1.plot(SiSiO2_extrap_lower_C, SiSiO2_extrap_lower_y, color='purple', linestyle='dashed', linewidth=1.5)
			figure_labels.append('SiSiO2 Extrapolated Lower')
		if SiSiO2_extrap_upper is not None:
			SiSiO2_extrap_upper_plt, = ax1.plot(SiSiO2_extrap_upper_C, SiSiO2_extrap_upper_y, color='purple', linestyle='dashed', linewidth=1.5)
			figure_labels.append('SiSiO2 Extrapolated Upper')


	#Define legend styles
	labels = [i for i in figure_labels]

	#Label axes
	plt.xlabel('Temperature ($^\circ$C)')
	plt.ylabel('log $f$O$_2$')

	ax1.legend(labels, loc='lower right')
	plt.title('Redox buffers at ' + str(pressure) + ' bar')
	return plt.show()

#TESTING
#plotfO2(1, 700, 1500, temperature_step=1, buffers=['NNO', 'IW'])
