import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from math import *
import random

# -------------CORE DEFINITIONS----------- #
class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

#-------------DEFINE SOME UNIVERSAL FUNCTIONS------------#
def get_curve_color(buffer):
	"""
	Return prefered color for plots for any buffer

	Parameters
	----------
	buffer: str
		Name of buffer

	Returns
	-------
	str
		color name
	"""
	color = {'NNO': 		'green',
			 'QFM': 		'darkorange',
			 'IW':  		'black',
			 'HM':    		'red',
			 'CoCoO': 		'blue',
			 'ReReO': 		'magenta',
			 'Graphite':  	'gray',
			 'QIF': 		'mediumaquamarine',
			 'SiSiO2':		'purple',
			 'CrCr2O3':		'teal',
			 'MoMoO2':		'olive',
			 'CaCaO':		'peru',
			 'AlAl2O3':		'chartreuse',
			 'KK2O':		'deeppink',
			 'MgMgO':		'maroon',
			 'MnMnO':		'midnightblue',
			 'NaNa2O':		'dodgerblue',
			 'TiTiO2':		'orangered'}

	return color[buffer]

def get_calibration_temp_range(buffer, units='C'):
	"""
	Get upper and lower T bounds for which a buffer is calibrated
	=============================================================

	Parameters
	----------
	buffer: str
		Name of buffer

	units: str
		Can be 'K' or 'C'. Specifies which units to return.

	Returns
	-------
	tuple 
		(lower_bound, upper_bound)
	"""
	if units == 'K' or units == 'C':
		pass 
	else:
		raise InputError('Units must be one of "K" or "C".')

	temp_ranges_C = {'NNO': 		(600, 1200),
				     'QFM': 		(400, 1200),
				     'IW':  		(565, 1200),
				     'HM':    		(300, 1100),
				     'CoCoO': 		(600, 1200),
				     'ReReO': 		(565, 1200),
				     'Graphite':  	(565, 1200),
				     'QIF': 		(150, 1200),
				     'SiSiO2':		(565, 1200),
				     'CrCr2O3':		(327, 1527),
				     'MoMoO2':		(951, 1311),
				     'CaCaO':		(25,  1527),
				     'AlAl2O3':		(25,  1527),
				     'KK2O':		(25,  1527),
				     'MgMgO':		(25,  1527),
				     'MnMnO':		(25,  1527),
				     'NaNa2O':		(25,  1527),
				     'TiTiO2':		(25,  1527)}

	temp_ranges_K = {k: (v[0]+273, v[1]+273) for k, v in temp_ranges_C.items()}

	if units == 'C':
		return temp_ranges_C[buffer]
	elif units == 'K':
		return temp_ranges_K[buffer]
	
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
def calc_buffer(P, T, buffer):
	"""
	Master function to calc any buffer given a name.

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	buffer: str
		Name of buffer

	Returns
	-------
	float or numpy array
		logfO2
	"""

	if buffer == 'NNO': 
		return calc_NNO(P, T)
	elif buffer == 'QFM': 
		return calc_QFM(P, T)
	elif buffer == 'IW': 
		return calc_IW(P, T)
	elif buffer == 'CrCr2O3': 
		return calc_CrCr2O3(P, T)
	elif buffer == 'SiSiO2': 
		return calc_SiSiO2(P, T)
	elif buffer == 'HM': 
		return calc_HM(P, T)
	elif buffer == 'CoCoO': 
		return calc_CoCoO(P, T)
	elif buffer == 'ReReO': 
		return calc_ReReO(P, T)
	elif buffer == 'Graphite': 
		return calc_Graphite(P, T)
	elif buffer == 'QIF': 
		return calc_QIF(P, T)
	elif buffer == 'MoMoO2':
		return calc_MoMoO2(P,T)
	elif buffer == 'CaCaO':
		return calc_CaCaO(P,T)
	elif buffer == 'AlAl2O3':
		return calc_AlAl2O3(P,T)
	elif buffer == 'KK2O':
		return calc_KK2O(P,T)
	elif buffer == 'MgMgO':
		return calc_MgMgO(P,T)
	elif buffer == 'MnMnO':
		return calc_MnMnO(P,T)
	elif buffer == 'NaNa2O':
		return calc_NaNa2O(P,T)
	elif buffer == 'TiTiO2':
		return calc_TiTiO2(P,T)
	else:
		raise InputError('Buffer name not recognized')


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
		logfO2

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
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = (8.699 + 0.01642*P -0.0003*P**2 + (2.7*10**(-6))*P**3 - (10**(-8))*P**4) + (-24205 + 444.73*P - 0.5929*P**2 + 0.00153*P**3)/T

	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append((8.699 + 0.01642*P -0.0003*P**2 + (2.7*10**(-6))*P**3 - 
								(10**(-8))*P**4) + (-24205 + 444.73*P - 0.5929*P**2 + 
								0.00153*P**3)/temp)
		log_fO2 = np.array(log_fO2_list)

	return log_fO2

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
		logfO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	# translate P to bars
	P_bars = P*10000

	if isinstance(T, float) or isinstance(T, int):
		if T<573:
			log_fO2 = (-26445.3/T) + 10.344 + 0.092 * (P_bars-1)/T
		if T>=573:
			log_fO2 = (-25096.3/T) + 8.735 + 0.11 * (P_bars-1)/T

	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			if temp<573:
				log_fO2_list.append((-26445.3/temp) + 10.344 + 0.092 * (P_bars-1)/temp)
			if temp>=573:
				log_fO2_list.append((-25096.3/temp) + 8.735 + 0.11 * (P_bars-1)/temp)
		log_fO2 = np.array(log_fO2_list)

	return log_fO2

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
		log_fO2

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
	log_fO2 = (6.54106+0.0012324*P) + (-28163.6+546.32*P-1.13412*P**2+0.0019274*P**3)/T
	return log_fO2

def calc_CrCr2O3(P, T):
	"""
	Cr-Cr2O3
	========
	Define the chromium-chromium oxide buffer value (calibrated only for 1 bar)
	Valid between 600–1800 K

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Holzheid and O'Neill (1994) The Cr-Cr2O3 oxygen bufer and the free
	energy of formation of Cr2O3 from high-temperature electrochemical
	measurements, Geochimica et Cosmochimica Acta 59, pp. 475–497k.
	"""
	if isinstance(T, int) or isinstance(T, float):
		log_fO2 = log10(exp((-758585 + 349.718*T - 25.5216*T*log(T) + 0.00935*T**2)/
						(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2 = log10(exp((-758585 + 349.718*temp - 25.5216*temp*log(temp) + 0.00935*temp**2)/
						(8.314*temp)))
			log_fO2_list.append(log_fO2)
		log_fO2 = np.array(log_fO2_list)
	return log_fO2


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
		log_fO2

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
		log_fO2_1bar = log10(exp((-911336 + 212.2423*T - 4.512*T*log(T))/(8.314*T)))
		log_fO2 = (1/T*0.4343*120*((0.00251*P**3 - 0.2044*P**2 + 10.32257*P)-(0.00251*0.0001**3 - 0.2004*0.0001**2 + 10.32257*0.001)))+log_fO2_1bar
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_1bar = log10(exp((-911336 + 212.2423*temp - 4.512*temp*log(temp))/(8.314*temp)))
			log_fO2_list.append((1/temp*0.4343*120*((0.00251*P**3 - 0.2044*P**2 + 10.32257*P)-(0.00251*0.0001**3 - 0.2004*0.0001**2 + 10.32257*0.001)))+log_fO2_1bar)
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

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
		log_fO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	# translate P to bars
	P_bars = P*10000

	if isinstance(T, float) or isinstance(T, int):
		if T<=573.0:
			log_fO2 = (-25497.5/T) + 14.330 + 0.019 * (P_bars-1)/T
		if T>573.0 and T<682:
			log_fO2 = (-26452.6/T) + 15.455 + 0.019 * (P_bars-1)/T
		if T>=682:
			log_fO2 = (-25700.6/T) + 14.558 + 0.019 * (P_bars-1)/T
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			if temp<=573:
				log_fO2_list.append((-25497.5/temp) + 14.330 + 0.019 * (P_bars-1)/temp)
			if temp>573 and temp<682:
				log_fO2_list.append((-26452.6/temp) + 15.455 + 0.019 * (P_bars-1)/temp)
			if temp>=682:
				log_fO2_list.append((-25700.6/temp) + 14.558 + 0.019 * (P_bars-1)/temp)
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

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
		log_fO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	# translate P to bars
	P_bars = P*10000

	log_fO2 = (-24332.6/T) + 7.295 + 0.052 * (P_bars-1)/T
	return log_fO2

def calc_ReReO(P, T):
	"""
	Rhenium-Rhenium Oxide (ReReO)
	=============================
	Define ReReO buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Pownceby and O'Neill (1994) Thermodynamic data from redox reactions at high temperatures. 
	IV. Calibration of the Re-ReO2 oxygen buffer from EMF and NiO+Ni-Pd redox sensor measurements
	"""
	log_fO2 = (-451020 + 297.595 * T - 14.6585 * T * np.log(T))/(8.31441 * T * log(10))
	return log_fO2

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
		log_fO2

	References
	----------
	French and Eugster (1965) Journal of Geophysical Research
	"""
	log_fO2 = (-20586/T) - 0.044 + np.log10(P) - 0.028 * (P-1)/T
	return log_fO2

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
		log_fO2

	References
	----------
	B. R. Frost in Mineralogical Society of America "Reviews in Mineralogy" Volume 25
	"""
	# translate P to bars
	P_bars = P*10000
	
	if isinstance(T, float) or isinstance(T, int):
		if T<=573.0:
			log_fO2 = (-29435.7/T) + 7.391 + 0.44 * (P_bars-1)/T
		if T>573.0:
			log_fO2 = (-29520.8/T) + 7.492 + 0.050 * (P_bars-1)/T
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			if temp<573:
				log_fO2_list.append((-29435.7/temp) + 7.391 + 0.44 * (P_bars-1)/temp)
			if temp>=573:
				log_fO2_list.append((-29520.8/temp) + 7.492 + 0.050 * (P_bars-1)/temp)
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_MoMoO2(P, T):
	"""
	Molybdenum-Molybdenum Oxide (MoMoO2)
	====================================
	Define MoMoO2 buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	J. Bygdén, Du Sichen, and S. Seetharaman (1994) "A Thermodynamic Study of
	the Molybdenum-Oxygen System" Metallurgical and Materials 
	Transactions B, vol 25B
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = log10(exp((-580563+173*T)/(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(log10(exp((-580563+173*temp)/(8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_CaCaO(P, T):
	"""
	Calcium-Lime (CaCaO)
	====================
	Define CaCaO buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Robie et al. (1979) Thermodynamic properties of minerals and related
	substances at 298.15 K (25 °C) and one bar (10**5 Pascals) pressure
	and at higher temperature
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = 2*log10(exp((-674401+373.5116*T-32.9484*T*log(T))/(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(2*log10(exp((-674401+373.5116*temp-32.9484*temp*log(temp))/(8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_AlAl2O3(P, T):
	"""
	Aluminum-Aluminum Oxide (Al-Al2O3)
	==================================
	Define AlAl2O3 buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Robie et al. (1979) Thermodynamic properties of minerals and related
	substances at 298.15 K (25 °C) and one bar (10**5 Pascals) pressure
	and at higher temperature
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = 0.666 * log10(exp((-1673251+255.6191*T + 
									 8.39424*T*log(T))/(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(0.666 * 
								log10(exp((-1673251+255.6191*temp + 
									 	   8.39424*temp*log(temp))/(8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_KK2O(P, T):
	"""
	Potassium-Potassium Oxide (K-K2O)
	=================================
	Define KK2O buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Robie et al. (1979) Thermodynamic properties of minerals and related
	substances at 298.15 K (25 °C) and one bar (10**5 Pascals) pressure
	and at higher temperature
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = 2 * log10(exp((-301008-705.086*T + 
								 114.6013*T*log(T)) / 
								(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(2 * log10(exp((-301008 - 705.086*temp + 
								 				114.6013*temp*log(temp)) / 
											  (8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_MgMgO(P, T):
	"""
	Magnesium-Periclase (Mg-MgO)
	=================================
	Define MgMgO buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Robie et al. (1979) Thermodynamic properties of minerals and related
	substances at 298.15 K (25 °C) and one bar (10**5 Pascals) pressure
	and at higher temperature
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = 2 * log10(exp((-552747-465.738*T +
								 75.88778*T*log(T)) /
								(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(2 * log10(exp((-552747-465.738*temp +
								 			   75.88778*temp*log(temp)) /
											  (8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_MnMnO(P, T):
	"""
	Manganese-Manganese Oxide (Mn-MnO)
	==================================
	Define MnMnO buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Robie et al. (1979) Thermodynamic properties of minerals and related
	substances at 298.15 K (25 °C) and one bar (10**5 Pascals) pressure
	and at higher temperature
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = 2 *log10(exp((-376969-4.19123*T +
							10.02409*T*log(T)) /
							(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(2*log10(exp((-376969-4.19123*temp +
											10.02409*temp*log(temp)) /
										  (8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_NaNa2O(P, T):
	"""
	Sodium-Sodium Oxide (Na-Na2O)
	=============================
	Define NaNa2O buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Robie et al. (1979) Thermodynamic properties of minerals and related
	substances at 298.15 K (25 °C) and one bar (10**5 Pascals) pressure
	and at higher temperature
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = 2 *log10(exp((-335840-905.443*T +
								139.3508*T*log(T)) /
								(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(2 *log10(exp((-335840-905.443*temp +
											  139.3508*temp*log(temp)) /
											 (8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

def calc_TiTiO2(P, T):
	"""
	Titanium-Titanium Oxide (Ti-TiO2)
	================================
	Define TiTiO2 buffer value at 1 bar

	Parameters
	----------
	P: float
		Pressure in GPa

	T: float or numpy array
		Temperature in degrees K

	Returns
	-------
	float or numpy array
		log_fO2

	References
	----------
	Barin (1993) Thermo database
	"""
	if isinstance(T, float) or isinstance(T, int):
		log_fO2 = log10(exp((-945822 + 219.6816*T -
							 5.25733*T*log(T)) /
							(8.314*T)))
	if isinstance(T, np.ndarray):
		log_fO2_list = []
		for temp in T:
			log_fO2_list.append(log10(exp((-945822 + 219.6816*temp -
							 			   5.25733*temp*log(temp)) /
										  (8.314*temp))))
		log_fO2 = np.array(log_fO2_list)
	return log_fO2

#--------------PLOTTING-------------#
def plot_log_fO2(pressure, temperature_min, temperature_max, temperature_step=1, buffers=['NNO', 'QFM']):
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

	# DRAW THE FIGURE
	fig, ax1 = plt.subplots()
	figure_labels = []

	for buffer_name in buffers:
		curve_color = get_curve_color(buffer_name)
		cal_temp_C = get_calibration_temp_range(buffer_name, 'C')
		cal_temp_K = get_calibration_temp_range(buffer_name, 'K')

		buffer_usable_temp_range_C, buffer_usable_temp_range_K, buffer_extrap_lower, buffer_extrap_upper =  set_calibration_temp_range(temperature_min,
																																	   temperature_max,
																																	   cal_temp_C[0],
																																	   cal_temp_C[1])
		if buffer_extrap_lower is not None:
			buffer_extrap_lower_C = buffer_extrap_lower - 273.15
		if buffer_extrap_upper is not None:
			buffer_extrap_upper_C = buffer_extrap_upper - 273.15

		y = calc_buffer(PGPa, buffer_usable_temp_range_K, buffer_name)
		if buffer_extrap_lower is not None:
			buffer_extrap_lower_y = calc_buffer(PGPa, buffer_extrap_lower, buffer_name)
		if buffer_extrap_upper is not None:
			buffer_extrap_upper_y = calc_buffer(PGPa, buffer_extrap_upper, buffer_name)

		# plot the curves
		# plot the non-extrapolated portion
		the_plt, = ax1.plot(buffer_usable_temp_range_C, y, color=curve_color, linewidth=1.5)
		figure_labels.append(buffer_name)

		# if applicable, plot the lower extrapolation
		if buffer_extrap_lower is not None:
			extrap_lower_plt, = ax1.plot(buffer_extrap_lower_C, buffer_extrap_lower_y, 
										 color=curve_color, linestyle='dashed',
										 linewidth=1.5)
			figure_labels.append('_nolegend_')

		# if applicable, plot the upper extrapolation
		if buffer_extrap_upper is not None:
			extrap_upper_plt, = ax1.plot(buffer_extrap_upper_C, buffer_extrap_upper_y, 
										 color=curve_color, linestyle='dashed',
										 linewidth=1.5)
			figure_labels.append('_nolegend_')

		# TODO
		# # add text annotation directly on plotted curve if desired
		# min_text_loc = 0.1 * (temperature_max-temperature_min) + temperature_min
		# max_text_loc = temperature_max - 0.1 * (temperature_max-temperature_min)

		# # create place holders
		# xloc_list = []
		# yloc_list = []

		# # create new location
		# xloc = random.randint(min_text_loc, max_text_loc)
		# yloc = calc_buffer(PGPa, xloc+273, buffer_name)

		# # TODO ensure this loc isn't too close to one already used

		# plt.annotate(buffer_name, xy=(xloc, yloc), xycoords='data')

	# Define legend styles
	labels = [i for i in figure_labels]

	# Set x limits based on user inputs
	ax1.set_xlim([temperature_min, temperature_max])

	# Label axes
	plt.xlabel('Temperature ($^\circ$C)')
	plt.ylabel('log $f$O$_2$')

	ax1.legend(labels, loc='lower right')
	plt.title('Redox buffers at ' + str(pressure) + ' bar')
	return plt.show()

#TESTING
#plotlog_fO2(1, 700, 1500, temperature_step=1, buffers=['NNO', 'IW'])
