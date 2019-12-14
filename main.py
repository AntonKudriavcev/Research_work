##=============================================================================
##
##=============================================================================
##
## Program for research probability characteristics 
## of 1I signal Beidou navigation system
##
##-----------------------------------------------------------------------------
##
## Creator: Kudriavcev Anton
## email  : Kudriavcev.Anton@yandex.ru 
##
##=============================================================================
##---------------------------------IMPORTS-------------------------------------
##=============================================================================

import numpy as np
import random
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

from graph_builder import Builder

##=============================================================================
##-------------------------------PARAMETERS------------------------------------
##=============================================================================

##------------------------------file parameters--------------------------------

file_name      = 'D://study//5_year//Course_project//SoftGNSS_python//data.dat'
data_type      = 'int8'
signal_perform = 'complex'
sampling_freq  = 10e6 ## Hz

##---------------------------simulation parameters-----------------------------

ms_to_process        = 1
samples_to_shift     = 1000
skip_number_of_bytes = 4660
dopp_freq_max        = int(5e3) ## Hz
dopp_freq_step       = int(500) ## Hz
# threashold           = 2.0
num_of_satellites    = 37

error                = 0

##-----------------------------signal parameters-------------------------------

intermediate_freq = 1.098e6 ## Hz

## Ranging Code

Rang_code_length = 2046
Rang_code_freq   = 2.046e6 ## Hz

G2_shift = {
   1: (1,3),    2: (1,4),    3: (1,5),    4: (1,6),
   5: (1,8),    6: (1,9),    7: (1,10),   8: (1,11),
   9: (2,7),   10: (3,4),   11: (3,5),   12: (3,6),
  13: (3,8),   14: (3,9),   15: (3,10),  16: (3,11),
  17: (4,5),   18: (4,6),   19: (4,8),   20: (4,9),
  21: (4,10),  22: (4,11),  23: (5,6),   24: (5,8),
  25: (5,9),   26: (5,10),  27: (5,11),  28: (6,8),
  29: (6,9),   30: (6,10),  31: (6,11),  32: (8,9),
  33: (8,10),  34: (8,11),  35: (9,10),  36: (9,11),
  37: (10,11)
}

## Navigation message 

navig_mess_freq = 50 ## Hz

## Neumann-Hoffman code

NH_code      = np.array([1, 1, 1, 1, 1,-1, 1, 1,-1,-1, 
						 1,-1, 1,-1, 1, 1,-1,-1,-1, 1])
NH_code_freq = 1e3 ## Hz

##=============================================================================
##--------------------------SIMULATION VARIABLES-------------------------------
##=============================================================================

time_of_sample         = 1 / sampling_freq
time_of_Rang_code_chip = 1 / Rang_code_freq
samples_per_ms         = int(round(sampling_freq / 1e3))

samples_per_Rang_code_chip = int(round(sampling_freq / Rang_code_freq))

samples_per_Rang_code = np.longlong(round(Rang_code_length * 
									sampling_freq / Rang_code_freq)) ## num of samples per one 
																	 ## Ranging code chip
samples_for_prosessing = int(round(samples_per_ms * ms_to_process))

phase_points           = np.arange(samples_for_prosessing) * time_of_sample
Ranging_code_period    = Rang_code_length / Rang_code_freq
Ranging_code_index     = np.longlong((np.floor(np.arange(0, samples_for_prosessing) * 
									Rang_code_freq / sampling_freq))) ## number of whole Rang_code_chips
																		      ## per one sample
															      
Ranging_code_index = Ranging_code_index % Rang_code_length	
															      
number_of_frq_bins = np.int(np.floor(2 * dopp_freq_max / dopp_freq_step) + 1)

results = np.zeros((number_of_frq_bins, samples_for_prosessing))

##=============================================================================
##-------------------------------FUNCTIONS-------------------------------------
##=============================================================================

def carrier_generator(sigma):
	
	carrier = np.sin(2 * np.pi * (intermediate_freq + random.gauss(0, sigma)) * phase_points)

	return carrier


def ranging_code_generator(satellite_num):

    g2_shift = G2_shift[satellite_num]

    g1 = np.zeros(Rang_code_length)
    g2 = np.zeros(Rang_code_length)

##---------------------------------for g1 sequence-----------------------------

    reg = np.array([1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1]) ## init g1 register 
                                                           ## as 01010101010
    for i in range(Rang_code_length):

        g1[i] = reg[-1]
        save_Bit = reg[0] * reg[6] * reg[7] * reg[8] * reg[9] * reg[10]
        reg[1:] = reg[:-1]
        reg[0] = save_Bit

##---------------------------------for g2 sequence-----------------------------

    reg = np.array([1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1]) ## init g2 register 
                                                           ## as 01010101010
    for i in range(Rang_code_length):

        g2_shift_0 = g2_shift[0] - 1
        g2_shift_1 = g2_shift[1] - 1

        g2[i]    = reg[g2_shift_0] * reg[g2_shift_1] ## XOR for different phase coefficients
        save_Bit = reg[0] * reg[1] * reg[2] * reg[3] * reg[4] * reg[7] * reg[8] * reg[10]
        reg[1:]  = reg[:-1]
        reg[0]   = save_Bit

##---------------------------the resulting Ranging code------------------------

    G1G2 = g1 * g2

    # print(len(G1G2))

    Ranging_code = G1G2[Ranging_code_index] ## resampling Ranging code
    # print(Ranging_code)
    # print(len(Ranging_code))


    return Ranging_code

def navigation_message_genrator():
	pass

def nh_code_generator():
	pass

def signal_generator():

	rang_codes          = np.zeros((num_of_satellites + 1, samples_for_prosessing))
	rang_codes_specrtum = np.zeros((num_of_satellites + 1, samples_for_prosessing), dtype = np.complex128)
	signals             = np.zeros((num_of_satellites + 1, samples_for_prosessing))

	carrier             = carrier_generator(sigma = 0)


	for satellit in range(1, num_of_satellites + 1): 
		rang_codes[satellit]          = ranging_code_generator(satellit)
		rang_codes_specrtum[satellit] = np.conj(np.fft.fft(rang_codes[satellit]))
		signals[satellit]             = rang_codes[satellit] * carrier

		signal1 = signals[satellit] [0 : samples_to_shift]
		signal2 = signals[satellit] [samples_to_shift : ]

		signals[satellit] = np.concatenate([signal2, signal1])


	return rang_codes, rang_codes_specrtum, carrier, signals



def acquisition(rang_codes, rang_codes_specrtum, carrier, signals, satellite_num, sigma, amplitude, Threashold):

	rang_code = rang_codes[satellite_num]
	rang_code_specrtum = rang_codes_specrtum[satellite_num]
	signal = signals[satellite_num]
	
	signal_with_noise = np.empty(len(signal))
##-----------------------------------------------------------------------------

	# plot_builder = Builder(num_of_satellites = satellite_num,
 #                           sampl_freq       = sampling_freq, 
 #                           samples_per_code = samples_for_prosessing,
 #                           sigma = sigma)

##-----------------------------------------------------------------------------

	for i in range(len(signal)):
		signal_with_noise[i]  = amplitude * signal[i] + random.gauss(0, sigma)

	for freq_index in range(number_of_frq_bins):

		IQ_signal = (signal_with_noise * 
					 np.exp(-2 * np.pi * 1j * (intermediate_freq - 
					 							dopp_freq_max + 
					 							freq_index * dopp_freq_step) * 
					 							phase_points))

		IQ_signal_specrtum = np.fft.fft(IQ_signal)
		correlation = abs(np.fft.ifft(rang_code_specrtum * IQ_signal_specrtum))**2

		results[freq_index, :] = correlation
##-----------------------------------------------------------------------------
		# plot_builder.add_to_plot(data           = results[freq_index, :], 
  #                                freq_deviation = (-dopp_freq_max + 
		# 			 								freq_index * dopp_freq_step))

	# plot_builder.show_plot()
##-----------------------------------------------------------------------------

	peak_size       = results.max(1).max()
	frequency_index = results.max(1).argmax()

	peak_size  = results.max(0).max()
	code_phase = results.max(0).argmax()

	max_correlation = results[frequency_index, :]


	exclude_Range_Index1 = code_phase - samples_per_Rang_code_chip
	exclude_Range_Index2 = code_phase + samples_per_Rang_code_chip

	# boundaries
	if exclude_Range_Index1 <= 0:
		code_phase_range = np.r_[exclude_Range_Index2 : samples_for_prosessing + exclude_Range_Index1]

	elif exclude_Range_Index2 >= samples_for_prosessing - 1:
		code_phase_range = np.r_[exclude_Range_Index2 - samples_for_prosessing : exclude_Range_Index1]
	else:
		code_phase_range = np.r_[0 : exclude_Range_Index1 + 1, exclude_Range_Index2 : samples_for_prosessing]

	# --- Find the second highest correlation peak in the same freq. bin ---

	limited_correlation = max_correlation[code_phase_range]

	second_peak = limited_correlation.max()

	correlation_ratio = peak_size / second_peak

	if (correlation_ratio) > Threashold:

		# print('The signal is detected')
		# print('%.d\t%.d\t%.f' %(satellite_num, (- dopp_freq_max + 
		# 			 							frequency_index * dopp_freq_step), 
		# 										correlation_ratio))

# ##-----------------------------------------------------------------------------		

		# fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows = 4, ncols = 1)
		# plt.xlabel('Time')
		# plt.ylabel('Amplitude')

		# ax1.plot(phase_points, rang_code)
		# ax2.plot(phase_points, carrier)
		# ax3.plot(phase_points, signal_with_noise)
		# ax4.plot(phase_points, (signal_with_noise * carrier))
		# plt.show()

		# fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1)

		# ax1.plot(max_correlation)
		# ax2.plot(limited_correlation)
		# plt.show()

# ##-----------------------------------------------------------------------------		

		return 1
	else:
		return 0

def acquisition_for_all_satellites():
	pass

def plot_prob_acqu_char(rang_codes, rang_codes_specrtum, carrier, signals, sigma,  
						num_of_repetititons, A_step, A_min, A_max, satellite_num, Threashold):

	amplitude = (np.arange(A_min, A_max + A_step, A_step))
	print(amplitude)

	Power = (amplitude**2) / 2
	E_div_N0_array = (sampling_freq / 1) * (Power) / (sigma**2)

	print(E_div_N0_array)

	probabilities = np.zeros(len(amplitude))
	print(probabilities)

	for i in range(len(amplitude)):

		result = 0

		for j in range(num_of_repetititons):

			result += acquisition(rang_codes, rang_codes_specrtum, carrier, signals,
													satellite_num,  sigma, amplitude[i], Threashold)

		acqu_probability = (result / num_of_repetititons)

		probabilities[i]  = acqu_probability

	plt.plot((10 * np.log10(E_div_N0_array)), probabilities)
	plt.title('Probability characteristics, treashold = %.2f' %Threashold)
	plt.xlabel('P/N0, dB*Hz')
	plt.ylabel('P_acquision')
	plt.grid()
	plt.show()

def plot_prob_acqu_char_dB(rang_codes, rang_codes_specrtum, carrier, signals, sigma,
							snrDb_min, snrDb_max, snrDb_step, 
							num_of_repetititons, satellite_num, Threashold):

	snrDb = np.arange(snrDb_min, snrDb_max, snrDb_step)

	snr        = 10**(0.1*snrDb)              # Отношение сигнал-шум в разах
	var        = sigma**2                     # Дисперсия шума в одной квадратуре
	amplitude = np.sqrt(2*snr*var/sampling_freq)        # Амплитуда сигнала

	probabilities = np.zeros(len(snrDb))

	for i in range(len(snrDb)):

		result = 0

		for j in range(num_of_repetititons):

			result += acquisition(rang_codes, rang_codes_specrtum, carrier, 
									signals, satellite_num, sigma, amplitude[i], Threashold)

		acqu_probability = (result / num_of_repetititons)

		probabilities[i]  = acqu_probability

	print(time_start - time.time())


	plt.plot(snrDb, probabilities)
	plt.title('Probability characteristics, treashold = %.2f' %Threashold)
	plt.xlabel('P/N0, dB*Hz')
	plt.ylabel('P_acquision')
	plt.grid()
	# plt.show()

	plt.savefig(fname = ('Acquisition characteristics(step = %.1f, snrDb_min = %.1f, snrDb_max = %.1f , repetition = %.1f, threashold = %.1f).jpg'
							%(snrDb_step, snrDb_min, snrDb_max, num_of_repetititons, Threashold)))
	plt.clf()

def plot_prob_false_char(rang_codes, rang_codes_specrtum, carrier, signals,  
						num_of_repetititons, sigma_step, sigma_min, sigma_max, satellite_num, Threashold):


	sigma = np.arange(sigma_min, sigma_max + sigma_step, sigma_step, dtype = np.float64)
	print(sigma)

	Dispersion = sigma**2
	print(Dispersion)

	one_div_N0_array = (sampling_freq / 1) / Dispersion
	print(one_div_N0_array)


	probabilities = np.zeros(len(sigma))
	print(probabilities)

	for i in range(len(sigma)):

		result = 0

		for j in range(num_of_repetititons):

			result += acquisition(rang_codes, rang_codes_specrtum, carrier, signals,
													satellite_num,  sigma[i], 0, Threashold)

		error_probability = (result / num_of_repetititons)

		probabilities[i]  = error_probability

	plt.plot((10 * np.log10(one_div_N0_array)[::-1]), probabilities[::-1])
	plt.title('Probability characteristics, treashold = %.2f' %Threashold)
	plt.xlabel('1/N0, dB*Hz')
	plt.ylabel('P_error')
	plt.grid()
	plt.savefig(fname = ('False characteristics(step = %.1f, sigma_min = %.1f, sigma_max = %.1f , repetition = %.1f, threashold = %.1f).jpg'
							%(sigma_step, sigma_min, sigma_max, num_of_repetititons, Threashold)))
	plt.clf()
	
if __name__ == '__main__':

	rang_codes, rang_codes_specrtum, carrier, signals = signal_generator()

##--Probability acquisition----------------------------------------------------
	# time_start = time.time()

	# plot_prob_acqu_char_dB(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1,
	# 						snrDb_min = 25, snrDb_max = 55, snrDb_step = 0.5, 
	# 						num_of_repetititons = 130, satellite_num = 1, Threashold = 1.8)

	# plot_prob_acqu_char_dB(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1,
	# 						snrDb_min = 25, snrDb_max = 55, snrDb_step = 0.5, 
	# 						num_of_repetititons = 130, satellite_num = 1, Threashold = 2.0)

	# plot_prob_acqu_char_dB(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1,
	# 						snrDb_min = 25, snrDb_max = 55, snrDb_step = 0.5, 
	# 						num_of_repetititons = 130, satellite_num = 1, Threashold = 2.2)

	# plot_prob_acqu_char_dB(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1,
	# 						snrDb_min = 25, snrDb_max = 55, snrDb_step = 0.5, 
	# 						num_of_repetititons = 130, satellite_num = 1, Threashold = 2.5)

	# plot_prob_acqu_char_dB(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1,
	# 						snrDb_min = 25, snrDb_max = 55, snrDb_step = 0.5, 
	# 						num_of_repetititons = 130, satellite_num = 1, Threashold = 3.0)


	# plot_prob_acqu_char(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1, 
	# 	num_of_repetititons = 130, A_step = 0.009, A_min = 0.005, A_max = 0.2, satellite_num = 1, Threashold = 1.8)

	# plot_prob_acqu_char(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1, 
	# 	num_of_repetititons = 130, A_step = 0.009, A_min = 0.005, A_max = 0.2, satellite_num = 1, Threashold = 2.0)

	# plot_prob_acqu_char(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1, 
	# 	num_of_repetititons = 130, A_step = 0.009, A_min = 0.005, A_max = 0.2, satellite_num = 1, Threashold = 2.2)

	# plot_prob_acqu_char(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1, 
	# 	num_of_repetititons = 40, A_step = 0.009, A_min = 0.005, A_max = 0.2, satellite_num = 1, Threashold = 2.5)

	# plot_prob_acqu_char(rang_codes, rang_codes_specrtum, carrier, signals, sigma = 1, 
	# 	num_of_repetititons = 130, A_step = 0.009, A_min = 0.005, A_max = 0.2, satellite_num = 1, Threashold = 3.0)

##--Probability false----------------------------------------------------------

	plot_prob_false_char(rang_codes, rang_codes_specrtum, carrier, signals, 
		num_of_repetititons = 1000, sigma_step = 10, sigma_min = 5, sigma_max = 250, satellite_num = 1, Threashold = 1.8)

	plot_prob_false_char(rang_codes, rang_codes_specrtum, carrier, signals, 
		num_of_repetititons = 1000, sigma_step = 10, sigma_min = 5, sigma_max = 250, satellite_num = 1, Threashold = 2.0)

	plot_prob_false_char(rang_codes, rang_codes_specrtum, carrier, signals, 
		num_of_repetititons = 1000, sigma_step = 10, sigma_min = 5, sigma_max = 250, satellite_num = 1, Threashold = 2.2)

	plot_prob_false_char(rang_codes, rang_codes_specrtum, carrier, signals, 
		num_of_repetititons = 1000, sigma_step = 10, sigma_min = 5, sigma_max = 250, satellite_num = 1, Threashold = 2.5)

	plot_prob_false_char(rang_codes, rang_codes_specrtum, carrier, signals, 
		num_of_repetititons = 1000, sigma_step = 10, sigma_min = 5, sigma_max = 250, satellite_num = 1, Threashold = 3.0)

##--acquisition----------------------------------------------------------------

	# acquisition(rang_codes, rang_codes_specrtum, carrier, signals, satellite_num = 1, sigma = 1, amplitude = 1, Threashold = 2.5)




















	






