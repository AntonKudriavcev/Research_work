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
skip_number_of_bytes = 4660
dopp_freq_max        = int(5e3) ## Hz
dopp_freq_step       = int(500) ## Hz
threashold           = 2.5
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
									samples_per_Rang_code_chip)) ## num of samples per one 
																	 ## Ranging code chip
samples_for_prosessing = samples_per_ms * ms_to_process
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

def carrier_generator():
	
	carrier = np.sin(2 * np.pi * intermediate_freq * phase_points)

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



rang_code          = ranging_code_generator(1)
rang_code_specrtum = np.conj(np.fft.fft(rang_code))
carrier            = carrier_generator()
signal             = rang_code * carrier

signal_with_noise = np.empty(len(signal))

def acquisition(satellite_num, sigma):

	# rang_code          = ranging_code_generator(satellite_num)
	# rang_code_specrtum = np.conj(np.fft.fft(rang_code))
	# carrier            = carrier_generator()
	# signal             = rang_code * carrier

	# signal_with_noise = np.empty(len(signal))

	for i in range(len(signal)):
		signal_with_noise[i]  = signal[i] + random.gauss(0, sigma)

	for freq_index in range(number_of_frq_bins):

		IQ_signal = (signal_with_noise * 
					 np.exp(-2 * np.pi * 1j * (intermediate_freq - 
					 							dopp_freq_max + 
					 							freq_index * dopp_freq_step) * 
					 							phase_points))

		IQ_signal_specrtum = np.fft.fft(IQ_signal)
		correlation = abs(np.fft.ifft(rang_code_specrtum * IQ_signal_specrtum))**2

		results[freq_index, :] = correlation

	peak_size       = results.max(1).max()
	frequency_index = results.max(1).argmax()

	peak_size  = results.max(0).max()
	code_phase = results.max(0).argmax()

	max_correlation = results[frequency_index, :]


	exclude_Range_Index1 = code_phase - samples_per_Rang_code_chip
	exclude_Range_Index2 = code_phase + samples_per_Rang_code_chip

	# boundaries
	if exclude_Range_Index1 <= 0:
		code_phase_range = np.r_[exclude_Range_Index2 : samples_for_prosessing + exclude_Range_Index1 + 1]

	elif exclude_Range_Index2 >= samples_for_prosessing - 1:
		code_phase_range = np.r_[exclude_Range_Index2 - samples_for_prosessing : exclude_Range_Index1]

	else:
		code_phase_range = np.r_[0 : exclude_Range_Index1 + 1, exclude_Range_Index2 : samples_for_prosessing]

	    # --- Find the second highest correlation peak in the same freq. bin ---
	limited_correlation = max_correlation[code_phase_range]

	second_peak = limited_correlation.max()

	correlation_ratio = peak_size / second_peak

	if (correlation_ratio) > threashold:

		print('the signal is detected')
		print('%.d\t%.d\t%.f' %(satellite_num, (- dopp_freq_max + 
					 							frequency_index * dopp_freq_step), 
												correlation_ratio))

# -----------------------------------------------------------------------------		

		fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows = 4, ncols = 1)
		plt.xlabel('Time')
		plt.ylabel('Amplitude')

		ax1.plot(phase_points[:100], rang_code[:100])
		ax2.plot(phase_points[:100], carrier[:100])
		ax3.plot(phase_points[:100], signal_with_noise[:100])
		ax4.plot(phase_points[:100], (signal_with_noise * carrier)[:100])
		plt.show()

		fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1)

		ax1.plot(max_correlation)
		ax2.plot(limited_correlation)
		plt.show()

##-----------------------------------------------------------------------------		

		return 1
			
	else:
		return 0


if __name__ == '__main__':

	acquisition(1,  1)



	# num_of_repetititons = 100
	
	# sigma_step   = 1
	# sigma_max    = 20

	# num_of_sigma = int(sigma_max / sigma_step)
	# print(num_of_sigma)


	# sigma = (np.arange(num_of_sigma) * sigma_step) ** 1
	# print(sigma)
	# print((np.arange(num_of_sigma) * sigma_step))
	# # E_div_dispersion_array = 1 / sigma_array**2

	# E = 1

	# probabilities = np.zeros(num_of_sigma)
	# # print(probabilities)

	# for i in range(num_of_sigma):

	# 	result = 0

	# 	for j in range(num_of_repetititons):

	# 		result += acquisition(1,  i * sigma_step)

	# 	error_probability = 1 - (result / num_of_repetititons)

	# 	probabilities[i]  = error_probability


	# plt.plot(sigma, probabilities)
	# plt.title('Probability characteristics')
	# plt.xlabel('Sigma')
	# plt.ylabel('P_error')
	# plt.grid()
	# plt.show()







	






