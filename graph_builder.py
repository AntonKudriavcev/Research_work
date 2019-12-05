

import numpy as np
import pylab
from mpl_toolkits.mplot3d import Axes3D

class Builder():

    def __init__(self, num_of_satellites, sampl_freq, samples_per_code, sigma):

        self.sampl_freq       = sampl_freq
        self.samples_per_code = samples_per_code
        self.ts               = 1.0 / sampl_freq
        self.phase_points     = np.arange(samples_per_code) * self.ts

        self.fig = pylab.figure()
        self.ax  = self.fig.add_subplot(111, projection='3d')
        pylab.title('Ambiguity function for %.d satellit, sigma = %.1f' %(num_of_satellites, sigma))
        pylab.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        self.ax.set_zlabel('Correlation coefficient')
        self.ax.set_ylabel('Time')
        self.ax.set_xlabel('Frequency')

    def add_to_plot(self, data, freq_deviation):

        self.ax.plot(freq_deviation * np.ones(len(self.phase_points)), ## create array for certain frequency deviation
                 	 self.phase_points, 
                 	 data)

    def show_plot(self):

        pylab.show()


