import matplotlib
import matplotlib.pyplot as plt
import pylab as p
import numpy as np
from filter import *



###
# Plotting the acceleration time histories 		#
# and stress-strain history with backbone curve #
###



def readfile (filename, coeffx, coeffy):
	"""This little subroutine reads the filename file
	and returns the first two columns of the file by 
	multiplying them with given coefficients """

	print 'Reading the file', filename
	print '...'

	x = np.genfromtxt(filename, usecols=0)* coeffx
	y = np.genfromtxt(filename, usecols=1)* coeffy

	return x, y



def plot_them (time,accel,gamma,sigma,BBgamma,BBsigma):

	# Plotting
	fig = p.figure(num=1,figsize=(12,8),dpi=300,facecolor=None,edgecolor='k'); 
	fig.subplots_adjust(left=0.15,right=0.90,bottom=0.1,top=0.90,wspace=0.2,hspace=0.2)

	# sns.set_style(style='white')
	# sns.set_context('poster')

	# Acceleration vs time
	ax = fig.add_subplot(121)
	ax.plot(time, accel, lw=1.0, color = 'red', label = 'Surface' )

	ax.set_xlabel('Time [s]', fontsize=18)
	ax.set_ylabel('Acceleration [m/s/s]', fontsize=18)


	# Stress-strain curve
	ax = fig.add_subplot(122)

	ax.set_xlim([-5e-2,5e-2])
	ax.plot(BBgamma, BBsigma, lw=1.0, color = 'black', label = 'Backbone')
	ax.plot(-BBgamma, -BBsigma, lw=1.0, color = 'black')
	ax.plot(gamma,   sigma,   lw=1.0, color = 'red',   label = 'Uniaxial loading')

	ax.set_xlabel('Strain [%]', fontsize=18)
	ax.set_ylabel('Stress [kPa]', fontsize=18)


	fig.savefig('plot.png', dpi=300)





### PROGRAM ###


# Reading acceleration file
accelfile   = '../EXAMPLES/P1/outputfiles/accelx001'
time, accel = readfile(accelfile, 1.0, 1.0)

# Filtering the acceleration below 10 Hz by Butterworth filter
dt   		= time[2]-time[1]
accel    	= lowpass(accel, 10.0,  df=1/dt, corners=2, zerophase=True)


# Reading stress-strain file
sigepsfile   = '../EXAMPLES/P1/outputfiles/StressStrainxz005'
gamma, sigma = readfile(sigepsfile, 1e2, 1e-3)


# Reading backbone curve
bbfile      = '../EXAMPLES/P1/outputfiles/SOILHYPER'
BBgamma, BBsigma = readfile(bbfile, 1e2, 1e-3)



# Plotting all 
plot_them(time,accel,gamma,sigma,BBgamma,BBsigma)

#