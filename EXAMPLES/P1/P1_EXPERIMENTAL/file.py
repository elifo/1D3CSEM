import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import pylab as p
import numpy as np




time    = np.genfromtxt("ricker1c", usecols=0)
veloc   = np.genfromtxt("ricker1c", usecols=1)


savefile = 'ricker101'

velocy = np.zeros(len(veloc))
print "Creating file: ", savefile 
output = np.column_stack((time.flatten(),veloc.flatten(),velocy.flatten(),veloc.flatten()))
np.savetxt(savefile, output,delimiter='   ')
print "Done!"