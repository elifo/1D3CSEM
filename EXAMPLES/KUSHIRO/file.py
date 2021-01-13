import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import pylab as p
import numpy as np




time    = np.genfromtxt("kp_xyz", usecols=0)
velocx  = np.genfromtxt("kp_xyz", usecols=1)
velocz  = np.genfromtxt("kp_xyz", usecols=3)

savefile = '/Users/elifo/Work/THESE/CODES/2D/EXAMPLES/SITE_EFFECTS/MANUAL/EXAMPLES/KUSHIRO_MODEL/kp_xz'


# 1D
# velocy = np.zeros(len(velocx))
# print "Creating file: ", savefile 
# output = np.column_stack((time.flatten(),velocx.flatten(),velocy.flatten(),velocz.flatten()))
# np.savetxt(savefile, output,delimiter='   ')
# print "Done!"


# 2D
print "Creating file: ", savefile 
output = np.column_stack((time.flatten(),velocx.flatten(),velocz.flatten()))
np.savetxt(savefile, output,delimiter='   ')
print "Done!"