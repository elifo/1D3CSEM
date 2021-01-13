import matplotlib.pyplot as plt
import pylab as p
import numpy as np
import seaborn as sns
from filter import *
from   math   import *



def readfile(dosya,col):

	print
	print 'Reading file ', dosya

	y = np.genfromtxt(dosya, usecols=col)

	return y
#

def cree_temps(pore,dosya):

  time = np.genfromtxt(dosya, usecols=0)
  dt   = time[1]-time[0]

  temps = np.zeros(len(pore))
  for i in np.arange(len(pore)):
    temps[i] = i* dt

  return temps
#


def read_sigeps(dosya):

  eps = np.genfromtxt(dosya, usecols=0)*1e2
  sig = np.genfromtxt(dosya, usecols=1)/1e3

  return eps, sig
#


def plotall(time, pore, gamma, sigma, G,r,S, m1,m2):
  
  lfailure = np.zeros((15))
  lphase   = np.zeros((15))
  array    = np.zeros((15))

  # Failure and phase transformation lines
  for i in np.arange(0,1.5,0.1):
    lfailure[i] = i* m1
    lphase  [i] = i* m2
    array   [i] = i



  # Plotting
  fig = p.figure(num=1,figsize=(12,8),dpi=300,facecolor=None,edgecolor='k'); 
  fig.subplots_adjust(left=0.15,right=0.90,bottom=0.1,top=0.90,wspace=0.2,hspace=0.2)

  sns.set_style(style='whitegrid')

  # Pore pressure excess vs time
  ax = fig.add_subplot(221)
  ax.plot(time, pore, lw=1.0, color = 'black' )

  ax.set_xlabel('Time [s]', fontsize=18)
  ax.set_ylabel('Pore pressure excess [kPa]', fontsize=18)


  # Stress vs Strain
  ax = fig.add_subplot(222)
  ax.plot(gamma, sigma, lw=1.0, color = 'black' )

  # ax.plot(time, gamma, lw=1.0, color = 'black' )


  ax.set_xlabel('Strain [%]', fontsize=18)
  ax.set_ylabel('Stress [kPa]', fontsize=18)


  # G modulus vs Time
  ax = fig.add_subplot(223)
  ax.plot(time, G/1e6, lw=1.0, color = 'black' )

  ax.set_xlabel('Time [s]', fontsize=18)
  ax.set_ylabel('Shear modulus [MPa]', fontsize=18)


  # Deviatoric plan
  ax = fig.add_subplot(224)
  ax.plot(S, r,    lw=1.0, color = 'royalblue' )
  ax.plot(array, lfailure, color = 'black')
  ax.plot(array, lphase,   color = 'black', linestyle='--')

  ax.set_xlabel('Normalized effective stress', fontsize=18)
  ax.set_ylabel('Normalized deviatoric stress', fontsize=18)


  fig.savefig('effective_stress_analysis.png',dpi=300)
#




### PROGRAM ###

# At station GL-4.0 m

file_accel      = '../EXAMPLES/WILDLIFE/EFFECTIVE_STRESS_ANALYSIS/outputfiles/accelx001'
file_effectives = '../EXAMPLES/WILDLIFE/EFFECTIVE_STRESS_ANALYSIS/outputfiles/PressEffectiveParams007'
file_modules    = '../EXAMPLES/WILDLIFE/EFFECTIVE_STRESS_ANALYSIS/outputfiles/PressSoilParams007'
file_sigeps     = '../EXAMPLES/WILDLIFE/EFFECTIVE_STRESS_ANALYSIS/outputfiles/StressStrainxz007'

m1 = 0.5299
m2 = 0.4067


pore_excess = readfile  (file_effectives,6)
time        = cree_temps(pore_excess, file_accel)


strain, stress = read_sigeps(file_sigeps)

Gmodule = readfile(file_modules,5)

r       = readfile(file_effectives,1)
S       = readfile(file_effectives,5)

plotall(time,pore_excess/1e3, strain,stress,Gmodule,r,S,m1,m2)


#