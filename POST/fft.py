import matplotlib.pyplot as plt
import pylab as p
import numpy as np
import seaborn as sns
from filter import *
from   math   import *



# FFT computation

def fourier(x,dtoutput):
   # Tapering the data

   # FFT
   print 'dt', dtoutput
   s = np.abs( dtoutput * np.fft.fft(x,n=4096) )
   f = faxis(s,dtoutput)
   print 'df', f[1]-f[0]

   return f[:2049], s[:2049]
#

def faxis(x,dt):
   n  = len(x)
   df = 1.0 / (n*dt)
   f  = []
   for i in range( n ):
      f.append( df*i )

   return f
  #

# Konno-Ohmachi smoothening

def ko(y,dx,bexp):
  nx      = len(y)
  fratio  = 10.0**(2.5/bexp)
  ylis    = range( nx )
  ylis[0] = y[0]

  for ix in range( 1,nx ):
     fc  = float(ix)*dx
     fc1 = fc/fratio
     fc2 = fc*fratio
     ix1 = int(fc1/dx)
     ix2 = int(fc2/dx) + 1
     if ix1 <= 0:  ix1 = 1
     if ix2 >= nx: ix2 = nx
     a1 = 0.0
     a2 = 0.0
     for j in range( ix1,ix2 ):
        if j != ix:
           c1 = bexp*log10(float(j)*dx/fc)
           c1 = (sin(c1)/c1)**4
           a2 = a2+c1
           a1 = a1+c1*y[j]
        else:
           a2 = a2+1.0
           a1 = a1+y[ix]
     ylis[ix] = a1 / a2

  for ix in range( nx ):
     y[ix] = ylis[ix]

  return y
#


def readfile(dosya):

	print
	print 'Reading file ', dosya

	x = np.genfromtxt(dosya, usecols=0)
	y = np.genfromtxt(dosya, usecols=1)

	return x, y
#


def filtrele (dt,array,alt,ust):

	array = lowpass  (array, ust, 1/dt, corners=2, zerophase=True)
	array = highpass (array, alt, 1/dt, corners=2, zerophase=True)

	return array
#



def plotelo (time,E,N,Z, freq, EFFT, NFFT, ZFFT, kaydet,savename):

  sns.set_style('whitegrid')
  fig = plt.figure(figsize=(11.69,8.27)) 
  plt.subplots_adjust(hspace=0.4, right=0.98, left=0.1, wspace=0.4)

  # Vitesse en temps - Composante EW
  ax = fig.add_subplot(321)
  ax.set_xlim([0,50])
  # ax.set_ylim([-0.06,0.06])


  ax.plot(time,E,  color='black', lw=1.0)

  ax.set_ylabel('Acceleration [m/s/s]', fontsize=18)
  ax.yaxis.labelpad = 24
  plt.setp(ax.get_xticklabels(), visible=False)
  # plt.title('Fault - Parallel component'+'\n', fontsize=20)



  # Vitesse en temps - Composante NS
  ax = fig.add_subplot(323)
  ax.set_xlim([0,50])
  # ax.set_ylim([-0.06,0.06])

  ax.plot(time,N,  color='black', lw=1.0)

  ax.set_ylabel('Acceleration [m/s/s]', fontsize=18)
  ax.yaxis.labelpad = 24
  plt.setp(ax.get_xticklabels(), visible=False)
  # plt.title('Fault - Normal component'+'\n', fontsize=20)



  # Vitesse en temps - Composante verticale
  ax = fig.add_subplot(325)
  ax.set_xlim([0,50])
  # ax.set_ylim([-0.06,0.06])

  ax.plot(time,Z,  color='black', lw=1.0)

  ax.set_ylabel('Acceleration [m/s/s]', fontsize=18)
  ax.set_xlabel('Time [s]', fontsize=18)

  ax.yaxis.labelpad = 24
  plt.setp(ax.get_xticklabels(), visible=True)
  # plt.title('Vertical component'+'\n', fontsize=20)



  # Transforme de Fourier vs frequence - Composante EW
  ax = fig.add_subplot(322)
  ax.set_xlim([0.1,10])
  # ax.set_ylim([1e-4,1e0])

  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.plot(freq,EFFT,  color='black', lw=1.0, label='EW')

  ax.set_ylabel('FFT [m/s]', fontsize=18)
  ax.yaxis.labelpad = 24
  plt.setp(ax.get_xticklabels(), visible=False)
  # plt.title('Fault-Parallel component'+'\n', fontsize=20)

  # Transforme de Fourier vs frequence - Composante NS
  ax = fig.add_subplot(324)
  ax.set_xlim([0.1,10])
  # ax.set_ylim([1e-4,1e0])

  ax.set_xscale('log')
  ax.set_yscale('log')  
  ax.plot(freq,NFFT,  color='black', lw=1.0, label='NS')

  ax.set_ylabel('FFT [m/s]', fontsize=18)
  ax.yaxis.labelpad = 24
  plt.setp(ax.get_xticklabels(), visible=False)
  # plt.title('Fault-Normal component'+'\n', fontsize=20)

  # Transforme de Fourier vs frequence - Composante verticale
  ax = fig.add_subplot(326)
  ax.set_xlim([0.1,10])
  # ax.set_ylim([1e-4,1e0])

  ax.set_xscale('log')
  ax.set_yscale('log')  
  ax.plot(freq,ZFFT,  color='black', lw=1.0, label='UD')

  ax.set_ylabel('FFT [m/s]', fontsize=18)
  ax.set_xlabel('Frequency [Hz]', fontsize=18)
  ax.yaxis.labelpad = 24
  plt.setp(ax.get_xticklabels(), visible=True)
  # plt.title('Vertical component'+'\n', fontsize=20)

  ax.legend()


  if kaydet:
    fig.savefig(savename,dpi=300)
  else:
    plt.show()
#






### PROGRAM ###

savename = 'WRLA_acceleration_fft.png'

time , accx = readfile('../EXAMPLES/WILDLIFE/EFFECTIVE_STRESS_ANALYSIS/outputfiles/accelx001')
time , accy = readfile('../EXAMPLES/WILDLIFE/EFFECTIVE_STRESS_ANALYSIS/outputfiles/accely001')
time , accz = readfile('../EXAMPLES/WILDLIFE/EFFECTIVE_STRESS_ANALYSIS/outputfiles/accelz001')


accx = filtrele(time[2]-time[1],accx,0.1,10.0)
accy = filtrele(time[2]-time[1],accy,0.1,10.0)
accz = filtrele(time[2]-time[1],accz,0.1,10.0)


freq, FFT_x = fourier(accx, time[2]-time[1])
freq, FFT_y = fourier(accy, time[2]-time[1])
freq, FFT_z = fourier(accz, time[2]-time[1])


FFT_x  = ko(FFT_x,(freq[1]-freq[0]),40.0)
FFT_y  = ko(FFT_y,(freq[1]-freq[0]),40.0)
FFT_Z  = ko(FFT_z,(freq[1]-freq[0]),40.0)



plotelo (time,accx,accy,accz, freq, FFT_x, FFT_y, FFT_z, True, savename)

#