import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import pylab as p
import numpy as np
import sys
import os 



# Creating datamesh file 
# name  = 'datamesh'
# title = 'WILDLIFE REFUGE LIQUEFACTION ARRAY borehole model: New mesh'


name  = raw_input("Enter the Mesh file name to be created  ")
title = raw_input("Enter the Model name for information    ")


print "*************************************************************************************"
print "Specify total element number please :) "
giri = raw_input("Integer only ")
nel  = int(giri) 
print "***"
print "Specify total material/domain number please (PML included): "
giri = raw_input("Integer only ")
ndom = int(giri)



ndiric = 0
nnode = nel+ 1
semfile = open(name, "w")

semfile.write("%s \n"  % (' #'))
semfile.write("%s \n"  % (' # Mesh file for '+str(title)))
semfile.write("%s \n"  % (' #     1D Wave Propagation '))
semfile.write("%s \n"  % (' #'))
semfile.write("%s \n"  % (' # Number of elements, Number of nodes, Number of domains'))
semfile.write(" \t %s \t %s \t %s  \n"  % ( str(nel), str(nnode), str(ndom)))
semfile.write("%s \n"  % (' # Nodes number, coordinates and type (D:Dirichlet,N:Neumann,I:Internal,B:Borehole)'))



coord = 0.0
node  = 1
print "*************************************************************************************"
print "Specify upper boundary condition please :) "
giri = raw_input("N for Neumann/ D for Dirichlet: ")
if (giri == 'N'):
	print "You want Neumann condition at top, OK!"
elif (giri == 'D'):
	print "You want Dirichlet condition at top, OK!"
	ndiric = ndiric+ 1
else:
	print "only N or D accepted...!"
	sys.exit()

print "***"
bc = str(giri)
semfile.write("\t %s \t %2f \t %s \n"  % (str(node), coord, bc ))
bc = 'I'

print "*************************************************************************************"
print "Setting element size types"
elems = []
while True:
    entered = raw_input("Please enter an element size or leave a blank line to quit: ")
    if not entered: break
    print "You entered ", entered
    elems.append(float(entered))
elems.sort()
print "*************************************************************************************"
print "You defined element types as: ", elems


domelem = {}
mindom  = {}
maxdom  = {}
print "***"
print "Let's define the element type numbers for each domain..."
for i in np.arange(ndom):
	domelem[i] = 0
	mindom[i] = max(elems)
	maxdom[i] = min(elems)

	print "*********************************************************************************"

	for j in np.arange(len(elems)):
		print "Domain ", i+1, " :"
		print "Element number for size :", elems[j]
		number = raw_input("Specify please: ")
		if not number:
			pass
		else:
			print "You entered ", number
			print "***"

			if (int(number) > 0):
				for k in np.arange(int(number)):
					node  =  node+ 1
					coord = coord+ float(elems[j])

					if (elems[j]< mindom[i]):
						mindom[i] = elems[j]

					if (elems[j]> maxdom[i]):
						maxdom[i] = elems[j]

					if (node < nnode):
						semfile.write("\t %s \t %2f \t %s \n"  % (str(node), coord, bc ))	





				domelem[i]	= domelem[i]+ int(number)


	print "*********************************************************************************"
	print "Domain ", i+1
	print "Element number: ", domelem[i]
	print "Min size: ", mindom[i]
	print "Max size: ", maxdom[i]


	if (node == nnode):
		print "Specify lower boundary condition please :) "
		giri = raw_input("N for Neumann/ D for Dirichlet/ B for Borehole: ")
		if (giri == 'N'):
			print "You want Neumann condition at bottom, OK!"
			bc = 'N'
		elif (giri == 'D'):
			print "You want Dirichlet condition at bottom, OK!"
			ndiric = ndiric+ 1
			bc = 'D'
		elif (giri == 'B'):
			print "You want Borehole condition at bottom, OK!"
			bc = 'B'
		else:
			print "only N or D or B accepted...!"
			sys.exit()

		semfile.write("\t %s \t %2f \t %s \n"  % (str(node), coord, bc ))






semfile.write("%s \n"  % (' # Elements number, corresponding domain and nodes number'))
elsayi = 0
for i in np.arange(ndom):
	for j in np.arange(int(domelem[i])):
		elsayi = elsayi+ 1
		semfile.write("\t %s \t %s \t %s \t %s \n"  % (str(elsayi), str(i+1), str(elsayi), str(elsayi+1) ))		

semfile.write("%s \n"  % (' # Domains number, corresponding material number, Number of elements, dxmin, dxmax'))
for i in np.arange(ndom):
	semfile.write("\t %s \t %s \t %s \t %1f \t %1f \n"  % (str(i+1), str(i+1), str(domelem[i]), float(mindom[i]), float(maxdom[i])))

semfile.write("%s \n"  % (' # Parameters for the Dirichlet conditions Left/Right (vmax,t0,tmax)'))
for i in np.arange(ndiric):
	if (ndiric==2 and i==0):
		semfile.write(" %s \t %1f \t %1f \t %1f \n"  % ( 'L', 0.0, 0.0, 0.0  ))
	else:
		semfile.write(" %s \t %1f \t %1f \t %1f \n"  % ( 'R', 0.0, 0.0, 0.0  ))

print "*************************************************************************************"
print "File ", name, " is ready!"