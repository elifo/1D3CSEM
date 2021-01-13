WL		# Title of simulation
50.0  	        # Duration of the simulation
datamesh 	# Mesh input file
material	# Material parameters file
.false.       	# Source 
source        	# Source parameters file
stations      	# Receiver parameters file
.true.       	# Imposed signal
wl_xyz          # Imposed signal file		 (if size bigger than 50000, recompile!) 
1.0 		# Coefficient for input                               
50000        	# Time interval to skip snapshots
.true.       	# Save traces 			 
250		# Output factor 
.true.        	# Run time loop or not
.true.	      	# Viscoelasticity on
.true.       	# Nonlinearity on
.true.		# Pressure-dependent model
3		# Water level node number 
