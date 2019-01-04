import numpy as np
import scipy as sp
from chain import *
from scipy.integrate import RK45 as rk

def createData(N, z, name):
	f = Chain(N)
	f.setState(z)
	stepper = rk(f,0,z,1000, 0.002)
	samples = 1001
	theta = np.zeros((samples,N))
	omega = np.zeros((samples, N))
	tension = np.zeros((samples, N))

	i = 0
	j = 0
	while stepper.t < 1000:
		if (i % 500)  == 0:
			theta[j] = f.theta
			omega[j] = f.omega
			tension[j] = f.tension
			j = j+1      
			print(j) 
		stepper.step()
		i = i + 1   
		
	np.save(name + "_theta.npy", theta)
	np.save(name + "_omega.npy", omega)
	np.save(name + "_tension.npy", tension)
	print(name + " done")


print("Starting to generate chains")

#Chain 1a
N = 10
angles = np.zeros(N)
omegas = np.zeros(N)
angles[0:N//2:1] = np.pi/4
angles[N//2::1] = 3*np.pi/4
z = np.concatenate([angles, omegas])
createData(N, z, "1a")

#Chain 1b
N = 100
angles = np.zeros(N)
omegas = np.zeros(N)
angles[0:N//2:1] = np.pi/4
angles[N//2::1] = 3*np.pi/4
z = np.concatenate([angles, omegas])
createData(N, z, "1b")

#Chain 1c
N = 1000
angles = np.zeros(N)
omegas = np.zeros(N)
angles[0:N//2:1] = np.pi/4
angles[N//2::1] = 3*np.pi/4
z = np.concatenate([angles, omegas])
createData(N, z, "1c")

#Chain 2a
N = 5
angles = np.full(N,(3*np.pi)/4)
omegas = np.zeros(N)
z = np.concatenate([angles, omegas])
# createData(N, z, "2a")

#Chain 2b
N = 5
angles = np.full(N,(3*np.pi)/4 ) + 0.001
omegas = np.zeros(N)
z = np.concatenate([angles, omegas])
createData(N, z, "2b")

# #Chain 2c
N = 5
angles = np.full(N,(3*np.pi)/4) - 0.001
omegas = np.zeros(N)
z = np.concatenate([angles, omegas])
createData(N, z, "2c")

#Chain 3a
N = 1000
angles = np.full(N,(3*np.pi)/4)
omegas = np.zeros(N)
z = np.concatenate([angles, omegas])
createData(N, z, "3a")

#Chain 3b
N = 1000
angles = np.full(N,(3*np.pi)/4) + 0.001
omegas = np.zeros(N)
z = np.concatenate([angles, omegas])
createData(N, z, "3b")

#Chain 3c
N = 1000
angles = np.full(N,(3*np.pi)/4) - 0.001
omegas = np.zeros(N)
z = np.concatenate([angles, omegas])
createData(N, z, "3c")

#Chain 4a
N = 10
angles = np.zeros(N)
omegas = np.zeros(N)
angles[0::2] = np.pi/2
angles[1::2] = 3*np.pi/2
z = np.concatenate([angles, omegas])
createData(N, z, "4a")

#Chain 4b
N = 100
angles = np.zeros(N)
omegas = np.zeros(N)
angles[0::2] = np.pi/2
angles[1::2] = 3*np.pi/2
z = np.concatenate([angles, omegas])
createData(N, z, "4b")

#Chain 4c
N = 1000
angles = np.zeros(N)
omegas = np.zeros(N)
angles[0::2] = np.pi/2
angles[1::2] = 3*np.pi/2
z = np.concatenate([angles, omegas])
createData(N, z, "4c")