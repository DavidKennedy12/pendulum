import numpy as np
import matplotlib.pyplot as plt
from chain import *
from scipy.integrate import RK45 as rk

def coordinates(theta):
        x, y = np.cumsum(np.sin(theta), axis = 1), -np.cumsum(np.cos(theta), axis = 1)
        return x, y
    
def totalEnergy(self):
        return self.getKineticEnergy() + self.getPotentialEnergy()
    
def potentialEnergy(theta):
        y = -np.cumsum(np.cos(theta), axis = 1) + len(theta[0])
        return np.sum(y, axis = 1)
    
def kineticEnergy(theta, omega):
        xVel = np.cumsum(omega * np.cos(theta), axis = 1)
        yVel = np.cumsum(omega * np.sin(theta), axis = 1)
        kinetic = (1/2) * (xVel**2 + yVel**2)
        return np.sum(kinetic, axis = 1)    

def loadData(file):
    theta = np.load("data/" + file + "_theta.npy")
    omega = np.load("data/" + file + "_omega.npy")
    tension = np.load("data/" + file + "_tension.npy")
    x,y = coordinates(theta)
    potential = potentialEnergy(theta)
    kinetic = kineticEnergy(theta, omega)
    return theta, omega, tension, x, y, potential, kinetic

def plot(x, y):
    N = len(x)
    fig, ax = plt.subplots(figsize=(10,10))
    plt.axis('off')
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    ax.set_xlim([-N-1,N+1])
    ax.set_ylim([-N-1,N+1])
    ax.scatter(0,0, c = "red", s=50)
    ax.scatter(x,y, c = "blue", s = 50 )  