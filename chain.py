import numpy as np
import scipy as sp
class Chain: 
    def __init__(self, n):
        self.time = 0
        self.n = n
        self.theta = np.zeros(n) #Angles between successive masses
        self.omega = np.zeros(n) #angular velocities
        self.alpha = np.zeros(n) #delta angle
        self.tension = np.zeros(n) #tension between masses
        self.delta = np.zeros(n-1) #Difference in angle between adjacent masses

        self.L = np.zeros((3,n))       #The L matrix used to compute Tension
        self.L[1,1:] = np.full(n-1,2)
        self.L[1,0] = 1
        self.e0 = np.zeros(self.n)
        self.e0[0] = 1
        
        self.L = np.zeros((3,n))
        self.L[1,1:] = np.full(n-1,2)
        self.L[1,0] = 1
    
    def __call__(self, t, z): #z is [theta, omega]
        
        time = t
        self.theta, self.omega = z[:self.n], z[self.n:]
        self.delta = np.diff(self.theta)  
        
        self.L[0,1:] = -np.cos(self.delta)
        self.L[2, 0:-1] = -np.cos(self.delta) 
        
        theta1 = np.append(self.delta, 0)
        theta2 = np.append(0, self.delta)
        lDiagonal = 2*np.ones(self.n)
        lDiagonal[0] = 1   
        D = sp.sparse.dia_matrix(([-np.sin(theta1),np.sin(theta2)],[-1,+1]),shape=(self.n,self.n))  
        
        RHS = self.omega**2
        RHS[0] += np.cos(self.theta[0]);
        T = sp.linalg.solve_banded((1,1), self.L, RHS)
        self.tension = T
        self.alpha = D.dot(T)
        self.alpha[0] -= np.sin(self.theta[0])
        
        return np.concatenate([self.omega, self.alpha])

    def setState(self, z):
        self.theta, self.omega = z[:self.n], z[self.n:]
    
    def createState(self, a, b):
        self.theta, self.omega = np.pi - np.arctan2(a, np.linspace(-b, 1-b, self.n)), np.zeros(self.n)
        return np.concatenate([self.theta, self.omega])
    
    def getCoordinates(self):
        x = np.cumsum(np.sin(self.theta))
        y = -np.cumsum(np.cos(self.theta))
        return x, y
       
    def getPotentialEnergy(self):
        y = -np.cumsum(np.cos(self.theta)) + self.n
        return np.sum(y)
    
    def getKineticEnergy(self):
        xVel = np.cumsum(self.omega * np.cos(self.theta))
        yVel = np.cumsum(self.omega * np.sin(self.theta))
        kinetic = (1/2) * (xVel**2 + yVel**2)
        return np.sum(kinetic)
    
    def getEnergy(self):
        return self.getKineticEnergy() + self.getPotentialEnergy()