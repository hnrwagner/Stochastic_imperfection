
"""
Created on 06.04.2021
###############################################################################
###############################################################################
@author: H.N.Ronald Wagner 
###############################################################################

###############################################################################
"""

import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from numpy import linalg as LA
from numpy.linalg import matrix_rank

###############################################################################
    
# Cylinder Geometry Identifier

myShellName = "IW1"

# File Names for the analysis 

myFileName = str(myShellName)+"_Nodes.txt"
myImageName = str(myShellName)+"_Shell.png"

# Length of the cylinder

myLength = 100.0

# Radius of the cylinder

myRadius = 33.0

# Wall thickness of the cylinder

myThickness = 0.01

# Number of Nodes in axial and circumferential direction

na = 61

nc = 240


# Displacement field

w = np.zeros((na,nc))


xyz = np.zeros((na*nc,3))

    
"""
###############################################################################
@equation
The following equations are described in
"Stochastic Modeling of Geometric Imperfections in Aboveground Storage Tanks for Probabilistic Buckling Capacity Estimation"
ASCE-ASME J. Risk Uncertainty Eng. Syst., Part A: Civ. Eng., C4015005
###############################################################################
###############################################################################
"""

    
#------------------------------------------------------------------------------    
    
# for Cylinder

#------------------------------------------------------------------------------
#number of panels in meridian and circumferential direction
P_m = 2
P_c = 10
n1 = int(1*P_m)
n2 = int(1*P_c)
# Fourier coefficient    
A = np.zeros((n1,n2))
B = np.zeros((n1,n2))
# normally distributed random variable - mean = 0 & stdev = 1  
mu = 0
sigma = 1
# non-negative constants
a = 1
b = 1
#decay coefficient
alpha = 0.1
beta = 0.1

for k in range(0,n1,1):
    for l in range(0,n2,1):

        A[k][l] = np.abs(np.random.normal(mu, sigma, 1)) * np.exp(-k*alpha-l*beta) * (a + abs(np.cos((k)*np.pi/P_m))) * (b + abs(np.cos((l)*np.pi/P_c)))
        B[k][l] = np.abs(np.random.normal(mu, sigma, 1)) * np.exp(-k*alpha-l*beta) * (a + abs(np.sin((k)*np.pi/P_m))) * (b + abs(np.sin((l)*np.pi/P_c)))
#------------------------------------------------------------------------------  
  

for i in range(0,na,1):
    for j in range(0,nc,1):
        x = myLength*(i)/(na)
        y = 2.0*np.pi*myRadius*(j)/(nc)
        
        for k in range(0,n1,1):
            for l in range(0,n2,1):
                # half wave cosine approach
                w[i][j] = w[i][j] + A[k][l]*np.cos((k)*np.pi*x/myLength)*np.cos((l)*y/myRadius-B[k][l])
 
       
        w[i][j] = w[i][j] * myThickness
        xyz[(i)*nc+j][0] = (myRadius-w[i][j])*np.cos(y/myRadius)
        xyz[(i)*nc+j][1] = (myRadius-w[i][j])*np.sin(y/myRadius)
        xyz[(i)*nc+j][2] = x*1.016666671
##    
##################################################################################
##   
        
myFileName 
myImageName       
nm = len(xyz)

ABAQUS_NODES = np.zeros((nm,4))
#
for i in range(0,nm,1):

    ABAQUS_NODES[i][0] = i+1
    ABAQUS_NODES[i][1] = xyz[i][0]
    ABAQUS_NODES[i][2] = xyz[i][1]
    ABAQUS_NODES[i][3] = xyz[i][2]
    

np.savetxt(myFileName,ABAQUS_NODES,fmt ='%i, %10.5f, %10.5f, %10.5f')

print("\n")
print("###############################################################################")
print("End of Calculation Number 1")
print("###############################################################################")

cmap = cm.get_cmap('jet')
p = plt.pcolor(w, cmap = cmap)
cb = plt.colorbar(p, orientation = 'vertical')
plt.xlabel("Circumferential Nodes",fontsize=14)
plt.ylabel("Axial Nodes",fontsize=14)
#plt.xticks(np.arange(0,405,45))
cb.set_label('Radial Displacement [mm]',fontweight='bold',fontsize=14)
plt.savefig(myImageName, dpi = 1000)
# cb.remove()
# p.remove()



