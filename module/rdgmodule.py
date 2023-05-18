#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Radon in Different Gases Module
 
Includes functions and calculation for the RDG simulation

Created on Wed Jun 30 13:17:21 2021

AUTHOR: ROBERT RENZ MARCELO GREGORIO
EMAIL: robert.gregorio@sheffield.ac.uk
"""

# =============================================================================
# MODULE DEPENDACIES
# =============================================================================

#RANDOM PACKAGES
from numpy.random import uniform
from numpy.random import random
from numpy.random import normal

#MATHS PAKAGES
from numpy import sqrt
from numpy import array
from numpy import pi
from numpy import cos
from numpy import sin
from numpy import log



#DATA ANALYSIS 
import pandas as pd
from scipy import spatial
from shapely.geometry import LineString

# PLOTTING
import matplotlib.pyplot as plt

# =============================================================================
# RANDOM RADON POSITION IN DETECTOR
# =============================================================================

def RandomCoordinates():    
    # Effectively y is h and r is x-z plane
    # Boundaries
    h_max = 6.1
    h_min = -5.1
    r_max = 5.1
    
    while True:
        # Provides a random h coordinate in detector volume including detector area
        h_ran=uniform(h_min,h_max)
        
         # Semicircle correction - If h is in the semicircle correct for random probability distribution 
        if h_ran < 0 : 
        
            # Weighted probability as a function of height eg smaller volume less likely for random generation
            prob_h_ran=sqrt(r_max**2-h_ran**2)/r_max            
            if random() > prob_h_ran: # If weighted probability does not allow start while loop again
                continue
            
            # Generate r value within spherical limits
            r_hmax=sqrt(r_max**2-h_ran**2)  # maxium radius for a given random height
            r_ran =r_hmax * sqrt(random()) # random r value in a circle (https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly/5838055
          
    
        #  h is in not in the semicircle equal probability no need for correction
        else: 
            r_ran=r_max*sqrt(random()) # Generate r value within cylindrical limit
        
          
        # If random coordinates are in PIPS detector repeat until new random position is made
        if (6.6<= h_ran <= 7.6 and -1.6 <= r_ran <=1.6) or (7.6<= h_ran <=8.9 and -0.5 <= r_ran <=0.5):
            print('WARNING: Random position is in PIPS detector')
            continue 
        
        # Allow for negative and positive r solutions 
        if random()<0.5:
           r_ran = -r_ran


        break
        
    return h_ran,r_ran

# =============================================================================
# CORRECTION FUNCTIONS
# =============================================================================
# Mean Square Displacement Correction
def msdCorrection(P, D, dt):
    msdsol=sqrt(3*2*D*dt) #spread from Mean square displacement in 3D cm/s 
    
    dhmsd=normal(loc=P[0],scale=msdsol,size=1) #Normal distribution about h - new P[0]
    drmsd=normal(loc=P[1],scale=msdsol,size=1) #Normal distribution about y  -new P[1]
            
    P = array([dhmsd,drmsd]) # position mean square disaplacement correction 

    return P


# Neutral Charged Species Correction
def NeutralSpeciesCorrection():
    if random() >=0.880:  #Probability of positive ion is 88% funtion is false 88% of the time
       ENDPOINT='NEUTRAL SPECIES'
       
    else:
        ENDPOINT = 'SUSPENDED'
       
    return ENDPOINT

# PIPS Correction
def PIPSdetection():
    if random() <0.5:  # isotropic alpha emission only half of the collected atoms on the PIPS can be detected
       ENDPOINT='DETECTED'
    else:
        ENDPOINT = 'COLLECTED'
    return ENDPOINT

# Neural Networks Correction Model

# Humidity Correction


# =============================================================================
# CARRIER GAS PARAMETER CALCULATIONS
# =============================================================================

# Constants
k_b= 1.380649e-23 #Boltzman constant in JK-1
amu=1.66054e-27 #1au to kg
ec=1.60217662e-19 #electron charge
N_A=6.02214086e23 #Avogadros Constant mol-1
R= 8.31446261815324 #GAS CONSTANT J⋅K−1⋅mol−1.
q=1*ec #ion charge C


# POLONIUM Properties 
m_po=218.008966 *amu #polonium mass in kg

# POLONIUM KINETIC RADIUS (NOT DIAMETER!) CALCULATED FROM MOBILITIES
#r_po=2.7877189569031766e-10 #Busgia 1981
r_po=3.143580223581969e-10 #Chu Hopke 1988

def CarrierGasParameters(M_gas,r_gas, T, P):
    #Ideal Gas Density for given T and P [kg/m3]
    rho_gas=(M_gas*1e-3*P)/(R*T) # COMPARE WITH DENSITIES: https://www.engineeringtoolbox.com/gas-density-d_158.html
    #Collision cross section
    sigma=pi*(r_po+r_gas)**2 #Collision cross sectional area of Po-ion with carrier gas 

    #Mobility in cm2/(V*S)
    mu=(q*M_gas*1e-3)/(2*rho_gas*N_A*sigma)*sqrt(3/(m_po*k_b*T))*1e4     
    
    #Diffusion coefficient cm2
    D=(1e-3*M_gas)/(2*sigma*rho_gas*N_A)*sqrt((3*k_b*T)/m_po)*1e4
    return rho_gas, mu, D






# =============================================================================
# EVENTS
# =============================================================================

# EVENT: Rn-222 alpha decay (5.49 MeV)
def Rn222Decay (P,SRIM_range, SRIM_straggle):
    # Random 3D direction - Spherical coordinates (recoil, theta, pi)
    recoil=float(normal(loc=SRIM_range,scale=SRIM_straggle,size=1))
    theta = uniform(-pi/2,pi/2)
    phi = uniform(0,2*pi)
    
    # Convert the recoil to r and h components (http://nbeloglazov.com/2017/04/09/random-vector-generation.html)
    h_recoil = recoil * sin(theta)
    r_recoil = recoil * cos(theta) * cos(phi)

    # Convert recoil vector to P(h,r)
    recoil = array([[h_recoil],[r_recoil]])
    # Update recoil position
    P = P + recoil
    return P

# Suspended Nuclei Event

# =============================================================================
# Smalll ion recombination 
# =============================================================================
def SIR(W, C_rn, mu, dt):
    e=1.602e-19 # electron charge
    E_rn= 5.40e6 # Alpha decay energy of Radon 222 in eV
    eps  = 8.85e-12 #permittivity of free space F/m
    mu = mu*1e-4
    
    p =  (e*E_rn*C_rn*mu)/(eps*W)*dt
    
    
    if random() <=p:  #Probability of SIR for carrier gas and radon concetration  - funtion is true p times of the time
        ENDPOINT = 'SIR'
    else:
        ENDPOINT = 'SUSPENDED'
    
    return ENDPOINT







# Surface Nuclei Event 

# =============================================================================
# NUCLEI TRACKING 
# =============================================================================
# Determine Nearest E-field Node 
path_to_merge = 'NODE_SOL_MERGE_RAD7.csv'

def NearestNodeEfield(h ,r):
    # LOAD ANSYS SOLUTION
    df=pd.read_csv(path_to_merge)
    
   
    df_mesh=df[['X','Y']]  #  NODAL COORDINATES
    A=df_mesh.values # Converting pandas to numpy array for scipy
    df_efsol=df[['NODE','EFX','EFY','EFSUM']]   #  NODAL ELECTRICFIELD SOLUTIONS
    
    
    #FIND NEAREST NODE NUMBER FOR GIVEN h and r
    position=[float(h),float(r)] #Position to find
    nearest=A[spatial.KDTree(A).query(position)[1]]
    distance,index = spatial.KDTree(A).query(position)
    NNODE=df.loc[(df['X']==nearest[0])&(df['Y']==nearest[1])].iloc[0]['NODE']
       
    #FIND CORRESPONDING E-FIELD FROM NODE NUMBER
    NEF_x=df_efsol.loc[df_efsol['NODE']==NNODE].iloc[0]['EFX']
    NEF_y=df_efsol.loc[df_efsol['NODE']==NNODE].iloc[0]['EFY']
    NEF=array([[NEF_x],[NEF_y]])
    
    return NEF



# Nuclei tracking boundaries walls, detector hit or detector case
def Boundaries(Pprev, P):
    
    # DETECTOR BOUNDARIES 
    # Condition for detector window count - path intersects window line and change in dh > 0
    dwindow = LineString([(3.9,-0.973),(3.9,0.973)])   #Define 1D Detector Window
    path = LineString([(Pprev[0],Pprev[1]), (P[0],P[1])]) # Create line path - Previous position P(h,r) or Ppre(h,r)
    count=dwindow.intersects(path) # Count=True if intersect between linepath and detector window
    dh = P[0] - Pprev[0] #change in dh
    
    
    # COLLECTED - WINDOW HIT
    if count==True and dh > 0: #if it passes through the detector AND dx is positve its a count
        ENDPOINT='COLLECTED'

    
    # FALSE COLLECTED
    elif count==True and dh < 0: #if it passes through the detector AND dx is negative therefore moving backwards its not a count
        ENDPOINT='DCASE'
        print('FALSECOLLECTED')
        
    # CASE - ELSE
    elif 3.9<= P[0] <=6.1 and -1.4 <= P[1] <=1.4 and count==False : #Inside the detector square and does not pass through detector window
        ENDPOINT='DCASE'

    
    # WALL BOUNDARIES P(h , r)
     # |r| > 5.1 - WALL SIDES 
     #  h > 8.9 - HV BOARD
     #  h < 5.1 - DOME
     
    elif abs(P[1]) > 5.1 or  P[0] > 8.9 or P[0] < - 5.1:
        ENDPOINT = 'WALLS'
    
    # h is inside semicircle and |r| > r max for given h - WALL DOME
    elif P[0] < 0 and P[1] > sqrt(5.1**2-P[0]**2):
         ENDPOINT = 'WALLS'
    
    
    else:
         ENDPOINT='SUSPENDED'

    
    return(ENDPOINT)




# Plotting
def TrackPlot(ENDPOINT, singledata):
    
    if ENDPOINT == 'WALLS' or 'DCASE':
        plt.plot(singledata['h (cm)'],singledata['r (cm)'], linestyle='-.', marker='o',markersize=0.5, markerfacecolor='white', color='#D01E2F',alpha=0.9, linewidth=0.5)

    
    if ENDPOINT == 'SIR':
        plt.plot(singledata['h (cm)'],singledata['r (cm)'], linestyle='-.', marker='o',markersize=0.5, markerfacecolor='white', color='yellow',alpha=0.9, linewidth=0.5)



    if ENDPOINT == 'COLLECTED':
        plt.plot(singledata['h (cm)'],singledata['r (cm)'], linestyle='-.', marker='o',markersize=0.5, markerfacecolor='white', color='black',alpha=0.9, linewidth=0.5)
        
    if ENDPOINT == 'NEUTRAL SPECIES':
         plt.plot(singledata['h (cm)'],singledata['r (cm)'], linestyle='-.',markersize=2, marker="o", markerfacecolor='orange', color='orange',alpha=0.9)
         
    return
        
        
        

# =============================================================================
# #   PLOTTING RAD7 MODEL 
# # =============================================================================



# =============================================================================
# OUTPUT CALCULATION FUNCTIONS
# =============================================================================

# Convert radon concentration in RAD7 to nuclei 
def RAD7CONCtoNUCLEI(ConcBq):
    #A = N*lambda; C = A/V
    thalf_Rn222 = 3.3035e5 #3.8235 days 
    dconst = log(2)/thalf_Rn222
    #RAD7 Volume (m3)
    V = 1e-3
    #Number of atoms for RAD7 concentration 
    N = int(V*ConcBq/dconst)
    
    return N

# Convert nuclei to radon concentration in RAD7
def NUCLEItoRAD7CONC(N):
    #A = N*lambda; C = A/V
    thalf_Rn222 = 3.3035e5 #3.8235 days 
    dconst = log(2)/thalf_Rn222
    #RAD7 Volume (m3)
    V = 1e-3
    #RAD7 concentration for number of nuclei 
    ConcBq = N*dconst/V
    
    return ConcBq


