#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MONTE CARLO SIMULATION OF ELECTROSTATIC RADON DETECTOR IN DIFFERENT GASES

Created on Wed Jun 30 11:18:24 2021

AUTHOR: ROBERT RENZ MARCELO GREGORIO
EMAIL: robert.gregorio@sheffield.ac.uk

"""




# =============================================================================
# SIMULATION INPUT - CARRIER GAS
# =============================================================================

# SRIM PARAMETERS 
# http://www.srim.org/ NOTE: IN SRIM SIMULATION MAKE SURE USE SPECIFIC DENSITY FOR P AND T
# SRIM_range - mean projected range (um x e-4)cm 
# SRIM_straggle - standard deviation in range cm
# SRIM_range = 7.3e-4
# SRIM_straggle = 1.1e-4


# ENVIROMENTAL PARAMETERS
# T - Temperature (K)
# P - Pressure (Pa)
# RH - Relative Humidity (%)
T = 293
Pressure = 101325
RH = 10

# NEURAL NETWORKS INPUT
# EA - Affinity Energy (eV) https://webbook.nist.gov/chemistry/name-ser/ (Gas phase ion energetics data)
# EI - Ionisation Energy (eV) https://webbook.nist.gov/chemistry/name-ser/ (Gas phase ion energetics data)
EA = 1
IE = 1

#CARRIER GAS INFORMATION
# r_gas - Kinetic Diameter: https://en.wikipedia.org/wiki/Kinetic_diameter
# M_gas - Molar Mass :https://webbook.nist.gov/chemistry/name-ser/
# rho_gas - Density: Calculated in IDEALGAS_DENSITY_V1.py  (Compare with https://www.engineeringtoolbox.com/gas-density-d_158.html)

# # N2
# GAS='N2'
# r_gas=364e-12/2 # Kinetic radii (NOT DIAMETER!) of carrier gas in (diameter pm/2 e-12) m
# M_gas= 14*2 # Molar mass of carrier gas g/mol
# W= 36.6 # W-value of gas in eV 
# SRIM_range = 73e-4 #cm
# SRIM_straggle = 11e-4 #cm



# SF6
GAS='SF6'
r_gas=550e-12/2 # Kinetic radii (NOT DIAMETER!) of carrier gas in (diameter pm/2 e-12) m
M_gas= 146.1 # Molar mass of carrier gas g/mol
W= 34 # W-value of gas in eV 
SRIM_range = 86e-4 #cm
SRIM_straggle = 17e-4 #cm

# # CF4
# GAS='CF4'
# r_gas=470e-12/2 # Kinetic radii (NOT DIAMETER!) of carrier gas in (diameter pm/2 e-12) m
# M_gas= 88 # Molar mass of carrier gas g/mol
# W= 34.4 # W-value of gas in eV 
# SRIM_range = 80e-4 #cm
# SRIM_straggle = 13e-4 #cm


# # He
# GAS='He-II'
# r_gas=260e-12/2 # Kinetic radii (NOT DIAMETER!) of carrier gas in (diameter pm/2 e-12) m
# M_gas= 4 # Molar mass of carrier gas g/mol
# W= 42.7 # W-value of gas in eV 
# SRIM_range = 44.5e-4 #cm
# SRIM_straggle = 33e-4 #cm




# CHEMICAL NUETRALSATION
C_rn = 114000 # radon concentration bq/m3 - enviroment




# =============================================================================
# SIMULATION INPUT - DETECTOR CONFIGURATION
# =============================================================================

# ANSYS ELECTRIC FIELD SOLUTION 
# NLIST and PRNSOL merged
path_to_merge = 'NODE_SOL_MERGE_RAD7.csv'



# =============================================================================
# SIMULATION CONFIGURATION
# Note to find out number of nuclei N converts to RAD7 concentration in (Bq/m3) use NUCLEItoRAD7CONC(N)
# =============================================================================
N = 100000 # Number of radon nuclei
dt =5e-4 # simulation step [s]
totalsimtime = 20  # simulation time [s]


PLOT = 'OFF' # PLOT ON or OFF


#************************ NO USER INPUT REQUIRED FROM HERE************************#


# Load modules here
import rdgmodule
import pandas as pd
from numpy import array
from numpy import append
import csv
import matplotlib.pyplot as plt
from numpy.random import random


# =============================================================================
# DATA ARRAY CONFIGURATION
# =============================================================================
#OUTPUT FILES: PARTICLE NUMBER, N, ENDPOINT, STARTING COORDINATES AND FINAL COORDINATES
outputfile = []


# =============================================================================
# CARRIER GAS PARAMETERS CALCULATIONS
# =============================================================================
# rho_gas - DENSITY OF GAS FOR GIVEN T and P [kg/m3]
# mu - ION MOBILITY IN CARRIER GAS  [cm2/(vs)]
# D = - DIFFUSION COEEFICIENT IN CARRIER GAS [cm2/s]
rho_gas, mu, D = rdgmodule.CarrierGasParameters(M_gas,r_gas, T, Pressure)

# Conversion from nuclei to radon conc in the RAD7 -  1 piC/l = 37 Bq/m3
RAD7conc = rdgmodule.NUCLEItoRAD7CONC(N) #RAD7 concentration in (Bq/m3)

# =============================================================================
# PLOT RAD7 MODEL
# =============================================================================


if PLOT == 'ON':
    plt.axes()
    #Circle
    circle = plt.Circle((0,0),5.1, fc='whitesmoke',ec="None")
    plt.gca().add_patch(circle)
    #Rectangle
    rectangle = plt.Rectangle((0,-5.1),6.1,10.2, fc='whitesmoke',ec="None")
    plt.gca().add_patch(rectangle)
    #Rectangle - detector
    rectangle = plt.Rectangle((3.9,-1.4),1.3,2.8, fc='white',ec="darkgrey")
    plt.gca().add_patch(rectangle)
    #Rectangle - detector holder
    rectangle = plt.Rectangle((5.2,-0.5),0.9,1, fc='white',ec="darkgrey")
    plt.gca().add_patch(rectangle) 
    plt.axis('scaled')
    #Rectangle - dectecto window
    rectangle = plt.Rectangle((3.9,-1.0),0.1,2.0, fc='whitesmoke',ec="darkgrey")
    plt.gca().add_patch(rectangle)



# =============================================================================
# MONTE CARLO SIMULATION 
# =============================================================================
# SINGLE PARTICLE
for i in range(1,N+1):      
    # Random Position Generation (h, r) 
    h_ran,r_ran = rdgmodule.RandomCoordinates()
    P = array([[h_ran],[r_ran]]) # Position (h,r) [cm]
    
    # Save intial position
    h0 = h_ran
    r0 = r_ran
    
    # intial single particle conditions
    V = array([[0], [0]]) # INITIAL VELOCITY (h, r) [cm/s] - assume staionary due to collisions 
    t = 0 # time in simulation
    
    # single particle tracking data with intial parameters 
    particletracking = [[h_ran,r_ran,t]] 
    
    
  # CORRECTION: Neutral Po-218 Species
    ENDPOINT = rdgmodule.NeutralSpeciesCorrection()

  # EVENT: Rn-222 alpha decay (5.49 MeV)
    if ENDPOINT !=  'NEUTRAL SPECIES': # only continue if nuclei is charged
        P = rdgmodule.Rn222Decay(P,SRIM_range, SRIM_straggle)
        particleupdate = append(P,0).tolist() #  # Single particle tracking data t= 0
        particletracking.append(particleupdate) #append to single particle tracking
        Pprev = P # previous position required for EndPoint()
        ENDPOINT = rdgmodule.Boundaries(Pprev, P)

   # SINGLE PARTICLE TRACKING
        if not ENDPOINT in ( 'WALLS', 'COLLECTED', 'DCASE', 'DETECTED'): # only continue if recoil doesnt meet boundary conditions
        
            while t < totalsimtime:
                # Determine Nearest E-field Node
                NEF = rdgmodule.NearestNodeEfield(P[0], P[1])
                # Update new velocity in dt
                Vgained = NEF * mu # Velocity gained from E field
                # V = Vgained + V # Updated total velocity 
                # ******** V IS NOT GAINED it is UPDATE!********* 20221111
                V = Vgained
                
                
                
                # Calculate new position in dt
                Pgained = V * dt # displacement in dt
                Pprev = P # previous position required for EndPoint()
                P = Pgained + P # update absolute position
                
                # Mean square displacement correction in dt
                P = rdgmodule.msdCorrection(P, D, dt) 
                
                # Single particle tracking data
                particleupdate = append(P,t).tolist() # update in dt
                particletracking.append(particleupdate) #append to single particle tracking
                t = t + dt #update simtime by dt
                
                # Nuclei tracking boundaries
                ENDPOINT = rdgmodule.Boundaries(Pprev, P)
                
                
                
# =============================================================================
#                 # Small ion recombination 
# =============================================================================
                # ENDPOINT = rdgmodule.SIR(W, C_rn, mu, dt)
                e=1.602e-19 # electron charge
                E_rn= 5.40e6 # Alpha decay energy of Radon 222 in eV
                eps  = 8.85e-12 #permittivity of free space F/m
                k = 312.5192277538609 # negative ion calibration Raabe (24)
                
                
               
                # half life
                t_hl = 0.693*eps/((E_rn/W*C_rn*k)*e*mu*1e-4)
                
                # tau
                tau = t_hl/0.693
                
                # Probability in dt
                p =  dt/tau
                                
                
            
                
                if random() <=p:  #Probability of SIR for carrier gas and radon concetration  - funtion is true p times of the time
                    ENDPOINT = 'SIR'

                
                
                
                
                if  ENDPOINT in ( 'WALLS', 'COLLECTED', 'DCASE', 'SIR'): #Stop tracking once boundary is reached
                    break
                
    
    # Alpha detection - isotropic alpha emission only half of the collected atoms on the PIPS can be detected
    # if ENDPOINT == 'COLLECTED':
        # ENDPOINT = rdgmodule.PIPSdetection()
    
    # Nuclei endpoint
    print(i, ENDPOINT)
    
    

    # df for single particle tracking
    singledata = pd.DataFrame(particletracking)
    singledata.columns = ['h (cm)', 'r (cm)', 'simtime (s)']
    
    
    # Plotting tracks 
    if PLOT == 'ON':
        # ENDPOINT ='COLLECTED'
        rdgmodule.TrackPlot(ENDPOINT, singledata)
        plt.show()
    
    # Simulation data - save simdatabackup.csv
    OutputUpdate = [i,ENDPOINT,t,h0,r0,h_ran,h_ran]
    outputfile.append(OutputUpdate) #append to outputfile
    with open('Output/'+GAS+'_'+str(N)+'_simdatabackup.csv', "w") as f: #csv backup
        writer = csv.writer(f)
        writer.writerows(outputfile)
        

 

  
    del singledata

# Simulation data output
simdata = pd.DataFrame(outputfile)
simdata.columns = ['Nuclei','Endpoint','flighttime [s]','h0 [cm]','r0 [cm]','hf [cm]','hf [cm]']


# plot legends
plt.plot([],[], linestyle='-.', marker='o',markersize=1, markerfacecolor='white', color='#D01E2F',alpha=0.9, linewidth=1, label = 'Lost')
plt.plot([],[], linestyle='-.', marker='o',markersize=1, markerfacecolor='white', color='black',alpha=0.9, linewidth=1, label = 'Collected')
plt.plot([],[], linestyle='none',markersize=2, marker="o", markerfacecolor='orange', color='orange',alpha=0.9, label = 'Neutral $^{218}$Po')
plt.plot([],[], linestyle='-.', marker='o',markersize=1, markerfacecolor='white', color='yellow',alpha=0.9, linewidth=1, label = 'Small Ion Recombination')
plt.ylabel('r dimension (cm)')
plt.xlabel('h dimension (cm)')

plt.legend()


print(simdata)  # print simdata

# Post simlation calculation 
N_DETECTED = (simdata['Endpoint'] == ('DETECTED') ).sum()
N_COLLECTED = (simdata['Endpoint'] == 'COLLECTED').sum() + N_DETECTED

N_SIR = (simdata['Endpoint'] == ('SIR') ).sum()
# CPM = rdgmodule.NUCLEItoRAD7CONC(N_DETECTED) * 1e-3 * 60 #Bq/m3 *m3 * 60s
# RAD7Sensitivity = CPM/RAD7conc
CollectEFF = N_COLLECTED/N
N_lost = N- N_COLLECTED


flight = simdata['flighttime [s]'].mean()

# Simulation config data output 
with open('Output/'+GAS+'_'+str(N)+'_simoutput.txt', 'w') as f: # save simconfig.txt
    print(
          
          '*****SIMULATION OUTPUT****',
          '\n \nCARRIER GAS CONFIGURATION',
          '\nGas:', GAS,
          '\nSRIM range (cm):',SRIM_range, '±', SRIM_straggle,
          '\nKinetic radii (cm):',r_gas,
          '\nMolar mass (g/mol):', M_gas,
          '\nIonisation energy (eV):',IE,
          '\nElectron affinity (eV):', EA,
          '\nCalculated mobility(cm2/(vs)):',mu,
          '\nCalculated density (kg/m3):', rho_gas,
          '\n \nENVIROMENTAL PARAMETERS',
          '\nTemp (K):',T,
          '\nPressure (pa):', Pressure,
          '\nRelative humidity (%):', RH,
          '\n \nSIMULATION CONFIGURATION',
          '\nNumber of Nuclei:',N,
          # '\nRAD7 Concentration(Bq/m3):',RAD7conc,
          '\nSimulation time step (s):', dt,
          '\nANSYS solution file:', path_to_merge,
          '\n \nPOSTSIMULATION CALCULATIONS',
          '\nNuclei collected:', N_COLLECTED,    
                    '\nNuclei lost:', N_lost,     
            '\nNuclei SIR:', N_SIR,    
           # '\nNuclei detected:', N_DETECTED,
          '\nCollection Efficiency):', CollectEFF,
          '\nAverage Time of Flight (s):', flight,
          # '\nRAD7 sensitivity CPM/(Bq/m3):', RAD7Sensitivity,
          
          file=f
          )

print(  # print simconfig.txt
          
          '*****SIMULATION OUTPUT****',
          '\n \nCARRIER GAS CONFIGURATION',
          '\nGas:', GAS,
          '\nSRIM range (cm):',SRIM_range, '±', SRIM_straggle,
          '\nKinetic radii (cm):',r_gas,
          '\nMolar mass (g/mol):', M_gas,
          '\nIonisation energy (eV):',IE,
          '\nElectron affinity (eV):', EA,
          '\nCalculated mobility(cm2/(vs)):',mu,
          '\nCalculated density (kg/m3):', rho_gas,
          '\n \nENVIROMENTAL PARAMETERS',
          '\nTemp (K):',T,
          '\nPressure (pa):', Pressure,
          '\nRelative humidity (%):', RH,
          '\n \nSIMULATION CONFIGURATION',
          '\nNumber of Nuclei:',N,
          # '\nRAD7 Concentration(Bq/m3):',RAD7conc,
          '\nSimulation time step (s):', dt,
          '\nANSYS solution file:', path_to_merge,
          '\n \nPOSTSIMULATION CALCULATIONS',
          '\nNuclei collected:', N_COLLECTED,   
          '\nNuclei lost:', N_lost,     
                 '\nNuclei SIR:', N_SIR,    
          # '\nNuclei detected:', N_DETECTED,
        '\nCollection Efficiency):', CollectEFF,
        '\nAverage Time of Flight (s):', flight,
          # '\nRAD7 sensitivity CPM/(Bq/m3):', RAD7Sensitivity
          
          )



    
