# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 12:28:29 2021
Please use the following code to generate the atomic velocity trajectory file:
#############################################################################
dump dump_id dump_group custom 1 dump.lammpstrj id type vx vy vz
dump_modify dump_id sort id
#############################################################################
@author: benward (hityingph@163.com)
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os
def find_pdos(v_all, Nc, dt, omega):    #Calculate the vacf from velocity data
    "The oringal code was written in Matlab by Fan Zheyong"
# v_all: all the velocity data in some format
# Nc: number of correlation steps
# dt: time interval between two frames, in units of ps
# omega: phonon angular frequency points you want to consider
    Nf = v_all.shape[0]                 # number of frames
    M  = Nf - Nc                        # number of time origins for time average
    vacf = np.zeros(Nc)                 # the velocity autocorrelation function (VACF)
    for nc in range(Nc):                # loop over the correlation steps
        ratio = (nc+1)/Nc * 100    
        print("Calculate PDOS Progress %s%%" %ratio)
        for m in range(M+1):            # loop over the time origins
            delta = np.sum(v_all[m + 0]*v_all[m + nc])
            # print(delta)
            vacf[nc] = vacf[nc] + delta

    
    vacf = vacf / vacf[0]                                       # normalize the VACF
    vacf_output = vacf                                          # copy the VACF before modifying it
    vacf = vacf*(np.cos(np.pi*np.arange(Nc)/Nc)+1)*0.5          # window function
    vacf = vacf*np.append(np.ones(1), 2*np.ones(Nc-1))/np.pi    # C(t) = C(-t)
    pdos = np.zeros(len(omega))                                 # the phonon density of states (PDOS)
    for n in range(len(omega)):                                 # Discrete cosine transform
        pdos[n] = dt * sum(vacf * np.cos(omega[n] * np.arange(Nc) * dt))
    return(vacf_output, pdos)
### set up some parameters related to the MD simulation
file_list = [
    r"E:\MOF\ALL USE UFF4MOF\VDOF_UFF4MOF\PDOS by python\UIO66\O_2.lammpstrj",
    r"E:\MOF\ALL USE UFF4MOF\VDOF_UFF4MOF\PDOS by python\UIO66\O_3_f.lammpstrj",
    r"E:\MOF\ALL USE UFF4MOF\VDOF_UFF4MOF\PDOS by python\UIO66\C_R.lammpstrj",
    r"E:\MOF\ALL USE UFF4MOF\VDOF_UFF4MOF\PDOS by python\UIO66\Zr8f4.lammpstrj",
]
save_list = [
    'C_2.lammpstrj',
    'O_2.lammpstrj',
    'C_R.lammpstrj',
    'Zn2f2.lammpstrj'
]
filepath  = r'E:\MOF\ALL USE UFF4MOF\VDOF_UFF4MOF\PDOS by python\Zn-MOF-74'
for save_file in save_list:
    fileName = filepath + '/' + save_file
    with open (fileName) as f:
        for j in range(3):                          # The first 9 lines should be excluded
            f.readline()
        N = f.readline().strip()
        N = int(N)

    #N = 1296                          # number of atoms
    Nc = 2000                           # number of correlation steps (a larger number gives a finer resolution)
    dt = 0.001                         # time interval in units of ps (its inverse is roughly the maximum frequency attainable)
    omega = np.arange(1, 380.5, 0.5)    # angular frequency in units of THz
    nu = omega / 2 / np.pi;             # omega = 2 * pi * nu, while nu is the frequency range
    Nf = Nc*10                          # The calculation numbers
    num_frame = Nf                      # Frame numbers equals to calculation numbers
    #fileName = r"E:\MOF\ALL USE UFF4MOF\VDOF_UFF4MOF\PDOS by python\UIO66\O_2.lammpstrj" #Read the dump file
    ### Read the atomic velocity data from the lammps dump file
    # Noted that the atomic ID should be sorted, set "dump_modify dump_id sort id" in lammps input file
    v_all = np.zeros((num_frame, N, 3))             # Initialize the velocity data
    fin = open(fileName, "r")
    for i in range(num_frame):
        print(i)
        ratio = (i+1)/num_frame * 100 
        print("Read Data Progress %s%%" %ratio)
        initial = i * (9 + N)
        for j in range(9):                          # The first 9 lines should be excluded
            fin.readline()
            t=1
        for k in range(N):
            line = fin.readline().split()[2:]
            line = [float(l) for l in line]
            v_all[i, k] = line

    vacf, pdos = find_pdos(v_all, Nc, dt, omega)    # Call the function and calculate the vacf and pdos
    t = np.arange(Nc)*dt                            # The discrete correlation time
    os.chdir('{}'.format(filepath))
    os.system('mkdir {}result'.format(save_file[:-10]))
    os.chdir('{}result/'.format(save_file[:-10]))
    vacf_put = pd.DataFrame({ 'Correlation Time (ps)' : t, 'Normalized VACF': vacf})
    vacf_put.to_csv(r"VACF.csv")
    pdos_put = pd.DataFrame({ '$\omega$ (THz)' : nu, 'PDOS': pdos})
    pdos_put.to_csv(r"PDOS.csv")
    os.chdir('{}'.format(filepath))



# ### Plot the figure
# plt.figure(figsize=(12,5))
# plt.subplot(1, 2, 1)
# plt.plot(t, vacf, linewidth = 2, color="C1")
# plt.xticks(fontsize = 14)
# plt.xlabel('Correlation Time (ps)', fontsize = 14)
# plt.yticks(fontsize = 14)
# plt.ylabel('Normalized VACF', fontsize = 14)

# plt.subplot(1, 2, 2)
# plt.plot(nu, pdos, linewidth = 2, color="C4")
# plt.xticks(fontsize = 14)
# plt.xlabel('$\omega$ (THz)', fontsize = 14)
# plt.yticks(fontsize = 14)
# plt.ylabel('PDOS', fontsize = 14)

# plt.subplots_adjust(wspace=0.3)
# plt.savefig("PDOS.pdf", bbox_inches='tight')
# plt.show()