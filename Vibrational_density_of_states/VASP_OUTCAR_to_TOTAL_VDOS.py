'''
This code is written by:   Bishal Bhattarai
                    Email: bishal_bhattarai@outlook.com
-------------------------------------------------------------------------------------------------------------------
Read phonon frequencies from OUTCAR and compute total density of states for amorphous materials
------------------------------------------------------------------------------------------------------------------
Pymatgen installation is required.


Gaussian broadening is used to calulate VDOS. Please read more about it on:
    B. Bhattarai and D. A. Drabold, “Vibrations in amorphous silica”, Journal of Non-Crystalline Solids,
    439, 6-14 (2016).
    B. Bhattarai, R. Thapa and D. A. Drabold, “Ab initio inversion of structure and the lattice dynamics of
    a metallic glass: The case of Pd40Ni40P20”, Modelling Simul. Mater. Sci. Eng, 27, 075002 (2019).

'''


from pymatgen.core.lattice import Lattice
from pymatgen.io.vasp.inputs import Poscar         #importing Pymatgen Poscar object
from pymatgen.io.vasp.outputs import Outcar
import numpy as np
import matplotlib.pyplot as plt
import math
import os,shutil
import glob
import io


#-----------------------------------------Reading Contcar file as pymatgen object -------------------

contcar = Poscar.from_file('CONTCAR')

atoms_count = contcar.natoms
natoms = sum(atoms_count)


#-----------------------------------------Reading outcar file for frequencies ---------------------
outcar_lines = []
with open('OUTCAR', 'r') as f:
        for line in f:
                outcar_lines.append(line)
f.close()

lookup= ' Eigenvectors and eigenvalues of the dynamical matrix\n' #Start point phrase
startpoint = 0
for i in range(0 , len(outcar_lines)):
        if(outcar_lines[i]== lookup):
                startpoint = int(i + 4)

endpoint = int(startpoint+ 3 * natoms*(natoms+ 3) - 1)

freq_cm_inv_values = []
j=0
for i in range(startpoint, endpoint,natoms+3):
     b= outcar_lines[i].split()
     if(j < 3 * natoms -3):
         try:
             freq_cm_inv_values.append(float(b[7]))   # READING FREQUENCY in cm-1 from OUTCAR
         except ValueError:
             print('\n **** ERROR---> Imaginary frequencies encountered -- system not relaxed well \n')
             print('\n ****           Re-Relax your system and run phonon calculations again')
             freq_cm_inv_values.append(float(b[6]))
     else:
          freq_cm_inv_values.append(float(b[6]))   # READING FREQUENCY in cm-1 from OUTCAR

     j+=1
#-----------------------------------------Calculation phonons -----------------------------------------------------------

filename = 'Total_phonon_vdos.dat'
if os.path.exists(filename):
    os.remove(filename)

outf = open('Total_phonon_vdos.dat', 'w')
outf.write('# Frequency (cm-1)         VDOS (arb. units) \n')


'''
     Gaussian function = 1/(sigma * sqrt(2* pi)) exp (- (x-mu)**2 / 2 *(sigma**2))
                       = coef * exp (- (x - mu)**2 /deno )
                       = coef * exp (- (numo * numo) / deno)

'''

sigma = 0.02 * max(freq_cm_inv_values)                 # Taking sigma as 2 percent of maximum frequency
coef = 1.0 / (sigma * math.sqrt(2.0 * math.pi))
deno = 2.0 * sigma * sigma

max_frequency = max(freq_cm_inv_values) + 0.20 * max(freq_cm_inv_values) # going upto 20 % above maximum frequency

k = 0
total_sum = 0
gaussian = 0

total_freq = []
total_vdos = []
while (k < max_frequency):
    for i in range(0, len(freq_cm_inv_values) - 3): #ignoring first 3 frequencies
        numo = k - freq_cm_inv_values [i]
        exponent = math.exp(-1.0 * (numo * numo) / deno)
        gaussian = coef * exponent
        total_sum += gaussian

    total_freq.append(k)
    total_vdos.append(total_sum)
    outf.write(' %15.10f   %15.10f \n'%(k,   total_sum/(3.0 * natoms)))
    total_sum = 0
    gaussian = 0
    k += 0.20

outf.close()


#----------------------------------------------------Plotting phonons ---------------------------------

# Some plot prperties ------------------------------------
scaling = 1.5

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.linewidth'] = 10*scaling
plt.rcParams['font.size'] =40*scaling
plt.rcParams['axes.labelsize'] = 52 *scaling
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 55 * scaling
plt.rcParams['axes.titleweight']='bold'
plt.rcParams['xtick.labelsize'] =35*scaling
plt.rcParams['ytick.labelsize'] =35* scaling
plt.rcParams['lines.linewidth']= 9.0* scaling
plt.rcParams['figure.figsize'] = 21*scaling,18*scaling
plt.rcParams['lines.markersize']=12 * scaling
plt.rcParams['axes.linewidth'] = 10.0* scaling

plt.rcParams['xtick.major.size']=20 * scaling
plt.rcParams['xtick.direction']= 'in'
plt.rcParams['xtick.major.width']=4.0 * scaling



plt.rcParams['ytick.major.size']=20 * scaling
plt.rcParams['ytick.direction']= 'in'
plt.rcParams['ytick.major.width']=4.0 * scaling

plt.rcParams['xtick.major.pad']='15'
plt.rcParams['ytick.major.pad']='15'


fig, ax = plt.subplots(1, 1)
ax.plot(total_freq, total_vdos,color='#f4a582',linestyle='-',label= "Total VDOS", linewidth = 7.0, zorder = 0)
ax.fill_between(total_freq, 0, total_vdos, color='#f4a582', alpha = 0.6)
ax.legend(loc='upper right', ncol=1, bbox_to_anchor=(0.60, 0.625, 0.375, 0.375),
         labelspacing=0.25, columnspacing=0.20,prop={'weight':'bold', 'size': 65})



x_lim = np.arange(0.0, max(total_freq) + 0.05 * max(total_freq), 200)
x_lim_list = x_lim.tolist()
ax.set_xticks(x_lim_list)
ax.set_ylim(0,np.max(total_vdos) + 0.05 * max(total_vdos))


ax.set_xlabel("Frequency $\mathbf{(cm^{-1})}$",rotation=0, fontweight='bold', fontsize = 72)
ax.set_ylabel("VDOS (arb. units)",rotation=90, fontweight='bold', fontsize = 72)
ax.grid(linestyle=':',linewidth='3.0')


plt.savefig("Total_VDOS.pdf",bbox_inches='tight',dpi=2000,pad_inches=0.25)
plt.close()
