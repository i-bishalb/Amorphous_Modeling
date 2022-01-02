'''
This code is written by:   Bishal Bhattarai
                    Email: bishal_bhattarai@outlook.com
-------------------------------------------------------------------------------------------------------------------
The main purpose of this code is to create a Random configuration to begin molecular dynamics simulation with VASP.
A POSCAR file is generated.


The code takes basic input from the user i.e. Atom symbol, number of atoms, number of each atom in the system and density
and creates random configuration.

Each atom in the final random configuration are separated by 1.80 Angstrom (can be changed by modifying MINIMUM)

------------------------------------------------------------------------------------------------------------------
Pymatgen installation is required.

'''




import random
import time
from decimal import Decimal, ROUND_HALF_UP
import math
import os
import shutil
from pymatgen.core.periodic_table import Element


MINIMUM = 1.80


ONES  = 1.00000000
ZEROS = 0.00000000


#seed the random number
random.seed(time.time())


TYPE = input ("Please input total NUMBER of ELEMENTS in the system (integer):   ")
print("\n")
TYPE = int(TYPE)

series_count = input ("Please input FOLDER COUNT (integer):   ")
print("\n")
series_count = int(series_count)


SYM = []
COUNT = []
density = 0.0

for i in range(0, TYPE):
    symbol = input("Please input SYMBOL of element %d:   "%(i+1))
    elementcount = input ("Please input COUNT of the element %s:   "%symbol)
    print("\n")
    SYM.append(symbol)
    COUNT.append(int(elementcount))

density = input("Please input Density of the material in g/cc ")
density = float(density)

MASS = []

for i in range(0, len(SYM)):
    mass = Element(SYM[i]).atomic_mass
    MASS.append(float(mass))



total_mass = 0.0
for i in range(0, len(COUNT)):
    total_mass += COUNT[i]*MASS[i]

Volume = total_mass *1.661/(density)

BOXSIZE = Volume**(1.0/3.0)


BOXSIZE = Decimal('%.3f' % (BOXSIZE * 1000 / 1000))


CURRENT_PATH = os.getcwd()

folder_name_string = '/' + 'S'+ str(series_count).zfill(2) + '_lat' + str(BOXSIZE) + '_' + 'd' + str(density)

NEW_FOLDER = CURRENT_PATH + folder_name_string

os.makedirs(NEW_FOLDER)
os.chdir(NEW_FOLDER)

NATOMS = sum(COUNT)

FINAL_X = []
FINAL_Y = []
FINAL_Z = []


TEMP_X = []
TEMP_Y = []
TEMP_Z = []


angle1 = 2.0* math.pi
angle2 = 4.0 * math.pi


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def CHECK_PBC(x,y,z):
    value = 0
    if( x <= float(BOXSIZE) - MINIMUM/2.0 and y <= float(BOXSIZE) - MINIMUM/2.0 and z <= float(BOXSIZE) - MINIMUM/2.0 and x >= MINIMUM/2.0 and y >= MINIMUM/2.0 and z >= MINIMUM/2.0):
        value = 1
    return value


def CHECK_DISTANCE(x,y,z, Xarray, Yarray, Zarray):

    distance = 0
    value_raised = 1

    DIST = []
    for j in range(0, len(Xarray)):
        distance = (Xarray[j] - x)** 2.0 + (Yarray[j] - y)** 2.0 + (Zarray[j] - z)**2.0
        distance = math.sqrt(distance)
        DIST.append(distance)

    for elmnts in DIST:
        if(elmnts < MINIMUM):
            value_raised = 0

    return value_raised


def SHUFFLE_XYZ(Xarray, Yarray, Zarray):


    combined = list(zip(Xarray, Yarray, Zarray))
    random.shuffle(combined)
    Xarray[:], Yarray[:], Zarray[:] = (zip(*combined))

    return Xarray, Yarray, Zarray



def MAKE_POSCAR(Xarray, Yarray, Zarray, filename):

    print("***** WRITING POSCAR FILE****** \n")

    start = 0
    end = 0

    f = open(filename, 'w')
    f.write(" Amorphous ")
    for t in range(0,len(SYM)):
        f.write("%2s%d"%(SYM[t],COUNT[t]))
    f.write(" Density= %f"%(density))
    f.write("\n")


    f.write("%8.5f\n"%(ONES))
    f.write("%8.3f %8.3f %8.3f\n"%(float(BOXSIZE), ZEROS, ZEROS))
    f.write("%8.3f %8.3f %8.3f\n"%(ZEROS, float(BOXSIZE), ZEROS))
    f.write("%8.3f %8.3f %8.3f\n"%(ZEROS, ZEROS, float(BOXSIZE)))

    for t in range(0,len(SYM)):
        f.write("%5s"%(SYM[t]))
    f.write("\n")

    for t in range(0,len(SYM)):
        f.write("%5d"%(COUNT[t]))
    f.write("\n")
    f.write("Direct\n")

    NATOMS = sum(COUNT)

    for q in range(0, NATOMS):
        f.write("%14.8f %14.8f %14.8f \n"%(FINAL_X[q]/float(BOXSIZE), FINAL_Y[q]/float(BOXSIZE), FINAL_Z[q]/float(BOXSIZE)))

    f.close()

    return None



#--------------------------------------------------------------------------------------------------------------------------------------



i=0
x = y = z = 0.0

while(i<NATOMS):
    x = y = z = 0.0
    r0=0.0

    test = 0
    double_test = 0

    test1 = 0
    test2 = 0
    distance = 0.0

    if(i==0):
        while(test==0):
            phi =  angle1* random.random()
            theta = angle2 * random.random()
            r0 = 25* float(BOXSIZE) * random.random()

            x = r0 * math.cos(phi) * math.sin(theta)
            y = r0 * math.sin(phi) * math.sin(theta)
            z = r0 * math.cos(theta)

            test = CHECK_PBC(x,y,z)

        FINAL_X.append(x)
        FINAL_Y.append(y)
        FINAL_Z.append(z)

    else:
        while(double_test==0):

            if(i<int(NATOMS* 0.60)):
                r0 = 5.00 * float(BOXSIZE) * random.random()
                phi =  angle1* random.random()
                theta = angle2 * random.random()

            else:
                r0 = 1.5 * float(BOXSIZE) * random.random()
                phi = 1.0* angle1* random.random()
                theta = 1.0* angle2 * random.random()


            x = r0 * math.cos(phi) * math.sin(theta)
            y = r0 * math.sin(phi) * math.sin(theta)
            z = r0 * math.cos(theta)

            test1 = CHECK_PBC(x,y,z)
            test2 = CHECK_DISTANCE(x,y,z, FINAL_X, FINAL_Y, FINAL_Z)

            double_test = test1 * test2

        FINAL_X.append(x)
        FINAL_Y.append(y)
        FINAL_Z.append(z)

    i+=1

    print(i)



filename='POSCAR'
MAKE_POSCAR(FINAL_X, FINAL_Y, FINAL_Z, filename)

print('The Random poscar has been created in the folder %s \n'%folder_name_string)
