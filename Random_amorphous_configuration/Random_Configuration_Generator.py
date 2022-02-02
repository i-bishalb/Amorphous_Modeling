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
import numpy as np


minimum = 1.80


ones  = 1.00000000
zeros = 0.00000000


#seed the random number
random.seed(time.time())


#------------------------------------------------Taking User inputs-----------------------------
type_input = input ("Please input total number of elements in the material (integer):   ")
print("\n")
type_input = int(type_input)



element_symbol = []
element_count = []
density = 0.0

for i in range(0, type_input):
    sym = input("Please input SYMBOL (eg: Si, O, Ga) of element %d:   "%(i+1))
    count = input ("Please input COUNT of the element %s:   "%sym)
    print("\n")
    element_symbol.append(sym)
    element_count.append(int(count))

density = input("Please input Density of the material in g/cc ")
density = float(density)
#--------------------------------------------------------------------------------------------------


mass = [float(Element(element_symbol[i]).atomic_mass) for i in range(0, len(element_symbol))] #individual atomic mass
total_mass = sum([element_count[i]*mass[i] for i in range(0, len(element_count))])            # total atomic mass
volume = total_mass *1.661/(density)   #volume of the system

boxsize = volume**(1.0/3.0)   #dimension of cubic boxsize
boxsize = round(boxsize, 2)   #rounding off boxsize to 2 floating point values
total_natoms = sum(element_count)


#final coordinates list
final_x = []
final_y = []
final_z = []



#angle theta and phi for spherical coordinates
angle1 = 1.0 * math.pi   #theta
angle2 = 2.0 * math.pi   # phi


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


def check_pbc(x,y,z):
    '''
    Takes coordinates and checks periodic boundary conditions
    eg:
    if minimum = 2.0 and boxsize = 10
    we want our coordinates between 0 + 1.0 and 10 - 1.0 i.e. 1 and 9
    such that the minimum distance of 2 is maintained between end elements
    i.e. 9 and 1 --> 9 to 10: 1 and 10 (0)--> 1
    '''

    condition1 = minimum/2.0
    condition2 = boxsize - minimum/2.0

    if (x <= condition1 or y <= condition1 or z <= condition1):
        return 0

    if (x >= condition2 or y >= condition2 or z >= condition2):
        return 0

    return 1



def check_distance(x,y,z, x_list, y_list, z_list):
    ''' Funtion to check cartesian distance between
       x ,y , z and previously valid coordinates in the list'''

    for j in range(0, len(x_list)):
        distance = 0
        distance = (x_list[j] - x)** 2.0 + (y_list[j] - y)** 2.0 + (z_list[j] - z)**2.0
        distance = math.sqrt(distance)
        if distance <= minimum:
            return 0

    return 1


def shuflle_xyz(x_list, y_list, z_list):
    '''
    This functions suffles the ordering of (x_list, y_list, z_list)
    using the zip function.
    '''
    combined = list(zip(x_list, y_list, z_list))
    random.shuffle(combined)
    x_list[:], y_list[:], z_list[:] = (zip(*combined))

    return x_list, y_list, z_list



def make_poscar(x_list, y_list, z_list, filename):

    '''
    This function writes the poscar file to the directory
    after all the coordinates for the atomic configuration has passed
    the conditions.
    '''

    print("***** WRITING POSCAR FILE****** \n")

    start = 0
    end = 0

    f = open(filename, 'w')
    f.write(" Amorphous ")
    for t in range(0,len(element_symbol)):
        f.write("%2s%d"%(element_symbol[t],element_count[t]))
    f.write(" Density= %f"%(density))
    f.write("\n")


    f.write("%8.5f\n"%(ones))
    f.write("%8.3f %8.3f %8.3f\n"%(float(boxsize), zeros, zeros))
    f.write("%8.3f %8.3f %8.3f\n"%(zeros, float(boxsize), zeros))
    f.write("%8.3f %8.3f %8.3f\n"%(zeros, zeros, float(boxsize)))

    for t in range(0,len(element_symbol)):
        f.write("%5s"%(element_symbol[t]))
    f.write("\n")

    for t in range(0,len(element_count)):
        f.write("%5d"%(element_count[t]))
    f.write("\n")
    f.write("Direct\n")

    natoms = sum(element_count)

    for q in range(0, total_natoms):
        f.write("%14.8f %14.8f %14.8f \n"%(x_list[q]/float(boxsize), y_list[q]/float(boxsize), z_list[q]/float(boxsize)))

    f.close()

    return None



#--------------------------------------------------------------------------------------------------------------------------------------



i=0
x = y = z = 0.0

while(i < total_natoms):
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
            r0 = 25* float(boxsize) * random.random()

            x = r0 * math.cos(phi) * math.sin(theta)
            y = r0 * math.sin(phi) * math.sin(theta)
            z = r0 * math.cos(theta)

            test = check_pbc(x,y,z)

        final_x.append(x)
        final_y.append(y)
        final_z.append(z)

    else:
        while(double_test==0):

            if(i<int(total_natoms* 0.60)):
                r0 = 5.00 * float(boxsize) * random.random()
                phi =  angle1* random.random()
                theta = angle2 * random.random()

            else:
                r0 = 1.5 * float(boxsize) * random.random()
                phi = 1.0* angle1* random.random()
                theta = 1.0* angle2 * random.random()


            x = r0 * math.cos(phi) * math.sin(theta)
            y = r0 * math.sin(phi) * math.sin(theta)
            z = r0 * math.cos(theta)

            test1 = check_pbc(x,y,z)
            test2 = check_distance(x,y,z, final_x, final_y, final_z)

            double_test = test1 * test2

        final_x.append(x)
        final_y.append(y)
        final_z.append(z)

    i+=1

    print(i)



filename='POSCAR'

for ii in range(0, int(total_natoms/5)):
    final_x, final_y, final_z = shuflle_xyz(final_x, final_y, final_z)


make_poscar(final_x, final_y, final_z, filename)
print('The Random poscar has been created !! \n')


#----------------------------------------For Testing--------------------------


# atoms_symbol = []
# atoms_symbol += [[str(element_symbol[i])]*element_count[i] for i in range(0, len(element_count))]
#
# total_atoms_symbol = []
#
# for x in atoms_symbol:
#     total_atoms_symbol.extend(x)
#
# f = open("XYZ.xyz", 'w')
# f.write("%3d"%total_natoms)
# f.write("\n\n")
# for ii in range(0, len(final_x)):
#     f.write("%3s %15.5f %15.5f %15.5f"%(total_atoms_symbol[ii], final_x[ii], final_y[ii], final_z[ii]))
#     if ii != len(final_x)-1:
#         f.write("\n")
#
# f.close()
