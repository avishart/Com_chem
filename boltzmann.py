#Made by Andreas Vishart
#Uses a list of a output files (*out)
import sys
import numpy as np


#-------Parameters-------
#The type of Energy 
Energy="thermal Free Ener"
#Temperature
Temperature=298.15
#The name of the txt file with all info
filename="boltzmann_files.txt"

#Constants
Har_to_J=4.3597482*10**(-18)
kB=1.380658*10**(-23)
kBT=kB/Har_to_J*Temperature


#------Program------

files_list=[]
energy_list=[]
#Get energy from all files and their names
for numsys in range(1,len(sys.argv)):
    energy=[]
    with open(sys.argv[numsys]) as thefile:
        content=thefile.readlines()
    for line in content:
        if Energy in line:
            newline=list(filter(None,line.split(" ")))
            energy.append(float(newline[7]))
    if len(energy)>0:
        energy_list.append(energy[-1])
        files_list.append(sys.argv[numsys])

#Get the order
order_list=sorted(zip(energy_list,files_list))
#Sort the energy from smallest to highest and the filenames
energy_order,file_order=[list(k) for k in zip(*order_list)]

#Calculate the boltzmann distribution
Q=0
for ener in energy_order:
    Q+=np.exp(-(ener-energy_order[0])/kBT)
distribution=[np.exp(-(ener-energy_order[0])/kBT)/Q for ener in energy_order]

#Generate a txt file with all details
newfile=open(filename,"w")
for file in range(len(file_order)):
    newfile.write(file_order[file]+"\n")
    newfile.write(str(energy_order[file])+"  ;  "+str(distribution[file]*100)+" %\n\n")
newfile.close()

#print the three files with lowest energies and highest populations
print(str(file_order[0])+" ; "+str(energy_order[0])+" A.U. ; "+str(distribution[0]*100)+" %")
try:
    print(str(file_order[1])+" ; "+str(energy_order[1])+" A.U. ; "+str(distribution[1]*100)+" %")
    print(str(file_order[2])+" ; "+str(energy_order[2])+" A.U. ; "+str(distribution[2]*100)+" %")
except:
    print("Only one or two conformers!")