#Made by Andreas Vishart
#Use one folder with output files with single point calculation, 
#so use python script in folder
#The two atoms in the bond distance is chosen for atom1 and atom2

#------------Package---------------------------
import numpy as np
import sys
import copy as cp
import commands as com
import matplotlib.pyplot as plt
#import subprocess

#------------Parameters---------------------------
#The first atom in the bond
atom1=2
#The second atom in the bond
atom2=4
#Data file name
data="Extract_dis_energy"


file_names=list(com.getstatusoutput("grep \"Normal termination\" *.out"))
#print(file_names)
file_names=file_names[1].split("\n")
#file_names=subprocess.call("grep \"Normal termination\" *.out", shell=True)
data_file_name=data+"_"+file_names[0].split(" ")[0][:-5]+".txt"
data_file=open(data_file_name,"w")
dis_list=[]
energy_list=[]
for line in file_names:
    ori=0
    xyz=[]
    name=line.split(" ")[0][:-1]
    with open(name) as thefile:
        content=thefile.readlines()
    for text in content:
        if "SCF Done:" in text:
             energy=float(text.split(" ")[7])
             energy_list.append(energy)
        elif "Input orientation:" in text:
            ori=1
        if ori==1:
            xyz.append(text)
        if "Distance matrix" in text:
            ori=0
    xyz=xyz[5:-2]
    xyz1=map(float,filter(None,xyz[atom1-1].split(" ")))
    xyz2=map(float,filter(None,xyz[atom2-1].split(" ")))
    dis2=(xyz2[3]-xyz1[3])**2+(xyz2[4]-xyz1[4])**2+(xyz2[5]-xyz1[5])**2
    dis=dis2**0.5
    dis_list.append(dis)
    data_file.write("{:.8f}".format(round(dis,8))+"; "+"{:.8f}".format(round(energy,8))+"\n")
data_file.close()

print("File "+data_file_name+" saved!")

plt.plot(dis_list,energy_list,"bo")
plt.xlabel("Displacement / [$\AA$]")
plt.ylabel("Energy / [Hartree]")
plt.show()
