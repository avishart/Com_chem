#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import sys


with open(sys.argv[1]) as file:
    content=file.readlines()

Iteration=[]
iter=0
Energy=[]
for i in content:
    if "SCF Done:" in i:
        con=list(filter(None,i.split(" ")))[4]
        Energy.append(float(con))
        Iteration.append(iter)
        iter+=1

plt.plot(Iteration,Energy,"bo-")
plt.xlabel("Iteration steps")
plt.ylabel("Energy / [A.U.]")
plt.show()
plt.clf()

#print(Energy,Iteration)
print("Number of iteration done: "+str(iter))
try:
  print("Energy change: "+str(100*abs(float(con)-Energy[-2])/Energy[-2])+" %")
except:
  print("Only 1 iteration is done!")

