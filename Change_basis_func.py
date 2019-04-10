#Imports 
import numpy as np
import copy as copy
import sys

#Parameters

title = (sys.argv[1])[:-4]


functionals = ["bp86", "rbvp86", "cam-b3lyp", "pbe1pbe"]
basis_set = ["6-31G", "6-31++G", "6-31G(d,p)", "6-31++G(d,p)", "6-311G", "6-311++G", "6-311G(d,p)", "6-311++G(d,p)"]
#functionals = ["BP86"]
#basis_set = ["6-31G", "6-31++G"]

linje = 2


with open(sys.argv[1]) as thefile:
    content = thefile.readlines()
length = len(content[linje])-1
for i in functionals:
    for j in basis_set:
        cpcontent = copy.copy(content)
        cpcontent[linje] = cpcontent[linje][:length]
        cpcontent[linje] = cpcontent[linje]+" "+i+"/"+j+'\n'
        if "(d,p)" in j:
            k = j.replace("(d,p)","dp")
            file = open(title+"_"+i+"_"+k+".com", "w")
        else:
            file = open(title+"_"+i+"_"+j+".com", "w")
        for k in range(0,len(cpcontent)):
            file.write(cpcontent[k])
        file.close()
            

