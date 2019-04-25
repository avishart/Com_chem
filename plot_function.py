#Made by Andreas L. Vishart
#Correct the list of functions by writting the function formula
#------------------------Packages------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import special


#------------------------Parameters----------------------
#X start value
x_start=-75
#X end value
x_end=75
#List of function formulas as string with x as "x"
function_formula=["(np.sin(x/2))**2/x**2"]
#Labels
label=["","$D^++A^-$","$E(D^+|D)+E(A^-|A)$"]
#The name of saved plot
title="delta_function_sin2.png"
#X-axis title 
x_label="$E_{if}t/\hbar$"
#Y-axis title
y_label="$D_tt^2/\hbar^2$"
#List of colors for plotting
color=["blue","green","red","orange","pink","purple","black","yellow"]
#plot y values in limit
y_limit=[0,0.26]
#Font
font = {'fontname':'Serif'}

#------------------------Functions-----------------------
def y_func(x_list,function):
    y_list=[eval(function) for x in x_list]
    return y_list



#------------------------Program-------------------------

x_conti=np.linspace(x_start,x_end,(x_end-x_start)/0.01)

for i in range(len(function_formula)):
    y_list=y_func(x_conti,function_formula[i])
    plt.plot(x_conti,y_list,color=color[i],linestyle="-",label=label[i])

plt.legend()
plt.xlabel(x_label,**font)
plt.ylabel(y_label,**font)
plt.ylim(y_limit[0],y_limit[1])
plt.xlim(x_start,x_end)
plt.tight_layout()
plt.savefig(title)
plt.show()




