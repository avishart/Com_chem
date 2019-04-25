#Made by Andreas L. Vishart
#Give it a .csv file with labels in first row and x data in first column and y data on the other columns
#------------------------Packages------------------------
import numpy as np
import matplotlib.pyplot as plt
import sys


#------------------------Parameters----------------------
#Seperation symbol
sep=","
#The name of saved plot
title="EC_DBA.png"
#X-axis title 
x_label="Distance / [Angstrom]"
#Y-axis title
y_label="EC / [mHartree]"
#List of colors for plotting
color=["red","blue","green","orange","pink","purple","black","yellow"]
#Want to Fit data (y/n)
fit_data="n"

#------------------------Functions-----------------------
def rm_none_values(x_list,y_list):
    x_list_new=[]
    y_list_new=[]
    for i in range(len(y_list)):
        if y_list[i]==None or x_list[i]==None:
            continue
        else:
            y_list_new.append(y_list[i])
            x_list_new.append(x_list[i])
    return x_list_new,y_list_new

def fit_exp(x_list,y_list):
    y_list_new=[np.log(y) for y in y_list]
    coef=np.polyfit(x_list,y_list_new,1)
    coef[1]=np.exp(coef[1])
    return coef

def y_func(x_list,funct,coeff):
    y_list=[eval(funct) for x in x_list]
    return y_list

def R_squared(y_list,y_ref):
    len_y_ref=len(y_ref)
    mean_ref=sum(y_ref)/len_y_ref
    S_tot=sum([(y_ref[i]-mean_ref)**2 for i in range(len_y_ref)])
    S_res=sum([(y_list[i]-y_ref[i])**2 for i in range(len_y_ref)])
    R2=1-S_res/S_tot
    return R2



#------------------------Program-------------------------
with open(sys.argv[1]) as thefile:
    content=thefile.readlines()

label=list(filter(None,content[0].split(sep)))[1:]
content=content[1:]
x_list=[]
len_list=len(list(filter(None,content[0].split(sep))))
y_list=[[] for i in range(len_list-1)]
for line in content:
    newline=list(filter(None,line.split(sep)))
    x_list.append(float(newline[0]))
    for i in range(1,len_list):
        try:
            y_list[i-1].append(float(newline[i]))
        except:
            y_list[i-1].append(None)

x_conti=np.linspace(x_list[0],x_list[-1],(x_list[-1]-x_list[0])/0.1)

for i in range(len_list-1):
    x_list_new,y_list_new=rm_none_values(x_list,y_list[i])
    plt.plot(x_list_new,y_list_new,color=color[i],marker="o",linestyle="",label=label[i])
    if fit_data.lower()=="y" or fit_data.lower()=="yes":
        coeff=fit_exp(x_list_new,y_list_new)
        R2=R_squared(y_func(x_list_new,"coeff[1]*np.exp(coeff[0]*x)",coeff),y_list_new)
        plt.plot(x_conti,y_func(x_conti,"coeff[1]*np.exp(coeff[0]*x)",coeff),color=color[i],linestyle="-",label=str(round(coeff[1],1))+"*exp("+str(round(coeff[0],3))+"*x),$R^2$="+str(round(R2,3)))

plt.legend()
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.savefig(title)
plt.show()




