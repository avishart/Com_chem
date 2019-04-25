#Made by Andreas L. Vishart
#Give it a .csv file with labels in first row and x data in first column and y data on the other columns
#------------------------Packages------------------------
import numpy as np
import matplotlib.pyplot as plt
import sys


#------------------------Parameters----------------------
#Seperation symbol
sep=";"
#The name of saved plot
title="EC_PET3.png"
#X-axis title 
x_label="Distance / [Angstrom]"
#Y-axis title
y_label="EC / [mHartree]"
#List of colors for plotting
color=["red","blue","green","orange","pink","purple","black","yellow"]
#List of markers for plotting
markers=["o","o","o","o","o","s","s","s"]
#Want to Fit data (y/n)
fit_data="n"
#Function to fit to
func_fit="coeff[1]*np.exp(coeff[0]*x)"
#X label is a string (y/n)
x_str="y"

#------------------------Functions-----------------------
#Remove all none values
def rm_none_values(x_list,y_list):
    x_list_new=[]
    y_list_new=[]
    for i in range(len(x_list)):
        if y_list[i]==None or x_list[i]==None:
            continue
        else:
            y_list_new.append(y_list[i])
            x_list_new.append(x_list[i])
    return x_list_new,y_list_new

#Fit data to an exponential function
def fit_exp(x_list,y_list):
    y_list_new=[np.log(y) for y in y_list]
    coef=np.polyfit(x_list,y_list_new,1)
    coef[1]=np.exp(coef[1])
    return coef

#Calculate the values of an function
def y_func(x_list,funct,coeff):
    y_list=[eval(funct) for x in x_list]
    return y_list

#Calculate the R2 value of the fit
def R_squared(y_list,y_ref):
    len_y_ref=len(y_ref)
    mean_ref=sum(y_ref)/len_y_ref
    S_tot=sum([(y_ref[i]-mean_ref)**2 for i in range(len_y_ref)])
    S_res=sum([(y_list[i]-y_ref[i])**2 for i in range(len_y_ref)])
    R2=1-S_res/S_tot
    return R2



#------------------------Program-------------------------
#Load the .csv
with open(sys.argv[1]) as thefile:
    content=thefile.readlines()

#Get the labels of y values and remove them without
label=content[0].replace("\n","").split(sep)[1:]
len_list=len(label)
content=content[1:]
x_list=[]
y_list=[[] for i in range(len_list)]
#Get the data to x and y values
for line in content:
    newline=line.split(sep)
    x_list.append(newline[0])
    for i in range(1,len_list+1):
        try:
            y_list[i-1].append(float(newline[i]))
        except:
            y_list[i-1].append(None)
#If strings are used for x values then make a new x list
if x_str.lower()=="y" or x_str.lower()=="yes":
    x_num_list=[i+1 for i in range(len(x_list))]
    plt.xticks(x_num_list,x_list,fontsize=6,rotation=-90)
#If values are used as x values then make them floats
else:
    x_list=map(float,x_list)
    x_conti=np.linspace(x_list[0],x_list[-1],(x_list[-1]-x_list[0])/0.1)

#Plot data
for i in range(len_list):
    if x_str.lower()=="y" or x_str.lower()=="yes":
        x_list_new,y_list_new=rm_none_values(x_num_list,y_list[i])
        plt.plot(x_list_new,y_list_new,color=color[i],marker=markers[i],alpha=0.5,linestyle="",label=label[i])
    else:
        x_list_new,y_list_new=rm_none_values(x_list,y_list[i])
        plt.plot(x_list_new,y_list_new,color=color[i],marker=markers[i],linestyle="",label=label[i])
        #Fit only if x data are values and required
        if fit_data.lower()=="y" or fit_data.lower()=="yes":
            coeff=fit_exp(x_list_new,y_list_new)
            R2=R_squared(y_func(x_list_new,func_fit,coeff),y_list_new)
            plt.plot(x_conti,y_func(x_conti,func_fit,coeff),color=color[i],linestyle="-",label=str(round(coeff[1],1))+"*exp("+str(round(coeff[0],3))+"*x),$R^2$="+str(round(R2,3)))

#Show the plot
plt.legend(loc=0,fontsize=8)
plt.yscale('log')
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.tight_layout()
plt.savefig(title)
plt.show()
plt.close()

