#Made by Andreas L. Vishart
#Give the script an input (.com) or output (.out/.log) file with fragments in it.
#Then the vector and distace between the donor and the acceptor will be printed.
#------------------------Packages------------------------
import sys
import numpy as np
import argparse


#------------------------Parameters----------------------
#Donor fragment
fd=1
#Acceptor fragment
fa=2
#Bridge/Chromophore
fb=3
#print ("y"/"n" or "yes"/"no")
print_crit="y"

#------------------------Parser--------------------------
parser = argparse.ArgumentParser()
parser.add_argument('Input', help="Input file is needed.", type=str)
parser.add_argument("-fa", "--Frag_A", help="The label of the fragment for the acceptor.", type=int)
parser.add_argument("-fd", "--Frag_D", help="The label of the fragment for the donor.", type=int)
parser.add_argument("-fc", "--Frag_C", help="The label of the fragment for the chromophore.", type=int)
args = parser.parse_args()
if args.Frag_A:
    FA = args.Frag_A
if args.Frag_D:
    FD = args.Frag_D
if args.Frag_C:
    FB = args.Frag_C


#------------------------Functions-----------------------
#Get the xyz coordinates for each fragment
def get_xyz(content,ifstatement):
    xyz=[]
    elements=[]
    for line in content:
        if ifstatement in line:
            con=list(filter(None,line.split(" ")))
            elements.append(str(con[0].split("(")[0]))
            xyz.append([])
            for j in range(1,4):
                xyz[-1].append(float(con[j]))
    return elements,np.array(xyz)

#Distance between two elements
def distance(coord1,coord2):
    distance=0
    for i in range(3):
        distance+=(coord2[i]-coord1[i])**2
    return np.sqrt(distance)

def run_program(fd,fa,fb,print_crit):
    #Import input file
    with open(sys.argv[1]) as thefile:
        content=thefile.readlines()

    #Set the coordinates into lists of fragments
    Elements1,Frag_1=get_xyz(content,"Fragment=1")
    Elements2,Frag_2=get_xyz(content,"Fragment=2")
    Elements3,Frag_3=get_xyz(content,"Fragment=3")

    #Identify the fragments of the donor and acceptor
    Frag=[Frag_1,Frag_2,Frag_3]
    Frag_D=np.array(Frag[fd-1])
    Frag_A=np.array(Frag[fa-1])

    #Number of elements in the donor and acceptor fragment
    len_Frag_D=len(Frag_D)
    len_Frag_A=len(Frag_A)

    #Find the center of the donor and acceptor fragment from the distances alone
    center_Frag_D=np.array([sum(Frag_D[:,j]) for j in range(3)])/len_Frag_D
    center_Frag_A=np.array([sum(Frag_A[:,j]) for j in range(3)])/len_Frag_A

    #Maximum Radii of the donor fragment and acceptor
    radii_D=[]
    for j in range(len_Frag_D):
        radii_D.append(distance(center_Frag_D,Frag_D[j]))
    radii_D=max(radii_D)

    radii_A=[]
    for j in range(len_Frag_D):
        radii_A.append(distance(center_Frag_A,Frag_A[j]))
    radii_A=max(radii_A)

    #Vector from donor to acceptor
    DA_vec=center_Frag_A-center_Frag_D
    #Distace between donor and acceptor
    R_DA=np.sqrt(sum(DA_vec**2))

    if print_crit.lower()=="y" or print_crit.lower()=="yes":
        print("Radius of the donor= "+str(radii_D))
        print("Radius of the acceptor= "+str(radii_A))
        print("Vector between the donor and acceptor= "+str(DA_vec))
        print("Distance between the donor and acceptor= "+str(R_DA))
    
    return radii_D,radii_A,DA_vec,R_DA

#------------------------Program-------------------------

if __name__ == "__main__":
    run_program(fd,fa,fb,print_crit)