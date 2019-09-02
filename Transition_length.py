#Made by Andreas L. Vishart
#Give the script an input (.com) or output (.out/.log) file with fragments in it.
#Then the vector and distace between the donor and the acceptor will be printed.
#------------------------Packages------------------------
import numpy as np
import argparse


#------------------------Parameters----------------------
#Donor fragment
fd=1
#Acceptor fragment
fa=2
#Bridge/Chromophore
fb=0
#print ("y"/"n" or "yes"/"no")
print_crit="y"
#VdW radii
vdw_radi=1.5
#Convergence factor, Angstrom to Bohr
conv_dist=1/(0.529177210903)

#------------------------Parser--------------------------
parser = argparse.ArgumentParser()
parser.add_argument('Input', help="Input file is needed.", type=str)
parser.add_argument("-fa", "--Frag_A", help="The label of the fragment for the acceptor.", type=int)
parser.add_argument("-fd", "--Frag_D", help="The label of the fragment for the donor.", type=int)
parser.add_argument("-fc", "--Frag_C", help="The label of the fragment for the chromophore.", type=int)
args = parser.parse_args()
if args.Frag_A:
    fa = args.Frag_A
if args.Frag_D:
    fd = args.Frag_D
if args.Frag_C:
    fb = args.Frag_C


#------------------------Functions-----------------------
#Atom info
def atom_info():
    Elements={}
    Elements["Element"]=["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
                            "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                            "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"]
    Elements["Number"]=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 
                            25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 
                            47, 48, 49, 50, 51, 52, 53, 54]
    Elements["VdW"]=[1.20,1.40,1.82,1.53,1.92,1.70,1.55,1.52,1.47,1.54,2.27,1.73,1.84,2.10,1.80,1.80,1.75,1.88,
                        2.75,2.31,None,None,None,None,None,None,None,1.63,1.40,1.39,1.87,2.11,1.85,1.90,1.85,2.02,
                        3.03,2.49,None,None,None,None,None,None,None,1.63,1.72,1.58,1.93,2.17,2.06,2.06,1.98,2.16]
    return Elements

#Get the xyz coordinates for each fragment of an input file
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

#Get the xyz coordinates for each fragment of an output file
def load_xyz(out_content):
    frag=False
    frag_content=[]
    out=False
    for line in out_content:
        #Get fragments for each element
        if "Charge" in line and "Multiplicity" in line:
            frag=True
        elif "Initial Parameters" in line or "Stoichiometry" in line:
            frag=False
            break
        if frag:
            newline=list(filter(None,line.replace("\n","").split(" ")))
            if len(newline)==4:
                frag_content.append(newline[0])
    for line in out_content:
        #Get the final coordinates
        if "Standard orientation:" in line:
            out=True
            l=1
            out_xyz=[]
        elif "Rotational constants" in line:
            out=False
        elif out:
            newline=list(filter(None,line.replace("\n","").split(" ")))
            if newline[0]==str(l):
                out_xyz.append(["{0:.8f}".format(round(float(newline[3]),8)),"{0:.8f}".format(round(float(newline[4]),8)),"{0:.8f}".format(round(float(newline[5]),8))])
                l+=1
    #If elements and coordinates are loaded coorectly then make the xyz input
    xyz_list=[]
    if len(frag_content)==len(out_xyz):
        for i in range(len(frag_content)):
            xyz_txt=""
            space=32-len(str(frag_content[i]))-len(str(out_xyz[i][0]))
            xyz_txt+=" "+str(frag_content[i])+" "*space
            for j in range(2):
                space=16-len(str(out_xyz[i][j+1]))
                xyz_txt+=str(out_xyz[i][j])+" "*space
            xyz_txt+=str(out_xyz[i][2])
            xyz_list.append(xyz_txt)
    return xyz_list

#Distance between two elements
def distance(coord1,coord2):
    distance=0
    for i in range(3):
        distance+=(coord2[i]-coord1[i])**2
    return np.sqrt(distance)

#Get the center of charge and the radii of the fragment
def center_radii_fragment(content,frag_num,Elements_prop,conv_dist):
    #Set the coordinates into lists of fragments
    Elements_frag,Frag_frag=get_xyz(content,"Fragment="+str(frag_num))
    #Number of elements in the fragment
    len_Frag=len(Frag_frag)
    #Coordinates of elements by charge weighted factor
    weight_frag=np.array([Elements_prop["Number"][Elements_prop["Element"].index(i)] for i in Elements_frag])
    Frag_frag_weight=np.array([Frag_frag[i]*weight_frag[i] for i in range(len_Frag)])
    #Find the center of the donor and acceptor fragment from the distances alone
    center_Frag_frag=np.array([sum(Frag_frag_weight[:,j]) for j in range(3)])/(sum(weight_frag))
    #Maximum Radii of the fragment with and without VdW radii
    radii_frag=0
    for j in range(len_Frag):
        dis=distance(center_Frag_frag,Frag_frag[j])
        if dis>=radii_frag:
            radii_frag=dis
            vdw=Elements_prop["VdW"][Elements_prop["Element"].index(Elements_frag[j])]
    radii_frag=radii_frag*conv_dist
    radii_frag_vdw=radii_frag+(vdw*conv_dist)
    return Frag_frag,center_Frag_frag,radii_frag,radii_frag_vdw


def run_program(fd,fa,fb,print_crit,vdw_radi):
    #Import input file
    with open(args.Input) as thefile:
        content=thefile.readlines()
    #If the file is an output file else it is an input file
    if args.Input[-3:]=="out" or args.Input[-3:]=="log":
        content=load_xyz(content)
    Elements_prop=atom_info()

    #Center and radii of donor
    Frag_D,center_Frag_D,radii_D,radii_D_vdw=center_radii_fragment(content,fd,Elements_prop,conv_dist)
    #Center and radii of acceptor
    Frag_A,center_Frag_A,radii_A,radii_A_vdw=center_radii_fragment(content,fa,Elements_prop,conv_dist)

    #Vector from donor to acceptor
    DA_vec=(center_Frag_A-center_Frag_D)*conv_dist
    #Distace between donor and acceptor
    R_DA=np.sqrt(sum(DA_vec**2))

    #All transition distances
    trans_dis=[]
    for i in range(len(Frag_D)):
        for j in range(len(Frag_A)):
            trans_dis.append(distance(Frag_A[j],Frag_D[i]))
    R_DA_min=min(trans_dis)*conv_dist

    #Center, radii, vector of bridge
    if fb!=0:
        try:
            Frag_B,center_Frag_B,radii_B,radii_B_vdw=center_radii_fragment(content,fb,Elements_prop,conv_dist)
            if radii_B>0:
                DB_vec=(center_Frag_B-center_Frag_D)*conv_dist
                R_DB=np.sqrt(sum(DB_vec**2))
                BA_vec=(center_Frag_A-center_Frag_B)*conv_dist
                R_BA=np.sqrt(sum(BA_vec**2))
                trans_dis_DB=[]; trans_dis_BA=[]
                for i in range(len(Frag_B)):
                    for j in range(len(Frag_D)):
                        trans_dis_DB.append(distance(Frag_B[i],Frag_D[j]))
                    for j in range(len(Frag_A)):
                        trans_dis_BA.append(distance(Frag_B[i],Frag_A[j]))
                R_DB_min=min(trans_dis_DB)*conv_dist
                R_BA_min=min(trans_dis_BA)*conv_dist
        except:
            fb=0

    if print_crit.lower()=="y" or print_crit.lower()=="yes":
        print("Radius of the donor without VdW = "+str(radii_D)+" Bohr")
        print("Radius of the acceptor without VdW = "+str(radii_A)+" Bohr")
        print("Radius of the donor with VdW = "+str(radii_D_vdw)+" Bohr")
        print("Radius of the acceptor with VdW = "+str(radii_A_vdw)+" Bohr")
        print("Vector between the donor and acceptor centers = "+str(DA_vec)+" Bohr")
        print("Minimum distance between the donor and acceptor = "+str(R_DA_min)+" Bohr")
        print("Distance between the donor and acceptor centers= "+str(R_DA)+" Bohr")
        if fb!=0 and radii_B>0:
            print("Radius of the bridge without VdW = "+str(radii_B)+" Bohr")
            print("Radius of the bridge with VdW = "+str(radii_B_vdw)+" Bohr")
            print("Minimum distance between the donor and bridge = "+str(R_DB_min)+" Bohr")
            print("Distance between the donor and bridge centers= "+str(R_DB)+" Bohr")
            print("Minimum distance between the bridge and acceptor = "+str(R_BA_min)+" Bohr")
            print("Distance between the bridge and acceptor centers= "+str(R_BA)+" Bohr")
            return radii_D,radii_A,DA_vec,R_DA,R_DA_min,radii_B,R_DB,R_BA
        else:
            return radii_D,radii_A,DA_vec,R_DA,R_DA_min

#------------------------Program-------------------------

if __name__ == "__main__":
    run_program(fd,fa,fb,print_crit,vdw_radi)