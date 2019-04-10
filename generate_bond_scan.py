#Made by Andreas Vishart
#Use one .com file with single point calculation
#The two atoms in the bond distance is labeled and a fragment 1 and 2 is denoted
#The two fragments are fixed compared to each self and only the distance in focus

#------------Package---------------------------
import numpy as np
import sys
import copy as cp

#------------Parameters---------------------------
#The first atom in the bond
atom1=2
#The second atom in the bond
atom2=4
#Smallest bond length 
bond_from=1.0
#Highest bond length 
bond_to=3.0
#Number of bond length in between
steps=20
bonds_list=np.linspace(bond_from,bond_to,steps)
#print(bonds_list)
#End
end=".com"  #.com or .xyz
#Keep fragments in the new files?
frag="n" #y/n


#------------Program------------------------------
#Download file
try:
    with open(sys.argv[1]) as thefile:
        content=thefile.readlines()
except:
    sys.exit("No file given!")

#Extract xyz coordinates
xyz_data=[]
start=0
for line in content:
    if "Fragment" in line:
        xyz_data.append(list(filter(None,line.split(" "))))
        if start==0:
            start=content.index(line)

length=len(xyz_data)

#The initial vector for the bond
dis_vec=[float(xyz_data[atom2-1][i])-float(xyz_data[atom1-1][i]) for i in range(1,4)]
dis_vec=np.array(dis_vec)
#The initial length of the bond
dis=np.linalg.norm(dis_vec,2)


for leng in bonds_list:
    #Copy of xyz list, so it is NOT overwritten
    xyz_data_file=cp.deepcopy(xyz_data)
    #Scaling the bond length
    scale=(leng/dis)
    vec=(scale-1)*dis_vec
    #The new files are created
    new_file=open(sys.argv[1][:-4]+"_"+str("{:.3f}".format(round(leng,3)))+end,'w')
    if end==".com":
        #Head part of .com file
        for i in range(0,start):
            new_file.write(content[i])
    elif end==".xyz":
        #xyz file:
        new_file.write(str(length)+"\n\n")
    #New XYZ coordinates
    for i in range(length):
        if "Fragment=1" in xyz_data_file[i][0]:
            xyz_data_file[i][3]=xyz_data_file[i][3][:-2]
            for j in range(1,4):
                xyz_data_file[i][j]="{:.8f}".format(round(float(xyz_data_file[i][j]),8))
                if float(xyz_data_file[i][j])>=0:
                    xyz_data_file[i][j]=" "+xyz_data_file[i][j]
        elif "Fragment=2" in xyz_data_file[i][0]:
            for j in range(1,4):
                value=float(xyz_data_file[i][j])+vec[j-1]
                xyz_data_file[i][j]="{:.8f}".format(round(value,8))
                if value>=0:
                    xyz_data_file[i][j]=" "+xyz_data_file[i][j]
        #Remove "Fragment=" if frag==n
        if frag=="n":
            try:
                indi=xyz_data_file[i][0].index("(")
                if indi:
                    xyz_data_file[i][0]=xyz_data_file[i][0][:indi]
            except:
                pass
        data_new=" "+xyz_data_file[i][0]+"     "+xyz_data_file[i][1]+"   "+xyz_data_file[i][2]+"   "+xyz_data_file[i][3]+"\n"
        new_file.write(data_new)
    #Tail part of .com file
    if end==".com": 
        for i in range(start+length,len(content)):
            new_file.write(content[i])
    new_file.close()





