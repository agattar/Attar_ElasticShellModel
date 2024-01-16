# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 12:14:50 2022

@author: Goktug
"""
import math 
import random
import numpy as np

# def Crosslinker(Chr_Num, prob,minb,maxb):

file = ("data.minimized")

file_open = open(file,"r")
lines = file_open.readlines()
file_open.close()


bond_num = int(lines[4][:6])
atom_num = int(lines[2][:5]) 

for k,i in enumerate(lines):
    if 'Bonds' in i:
        bondline = k+2
        
for k,i in enumerate(lines):
    if 'Atoms' in i:
        atomline = k+2 

l_list = []
for i in lines[atomline:atomline+atom_num]:
    i = i.split(' ')
    if i[1] == "1":
        l_list.append(i[:6])

b_list = np.zeros((len(l_list)*7,4))
allbonds = []
count = 0
count_1 = 0

for k in (lines[bondline:bondline+bond_num]):
    k = k.split()
    if k[1] == '1':
        b_list[count,0] = k[0]
        b_list[count,1] = k[1]
        b_list[count,2] = k[2]
        b_list[count,3] = k[3]
        count += 1
    else:
        allbonds.append(str(k[0])+"\t"+str(k[1])+"\t"+str(k[2])+"\t"+str(k[3])+"\n")

del (lines[bondline:bondline+bond_num])
       
ind = bond_num+1

ind = 1

print("Cap Bonds are being generated...")

minbond=5; maxbond=6;factor=1.5


bonds = np.zeros((len(l_list)*4,2))        
bondnumb = np.random.randint(minbond,maxbond,size=(len(l_list),1))

bondnumbers = np.linspace(1,len(l_list),len(l_list))
bondnumbers = np.reshape(bondnumbers,(len(l_list),1))

bondnumbers = np.hstack((bondnumbers,bondnumb))

l_list = np.array(l_list)

l_list = l_list.astype(float)

l_list = l_list[l_list[:, 0].argsort()]

l_list = l_list[:,3:]

bead_num = np.linspace(1,len(l_list),len(l_list))
bead_num = bead_num.astype(int);bead_num = bead_num.astype(str);bead_num=np.reshape(bead_num,(len(l_list),1))
l_list = np.column_stack((bead_num,l_list))

bondarray = l_list[:,:].astype(float)
       
bnum = []
for n in range(len(bondarray)):
    r_inc = 0
    neighbor_list = []
    Extracted = np.subtract(bondarray[n,1:],bondarray[:,1:])
    Multiplied = np.multiply(Extracted,Extracted)
    Multiplied = np.sum(Multiplied,axis=1)  
    Distances = np.sqrt(Multiplied)
    
    Result = np.where((Distances[:] < 0.80*factor) & (Distances[:] >= 0.60*factor))
    
    selected_num = int(bondnumb[n]) - int(bondnumbers[n][1])
                
    for k in Result[0]:
        if selected_num < 8 and ([k+1,n+1]==bonds[:,None]).all(-1).any() == False:
            neighbor_list.append(k)
                
                    
    while (int(bondnumbers[n,1])) !=0:
            
        selected_num = int(bondnumb[n]) - int(bondnumbers[n][1])

        if len(neighbor_list) != 0:
            chosen = random.choice(neighbor_list)
            if ([chosen+1,n+1]==bonds[:,None]).all(-1).any() == False:
                bondnumbers[int(chosen),1] -= 1
                bondnumbers[int(n),1] -= 1
                bnum.append(str(ind)+"\t"+"1"+"\t"+str(n+1)+"\t"+str(chosen+1)+"\n")
                bonds[ind-1,0] = str(n+1); bonds[ind-1,1] = str(chosen+1)
                ind += 1
                neighbor_list.remove(chosen)              
        
        elif r_inc > 0.50*factor:
            break
        
        else:
            
            neighbor_list_2 = []

            Result = np.where((Distances[:] >= 0.80*factor+r_inc) & (Distances[:] < 0.90*factor+r_inc))
        
            for k in Result[0]:
                selected_num_n = int(bondnumb[k]) - int(bondnumbers[k][1])
                if selected_num < 8 and selected_num_n < 8  and ([k+1,n+1]==bonds[:,None]).all(-1).any() == False and ([n+1,k+1]==bonds[:,None]).all(-1).any() == False:
                    neighbor_list_2.append(k)
                                
            
            if len(neighbor_list_2) != 0:
                chosen = random.choice(neighbor_list_2)
                
                if ([chosen+1,n+1]==bonds[:,None]).all(-1).any() == False and ([n+1,chosen+1]==bonds[:,None]).all(-1).any() == False:
                    bondnumbers[int(chosen),1] -= 1
                    bondnumbers[int(n),1] -= 1
                    bnum.append(str(ind)+"\t"+"1"+"\t"+str(n+1)+"\t"+str(chosen+1)+"\n")
                    bonds[ind-1,0] = str(n+1); bonds[ind-1,1] = str(chosen+1)
                    ind += 1
                    neighbor_list_2.remove(chosen)     
                    
            else:
                r_inc += 0.1
                  
mean_bondnum_ = np.subtract(bondnumb[:,0],bondnumbers[:,1])   
mean_bondnum = np.mean(mean_bondnum_)
print(mean_bondnum,"Lamin Mean Bond Number")


for t in bnum:
    allbonds.append(t)

bnum = []
for k,s in enumerate(allbonds):
    s = s.split('\t')
    s[0] = str(k+1)
    bnum.append(str(s[0])+"\t"+str(s[1])+"\t"+str(s[2])+"\t"+str(s[3]))

indstr = str(str(len(bnum))+" "+"bonds"+"\n")
bndstr = str(str(3)+" "+"bond types"+"\n")

lines[4] = indstr
lines[8] = bndstr
    
bnum.append("\n")
lines.insert(bondline,bnum)

fl = open("data.minimized",'w+')
     
for line in lines:
    fl.write("".join(line))
    
fl.close()
