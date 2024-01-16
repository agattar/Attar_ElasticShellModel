# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 01:00:51 2021

@author: Asus
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 00:30:58 2021

@author: Asus
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 22:17:06 2021

@author: Asus
"""
import timeit
    

def Nucleus_Generate(beads,chnum,place_factor,minbond,maxbond,factor,mut_ratio,crosslinker=False,mutation=False):
    
    import math
    import numpy as np
    import random
    
    start = timeit.default_timer()
    def Cap(factor,beads):
        global r
        global Sd
        global cap; global cap1
        global index
        
        cappoint = np.array([])
        r = 28
        r = r*factor
        broke = False

        # beads = int(beads/factor)
        additional = 0
        
        info = np.arange(1,beads+1+additional,1,dtype=int)
        info = info.reshape(beads+additional,1)
        typelist = np.ones((beads+additional,1))
        caplist = np.hstack((info,typelist))
        caplist = np.hstack((caplist,typelist))
        caplist = np.char.mod('%d', caplist)

        index = 0
        while not broke:
            xpos = np.random.uniform(-r,r,100)
            ypos = np.random.uniform(-r,r,100)
            zpos = np.random.uniform(-r,r,100)
            X2 = np.multiply(xpos,xpos)
            Y2 = np.multiply(ypos,ypos)
            Z2 = np.multiply(zpos,zpos)
                
            Sd = X2 +  Y2[:, np.newaxis] 
            Sd = Sd +  Z2[:, np.newaxis]

            wb = np.where(np.logical_and(np.greater(Sd,r**2-5),np.less(Sd,r**2+5)))
            if len(wb) == 0:
                break
            wb = np.array(wb); wb = np.swapaxes(wb,0,1)

            if len(cappoint) < beads*3:
                selected = np.zeros([1,3])
                for i in range(len(wb)):
                    x = random.choice(wb)       
                    if index > 0 and len(cappoint) < beads*3:
                        cappoint = cappoint.reshape(index,3)
                        selected[0,0] = xpos[x[1]]; selected[0,1] = ypos[x[0]]; selected[0,2] = zpos[x[0]]
                        extracted = np.subtract(cappoint[:,:],selected)      
                        multiplied = np.multiply(extracted,extracted)
                        summ = np.sum(multiplied,axis=1)
                        dist = np.sqrt(summ)
                        if np.any(dist[:] < 0.66*factor) == True and np.any(dist[:] > 1.0*factor) == True:
                            break 
                        else:
                            index +=1
                            cappoint = np.append(cappoint,[xpos[x[1]],ypos[x[0]],zpos[x[0]]])
                            cappoint = cappoint.astype('float')
                            
                    elif index == 0:
                        index += 1
                        cappoint = np.append(cappoint,[xpos[x[1]],ypos[x[0]],zpos[x[0]]])
                        cappoint = cappoint.astype('float')
                        
            elif len(cappoint) >= beads*3 and len(cappoint) < beads*3+additional*3:

                selected = np.zeros([1,3])
                for i in range(len(wb)):
                    x = random.choice(wb)       
                    cappoint = cappoint.reshape(index,3)
                    selected[0,0] = xpos[x[1]]; selected[0,1] = ypos[x[0]]; selected[0,2] = zpos[x[0]]
                    extracted = np.subtract(cappoint[:,:],selected)      
                    multiplied = np.multiply(extracted,extracted)
                    summ = np.sum(multiplied,axis=1)
                    dist = np.sqrt(summ)
                    if np.any(dist[:] < 0.6*factor) == True and np.any(dist[:] > 0.8*factor) == True:
                        break 
                    else:
                        index +=1
                        cappoint = np.append(cappoint,[xpos[x[1]],ypos[x[0]],zpos[x[0]]])
                        cappoint = cappoint.astype('float')
            else:
                broke = True
                
                
        timediff = timeit.default_timer()-start    
        print(str(int(timediff/60))+" "+str("minutes") + " " + str(float('{:.2f}'.format(timediff))%60)+" " + str("seconds"))   
                    
        cappoint = cappoint.reshape(index,3)
        cap = np.hstack((caplist,cappoint))
        cap1 = 0
        
    Cap(factor,beads)
    
    def Chromosome(chnum,index,place_factor):
        file = np.load("mm9_domains_200k.npy") 
        file_array = np.array([])
        for i in file:
            file_array = np.append(file_array,i)
        
        global dnaloc
        global xcap
        global DNAl
        global linker
        global index1
        
        dnaloc = np.array([])
            
        typ1 = "2"
        typ2 = "3"
        typ3 = "4"
        typ4 = "5"
        typ5 = "6"
        typ6 = "7"    
        typ7 = "8"  
        
        dna = np.array([])
        
        index += 1
        
        xcap = cap
    
         
        turn = 20
        block = 200
        DNAl = 6000
        chcount = 0
        
        xcap = xcap.astype('float')
        xcap[:,3:] = np.true_divide(xcap[:,3:],place_factor)
        xcap = xcap.astype('str')

        x1 = random.choice(xcap)
        
        chtry = 0
        counter = 0
        gap = 5000
        
        xc = float(x1[3])
        yc = float(x1[4])
        zc = float(x1[5])
        
        angledotlist = []
        
        linker = 1
        
        while chcount<int(chnum):
                
                index = len(cap)+ 1 + chcount*6002
                
                x1 = random.choice(xcap)          
                        
                xc = float(x1[3])
                yc = float(x1[4])
                zc = float(x1[5])
        
                if chcount >= int(chnum/2):
                    counter = 5000
                else:
                    counter = 0
                
                anglefound = False
                
                while anglefound == False:                
                    for i in xcap:
                        zsum = math.sqrt((float(xc-float(i[3]))**2)+(float(yc-float(i[4]))**2)+(float(zc-float(i[5]))**2))
                        if zsum >= 30:
                            angledotlist.append((float(i[3]),float(i[4]),float(i[5])))        
                            anglefound = True
                                
                angledot = (random.choice(angledotlist))
                                
                xinc = abs((abs(xc)+abs(angledot[0]))/(2*r))+0.02
                yinc = abs((abs(yc)+abs(angledot[1]))/(2*r))+0.02
                zinc = abs((abs(zc)+abs(angledot[2]))/(2*r))+0.02
                    
                
                tup = (xinc,yinc,zinc)
                
                dirnum=(max(tup))
                
                
                cdir = -1
                for i in tup:
                    cdir +=1
                    if abs(i) == dirnum:
                        break
                
                if cdir == 0:
                    direction = ("xinc")
                elif cdir == 1:
                    direction = ("yinc")
                else:
                    direction = ("zinc")
                
                
                if angledot[0] < xc:
                    xinc = -xinc
                if angledot[1] < yc:
                    yinc = -yinc
                if angledot[2] < zc:
                    zinc = -zinc
                    
                
                if yc<0:
                    yc +=1.6
                    ycinc = 0.9
                else:
                    yc -=1.6
                    ycinc =-0.9
                    
                if xc<0:
                    xc +=1.6
                    xcinc = 0.9
                else:
                    xc -=1.6
                    xcinc =-0.9
                    
                if zc<0:
                    zc +=1.6
                    zcinc = 0.9
                else:
                    zc -=1.6
                    zcinc =-0.9
                    
                linker = 1
                for i in range(linker):
                    if i == 0:
                        dna= np.append(dna,[index,typ4, typ4,xc,yc,zc])
                        index +=1         
                        xc += xinc
                        yc += yinc
                        zc += zinc
                
            
                for i in range(DNAl+linker):
                        
                    if i % turn == 0 and direction == "xinc" and i!=0 and i!=DNAl+linker*2-2:
                        if i % block == 0 and i!=0:
                              ycinc = -ycinc
                              zc += zcinc + int(zcinc/10)
                              yc -= ycinc
                            
                        yc += ycinc  
                        
                    elif i % turn == 0 and direction == "yinc" and i!=0 and i!=DNAl+linker*2-2:
                        if i % block == 0 and i!=0:
                            xcinc = -xcinc
                            zc += zcinc + int(zcinc/10)
                            xc -= xcinc
                            
                        xc += xcinc  
                        
                    elif i % turn == 0 and direction == "zinc" and i!=0 and i!=DNAl+linker*2-2:
                        if i % block == 0 and i!=0:
                            ycinc = -ycinc
                            xc += xcinc + int(xcinc/10)
                            yc -= ycinc
                                           
                        yc += ycinc  
                    
                    if i != 0 and i % turn ==0:
                        yinc = -yinc  
                        zinc = -zinc
                        xinc = -xinc
                            
                
                        xc += xinc
                        zc += zinc
                        yc += yinc
                    
                    
                    if i < 5:
                        dna= np.append(dna,[index,typ5, typ5,xc,yc,zc])
                        zc += zinc
                        yc += yinc
                        xc += xinc
                        index +=1 
                        
                    elif i >= 5 and i < 1000:
                        dna= np.append(dna,[index,typ3, typ3,xc,yc,zc])
                        zc += zinc
                        yc += yinc
                        xc += xinc
                        index +=1 
                        
                    elif i>=1000 and i<DNAl:
                        point = file_array[counter:counter+1]
                        point = float(point)
                        counter += 1
                        if point == 0:
                            dna= np.append(dna,[index,typ1, typ1,xc,yc,zc])
                            zc += zinc
                            yc += yinc
                            xc += xinc
                            index +=1 
                        
                        elif point == 1:
                            dna= np.append(dna,[index,typ2, typ2,xc,yc,zc])
                            zc += zinc
                            yc += yinc
                            xc += xinc
                            index +=1                         
                        
                for i in range(linker):
                    if i == 0:
                        xc -= xinc
                        yc -= yinc
                        zc -= zinc
                        dna= np.append(dna,[index,typ4, typ4,xc,yc,zc])
                        index +=1 
                        
                dna = dna.reshape(6002,6)
                
                for i in dna[len(dna)-2:len(dna)-7:-1]:

                    if i[1] == "2":
                        i[1] = "8"
                        i[2] = "8"
                    else:
                        i[1] = "7"
                        i[2] = "7"
                
                if chcount==0:
                    broke = True
                    for i in dna:
                        raddis = math.sqrt(float(i[3])**2+float(i[4])**2+float(i[5])**2)

                        if raddis >= (abs(r)/place_factor):
                            broke = False
                            dna = np.array([])
                            chtry += 1
                            counter -= gap
                            break              
                
                    if  broke == True:
                        dna = dna.reshape(6002,6)
                        dnaloc = np.array(dna)
                        chcount +=1
                        chtry +=1
                        print("Found and Placed: Chromosome"+"\t"+str(chcount))
    
                        dna = np.array([])

                else:
                    broke = False
                    Raddists = dna[:,3:].astype(float)
                    Raddist_1 = np.multiply(Raddists[:,:],Raddists[:,:])
                    Raddist_1 = np.sum(Raddist_1,axis=1); Raddist_1 = np.sqrt(Raddist_1) 
                    
                    Raddists_all = dnaloc[:,3:].astype(float)
                    Raddist_2 = np.multiply(Raddists_all[:,:],Raddists_all[:,:])
                    Raddist_2 = np.sum(Raddist_2,axis=1); Raddist_2 = np.sqrt(Raddist_2) 
                
                    for k in range(6002):
                    
                        Extracted = np.subtract(Raddists_all[:,:],Raddists[k,:])
                        Multiplied = np.multiply(Extracted,Extracted)
                        Multiplied = np.sum(Multiplied,axis=1)  
                        Distances = np.sqrt(Multiplied)
                        if (np.any(Distances[:] < 1.5) == True or np.any(Raddist_1[:] > (abs(r)/place_factor)) == True or
                                    np.any(Raddist_2[:] > (abs(r)/place_factor)) == True):
                            broke = True
                            dna = np.array([])
                            chtry +=1
                            break
                        
                    if broke == True:
                        counter -= gap

                    if broke == False:
                        dnaloc = np.append(dnaloc,dna)
                        chcount +=1
                        chtry +=1
                        print("Found and Placed: Chromosome"+"\t"+str(chcount))
                        dnaloc = dnaloc.reshape((int(len(dnaloc)/6),6))
                        dna = np.array([])

                        
                       
        dnaloc = dnaloc.reshape(chcount*6002,6)
        index1 = index
        return (chnum,DNAl,linker)
    
    Chromosome(chnum, index,place_factor)
    
    global Cappoint
    global bnum

    bnum = []
    Cappoint = cap[:,3:].astype(float)

    
    def Bond_Generator(beads,cap,cap1,minbond,maxbond,factor):
        
        global bonds
        global bondnumbers
        global bondnumb
        global mean_bondnum
        global neighbor_list
        global chosen 
        global Distances
        global ind
        global missing_bonds
        
        ind = 1

        print("Cap Bonds are being generated...")
        
        bonds = np.zeros((beads*4,2))        
        bondnumb = np.random.randint(minbond,maxbond,size=(len(cap),1))

        bondnumbers = np.linspace(1,len(cap),len(cap))
        bondnumbers = np.reshape(bondnumbers,(len(cap),1))

        bondnumbers = np.hstack((bondnumbers,bondnumb))

        bondarray = cap[:,:].astype(float)
               

        for n in range(len(bondarray)):
            r_inc = 0
            neighbor_list = []
            Extracted = np.subtract(bondarray[n,3:],bondarray[:,3:])
            Multiplied = np.multiply(Extracted,Extracted)
            Multiplied = np.sum(Multiplied,axis=1)  
            Distances = np.sqrt(Multiplied)
            
            Result = np.where((Distances[:] < 0.9*factor) & (Distances[:] >= 0.66*factor))
            
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
                
                elif r_inc > 0.30*factor:
                    break
                
                else:
                    
                    neighbor_list_2 = []

                    Result = np.where((Distances[:] >= 0.9*factor+r_inc) & (Distances[:] < 1.0*factor+r_inc))
                
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
        
        missing_bonds = []
        for k,s in enumerate(mean_bondnum_):
            if s < 5:
                missing_bonds.append(k)
        
        
        count = 1
        for j in range(len(cap)+1,len(cap) + chnum*DNAl+1+chnum*linker*2):

            if count // (DNAl+2*linker) == 1:
                count = 1   
                pass
            
            elif j+1<len(cap) + (chnum*DNAl+1+chnum*linker*2):
                bnum.append(str(ind)+"\t"+"2"+"\t"+str(j)+"\t"+str(j+1)+"\n")
                ind +=1
                count +=1
                    
    Bond_Generator(beads,cap,cap1,minbond,maxbond,factor)
    
    print("Bonds are generated...")
    
    
    def Crosslinker(chnum, prob,minbond,maxbond):
        
        global b_list
        global Each_Chr
        global c_het
        global ind
        
        global bondnumb
        global bondarray
        global beadname
        global bondnumbers
        
        finished = True
        gap = 0
        Each_Chr = np.array([])
        c_het = np.array([])
        b_list = np.zeros((chnum*6001+10000,2))
        count = 0
        
        while finished != False: 
            for m in range(chnum):
                for k,i in enumerate(dnaloc[gap+1:gap+1001]):
                    Each_Chr = np.append(Each_Chr,[i[0],i[1],i[2],i[3],i[4],i[5]])
               
                if m==0:
                    c_het = np.array(Each_Chr)
                else:
                    c_het = np.append(c_het,Each_Chr)
                    
                Each_Chr = np.array([])
                gap += 6002
                
            finished = False
        
        c_het = np.reshape(c_het, (chnum*1000,6))
                    
        bondnumb = np.random.randint(minbond,maxbond,size=(len(c_het),1))
        bondnumbers = np.linspace(1,len(c_het),len(c_het))
        bondnumbers = np.reshape(bondnumbers,(len(c_het),1))
        bondnumbers = np.hstack((bondnumbers,bondnumb))
        

        
        gap = -1000
        b_index = 0
        
        for k in range(chnum):
            gap += 1000
            bondarray = c_het[gap:1000+gap,:].astype(float)
            beadname = c_het[gap:1000+gap,:]
            
            for n in range(len(bondarray[:,:])):
                
                r = random.random()
                if r <= prob:
                    neighbor_list = []
                    Extracted = np.subtract(bondarray[n,3:],bondarray[:,3:])
                    Multiplied = np.multiply(Extracted,Extracted)
                    Multiplied = np.sum(Multiplied,axis=1)  
                    Distances = np.sqrt(Multiplied)
                
                    Result = np.where((Distances[:] < 2.00) & (Distances[:] > 1.00))
                    
                    for k in Result[0]:
                        
                        selected_num = int(bondnumbers[k+gap][1])
                        
                        if ([int(beadname[k][0]),int(beadname[n][0])]==b_list[:,None]).all(-1).any() == False and selected_num > 0:
                            diff = k - n
                            if (diff == 1 or diff == -1):
                                pass
                            else:
                                neighbor_list.append(k)
                        else:
                            pass
                                    
                    while (int(bondnumbers[n+gap,1])) > 0:
                                                      
                        if len(neighbor_list) != 0:
                            chosen = random.choice(neighbor_list)
                            if ([int(beadname[chosen][0]),int(beadname[n][0])]==b_list[:,None]).all(-1).any() == False and ([int(beadname[n][0]),int(beadname[chosen][0])]==b_list[:,None]).all(-1).any() == False:
                                bondnumbers[int(chosen)+gap,1] -= 1
                                bondnumbers[int(n)+gap,1] -= 1
                                bnum.append(str(ind)+"\t"+"3"+"\t"+str(beadname[n][0])+"\t"+str(beadname[chosen][0])+"\n")
                                b_list[count+b_index,0] = str(beadname[n][0]);b_list[count+b_index,1] = str(beadname[chosen][0])
                                neighbor_list.remove(chosen)
                                ind += 1
                                b_index +=1
                            else:
                                neighbor_list.remove(chosen)
                                            
                        else:
                            break
                                  
                else:
                    pass
                
                        
        mean_bondnum = np.subtract(bondnumb,bondnumbers[:,1])   
        mean_bondnum = np.mean(mean_bondnum)
        print(mean_bondnum,"Crosslink Mean Bond Number")

    
    if crosslinker:        
        Crosslinker(chnum, 0.1, 1, 3)
        print("Crosslinked")
    else:
        pass
    
    def Mutation(mut_ratio):
        global mutated
        global tb_mut
        
        mutated = []
        tb_mut = []
        
        for i in dnaloc:
            if i[1] == '2' or i[1] == '8':
                tb_mut.append(int(i[0])-beads-1)
        
        het_rat = len(tb_mut)/len(dnaloc)
        
        while het_rat > mut_ratio:
            
            chosen = random.choice(tb_mut)
            
            if dnaloc[chosen][1] == '2':            
                dnaloc[chosen][1] = '3';dnaloc[chosen][2] = '3'
                mutated.append(chosen)

            elif dnaloc[chosen][1] == '8':
                dnaloc[chosen][1] = '7';dnaloc[chosen][2] = '7'
                mutated.append(chosen)
                
                
            het_rat = len(tb_mut)/len(dnaloc)
            tb_mut.remove(chosen)
            
     
    if mutation:        
        Mutation(mut_ratio)
        print("Mutated")
    else:
        pass                  
    

    def Angle_Generator(dnaloc,chnum,atomnum = 0, ind = 1, count = 1):
        global annum
        
        annum = []   
        ##Angle for Hetero
        for s in dnaloc:
            if count != chnum*6002-1:
                if (s[1] == "2" or s[1] == "8")  and (dnaloc[count+1][1] == "2" or dnaloc[count+1][1] == "8"):
                    atomnum = int(s[0])
                    annum.append(str(ind)+"\t"+"1"+"\t"+str(atomnum)+"\t"+str(atomnum+1)+"\t"+str(atomnum+2)+"\n")
                    ind +=1
                    count +=1
                else:
                    count += 1
                    pass
                
    Angle_Generator(dnaloc, chnum)
    
    print("Angles are generated...")

    
    def File_Generator(cap,dnaloc):

        
        print("File is being generating")   

        fin= np.concatenate((cap,dnaloc))
        
        fl = open("Chromatin1.txt",'w+')
        fl.write("\n\n")
        fl.write(str(index1-1) +  " atoms\n")
        fl.write(str(len(bnum)) + " bonds\n")
        
        if chnum != 0:
            fl.write(str(len(annum)) +  " angles\n\n")
            
        fl.write("8 atom types\n\n")
        fl.write("3 bond types\n\n")
        fl.write("1 angle types\n\n")
        xb = str(int(-60)) + " " + str(int(60)) + " "+ "xlo" + " "+ "xhi\n"
        yb = str(int(-60))+ " " + str(int(60))+ " "+ "ylo" + " "+ "yhi\n"
        zb = str(-60) + " " + str(60) + " "+ "zlo" + " "+ "zhi\n\n"
        fl.write(xb + yb + zb)
        fl.write("Masses\n\n")
        fl.write("1 1\n")
        fl.write("2 1\n")
        fl.write("3 1\n")
        fl.write("4 1\n")
        fl.write("5 1\n")
        fl.write("6 1\n")
        fl.write("7 1\n")
        fl.write("8 1\n\n")
        
        fl.write("Atoms\n\n")
        
        
        for line in fin:
            fl.write(" ".join(line) + "\n")
            
        
        
        fl.write("\n")
        fl.write("Bonds\n\n")
        
        
        
        for count,k in enumerate(bnum):
            fl.write(k)    
            
                
        fl.write("\n")
        if chnum != 0:
            fl.write("Angles\n\n") 
        
        for count,k in enumerate(annum):
            fl.write(k)   
          
        fl.close()
        
        file = open("Chromatin1.txt","r")
        file1 = open("data.Chromatin","w")      
          
        x = file.readlines()
        file1.writelines(x)
        file1.close()
        timediff = timeit.default_timer()-start
        print(str(int(timediff/60))+" "+str("minutes") + " " + str(float('{:.2f}'.format(timediff))%60)+" " + str("seconds"))   
        print(str(int(timediff/3600))+" "+str("hours"))
        
    File_Generator(cap, dnaloc)
          
        
Nucleus_Generate(beads=15000,chnum=8,place_factor=1.3,minbond=5,maxbond=7,factor=1.5,mut_ratio=0.3,crosslinker=False,mutation=False)
