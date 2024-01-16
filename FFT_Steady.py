# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 01:05:29 2022

@author: Asus
"""
import numpy as np
import math

files = ["dump.final","dump.final2"]

filenum = 0
mode = 230

while filenum != len(files):
    
    i = files[filenum]    

    filenum += 1
    dump = open(i,"r")
    lines = dump.readlines()
    dump.close()
    
    atom_num = int(lines[3])
    
    
    time_step = int((len(lines))/atom_num)
    half_step = int(time_step/2)
    count_step = half_step
    step = (atom_num + 9)*half_step + 9

    
    dev_array = np.array([])
    correlation = np.array([])
    amp_array = np.zeros((mode,1))

    while count_step != time_step:
 
        l_list = []
        
        for i in lines[step:atom_num+step+1]:
            i = i.split(' ')
            if i[1] == "1":
                l_list.append(i[0]+"\t"+i[2]+"\t"+i[3]+"\t"+i[4])    
        
        rad_dists = np.array([])

        ##Defining a slice
        for s,i in enumerate(l_list):
            i = i.split("\t")
            if float(i[1]) <= 0.5 and float(i[1]) >=-0.5:
                rad_dists = np.append(rad_dists,i)

            else:
                pass

        rad_dists = np.reshape(rad_dists,(int(len(rad_dists)/4),4))
        rad_dists = rad_dists.astype("float")
        
        
        angle_array = np.array([])

        ##Finding angles
        for k in rad_dists:
            
            angle = math.atan2(k[2] - 0, k[3] - 0)
            angle = math.degrees(angle)
            angle_array = np.append(angle_array,angle)
        
        angle_array = np.reshape(angle_array,(int(len(rad_dists)),1))
                
        ##Radial distances of beads on a slice
        rad_dists = rad_dists[:,1:]**2
        rad_dists = np.sum(rad_dists,axis=1)
        rad_dists = np.sqrt(rad_dists)
        
        ##Average Radial Distance
        rad_average = np.average(rad_dists)    
        
        ##R_i - R_avg
        rad_diff = np.subtract(rad_dists,rad_average)
        rad_diff = np.reshape(rad_diff,(len(rad_diff),1))
        rad_diff = np.hstack((rad_diff,angle_array))
        
        ##Sorting by angle
        
        random_bead = rad_diff[np.random.choice(rad_diff.shape[0], mode, replace=False)]

        rad_diff = random_bead[random_bead[:,1].argsort()[::-1]]
    
        ##Correlation
        # mean_rad = np.mean(rad_diff)
        # correlation = np.append(correlation,mean_rad)
        
        ##STD by time
        deviation = np.std(rad_diff[:,0])
        dev_array = np.append(dev_array,deviation)
        
        ##FFT
        fft = np.fft.fft(rad_diff[:,0])
        amplitude = np.abs(fft)
        amplitude = np.reshape(amplitude,(len(amplitude),1))
        amp_array += amplitude
        
        count_step += 1
        step += atom_num + 9      
        

    dev_average = np.average(dev_array)

    amp_array = amp_array / time_step
    
    amp_array *= amp_array
    
    amp_array = amp_array[:int(len(amp_array)/2)]

    # Mean
    # mean = np.mean(correlation)
    
    # # Variance
    # var = np.var(correlation)
    
    # # # Normalized data
    # ndata = correlation - mean
    
    # acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:] 
    # acorr = acorr / var / len(ndata)
    
    wavenum = np.linspace(1,int(mode/2),int(mode/2))
    wavenum = np.reshape(wavenum,(int(mode/2),1))
    amp_array = np.hstack((wavenum,amp_array))
    amp_array = amp_array[2:]
    
    f = open("dev_average_half"+str(filenum)+".txt", "a")
    f.write(str(dev_average))
    f.close() 

    file = np.savetxt("amp_array_half"+str(filenum)+".csv",amp_array,delimiter=',')    
    # file1 = np.savetxt("AutoCorr"+str(filenum)+".csv",acorr, delimiter = ',')
