# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 23:06:35 2021

@author: Asus
"""
import numpy as np
import math
import csv


def DumpRead():
    files = ["dump.final","dump.final2"]
    
    filenum = 0
    while filenum != len(files):
        i = files[filenum]    
        rad_dist_h = []
        rad_dist_e = []
        rad_dist_c = []
        rad_dist_l = []
        rad_dist_l_x = []
        rad_dist_l_y = []
        rad_dist_l_z = []
        volume_frac = []
    
        filenum += 1
        dump = open(i,"r")
        lines = dump.readlines()
        dump.close()
        
        atom_num = int(lines[3])
        
        time_step = int((len(lines))/atom_num)
        count_step = 0
        step = 9
    
        while count_step != time_step:
            h_list = []
            e_list = []
            c_list = []
            l_list = []
            
            for i in lines[step:atom_num+step+1]:
                i = i.split(' ')
                if i[1] == "2" or i[1] == "8":
                    h_list.append(i[2]+"\t"+i[3]+"\t"+i[4])
                elif i[1] == "3":
                    e_list.append(i[2]+"\t"+i[3]+"\t"+i[4])
                elif i[1] == "4" or i[1] == "6":
                    c_list.append(i[2]+"\t"+i[3]+"\t"+i[4])
                elif i[1] == "1":
                    l_list.append(i[2]+"\t"+i[3]+"\t"+i[4])    
            
            hdist = 0
            for i in h_list:
                i = i.split("\t")
                dist1 = math.sqrt((float(i[0]))**2+(float(i[1]))**2+((float(i[2])))**2)
                hdist += dist1
            
            hdist = hdist/len(h_list)        
            rad_dist_h.append(hdist)
            
            
            dist = 0
            for i in e_list:
                i = i.split("\t")
                dist1 = math.sqrt((float(i[0]))**2+(float(i[1]))**2+((float(i[2])))**2)
                dist += dist1
            
            dist = dist/len(e_list)        
            rad_dist_e.append(dist)
        
            
            dist = 0
            for i in c_list:
                i = i.split("\t")
                dist1 = math.sqrt((float(i[0]))**2+(float(i[1]))**2+((float(i[2])))**2)
                dist += dist1
                
            dist = dist/len(c_list)
                
            rad_dist_c.append(dist)
            
            dist = 0
            dist_x = 0
            dist_y = 0
            dist_z = 0
            for i in l_list:
                i = i.split("\t")
                dist_x += abs(float(i[0]))
                dist_y += abs(float(i[1]))
                dist_z += abs(float(i[2]))
                dist1 = math.sqrt((float(i[0]))**2+(float(i[1]))**2+((float(i[2])))**2)
                dist += dist1
            
            
            dist = dist/len(l_list)
            
            Vol_fraction = float(6002/(dist**3))
            
            dist_x = dist_x/len(l_list)
            dist_y = dist_y/len(l_list)
            dist_z = dist_z/len(l_list)
            
            volume_frac.append(Vol_fraction)
            rad_dist_l.append(dist)
            rad_dist_l_x.append(dist_x)
            rad_dist_l_y.append(dist_y)
            rad_dist_l_z.append(dist_z)
            
            
            count_step += 1
            step += atom_num + 9
        
        raddists = np.array([rad_dist_h,rad_dist_e,rad_dist_c,rad_dist_l,volume_frac])
        raddists = np.savetxt("Raddists"+str(filenum)+".csv",raddists,delimiter=',')
    
    
        densityh = np.array([])
        densitye = np.array([])
        densityc = np.array([])
        step = 0.5
        maxrad = 30
        maxvol = (4/3)*math.pi*maxrad**3
        totaldna = len(h_list)+len(e_list)+len(c_list)
    
        for k in np.arange(0,maxrad+step,step):
            atomc = 0
            for i in h_list:
                i = i.split("\t")
                dist1 = math.sqrt((float(i[0]))**2+(float(i[1]))**2+((float(i[2])))**2)
                if dist1<k and dist1>k-step:
                    atomc += 1       
            relative = float(atomc/float(len(h_list)))
            densityh = np.append(densityh,(relative*100,100*((4/3)*math.pi*k**3)/maxvol))
            
            atomc = 0
            for i in e_list:
                i = i.split("\t")
                dist1 = math.sqrt((float(i[0]))**2+(float(i[1]))**2+((float(i[2])))**2)
                if dist1<k and dist1>k-step:
                    atomc += 1        
            relative = float(atomc/float(len(e_list)))
            densitye = np.append(densitye,(relative*100,100*((4/3)*math.pi*k**3)/maxvol))
            
            atomc = 0
            for i in c_list:
                i = i.split("\t")
                dist1 = math.sqrt((float(i[0]))**2+(float(i[1]))**2+((float(i[2])))**2)
                if dist1<k and dist1>k-step:
                    atomc += 1        
            relative = float(atomc/float(len(c_list)))
            densityc = np.append(densityc,(relative*100,100*((4/3)*math.pi*k**3)/maxvol)) 
    
        densityh = densityh.reshape(int((maxrad)/step)+1,2)
        densitye = densitye.reshape(int((maxrad)/step)+1,2)
        densityc = densityc.reshape(int((maxrad)/step)+1,2)
    
    
        densities = np.array([densityh[:,1],densityh[:,0],densitye[:,0],densityc[:,0]])
        densities = np.savetxt("densities"+str(filenum)+".csv",densities,delimiter=',')
    
        xyz_change = np.array([rad_dist_l_x,rad_dist_l_y,rad_dist_l_z])
        xyz_change = np.savetxt("xyz_change"+str(filenum)+".csv",xyz_change,delimiter=',')


def tensor():
    import numpy as np
    from numpy import linalg as LA
    
    files = ["Con_Rg.txt","Inv_Rg.txt"]
    
    filenum = 0
    while filenum != len(files):
        i = files[filenum]   
        dump = open(i,"r")
        lines = dump.readlines()
        lines = lines[3:]
        dump.close()
        
        step = 0
        tensor_list = []
        while step < len(lines):
            for k in lines[step+1:step+2]:
                tensor_list.append(k)
            step += 2
        
        whole_data = []
        for t in tensor_list:
            tensor_array = np.zeros((3,3))
            t = t.split(" ")
            tensor_array[0,0] = t[1]
            tensor_array[0,1] = t[4]
            tensor_array[0,2] = t[5]
            tensor_array[1,0] = t[4]
            tensor_array[1,1] = t[2]
            tensor_array[1,2] = t[6]
            tensor_array[2,0] = t[5]
            tensor_array[2,1] = t[6]
            tensor_array[2,2] = t[3]
            
            eigvals = LA.eigvals(tensor_array)
            whole_data.append(eigvals)
        
        tensor = np.savetxt("Tensor_"+str(filenum)+".csv",whole_data,delimiter=',')    
        filenum += 1

DumpRead()
tensor()

