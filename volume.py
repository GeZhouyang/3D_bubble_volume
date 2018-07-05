import sys
import time
import numpy as np
from scipy import ndimage


#--- Input a list of pairs, output a redundent list of connections.
#    (Some components may be a subset of others.)

def mm1(merge_list):

    merge_list_1 = []    
    lm = len(merge_list)
    
    if lm > 1:
        for p in range(lm): 
            ap = merge_list[p]
            buf = ap.copy()
            for q in range(lm):
                if q != p:
                    aq = merge_list[q] # aq can only have 2 elements
                    if aq[0] in buf and aq[1] not in buf: buf.append(aq[1])
                    if aq[1] in buf and aq[0] not in buf: buf.append(aq[0])
                    
            merge_list_1.append(buf)

    return merge_list_1


#--- Input a list of possibly incomplete, connected groups,
#    output the list of complete, unique connections.
#    (Double checking missing pairs due to forward sweeping in mm1.
#    Two sweeps guarentee the grouping of all pairs.)

def mm2(merge_list_1):

    merge_list_2,merge_list_3 = [],[]
    lm = len(merge_list_1)
    
    if lm > 1:            
        for p in range(lm):
            ap = merge_list_1[p]
            lp = len(ap)
            buf = ap.copy()
            for q in range(lm):
                if q != p:
                    aq = merge_list_1[q] # aq can have more than 2 elements
                    lq = len(aq)
                    for r in range(lq):
                        if aq[r] in ap:
                            buf = ap if lp >= lq else aq

            merge_list_2.append(buf)
            
        for i in range(len(merge_list_2)):
            a = sorted(merge_list_2[i])
            if a not in merge_list_3:
                merge_list_3.append(a)

    return merge_list_3


#--- Input a list of possibly multi-connected groups,
#    output the list of uniquely merged components.
#    (Some components may be a subset of others.)

def mm3(merge_list):

    merge_list_3,merge_list_4 = [],[]
    lm = len(merge_list)

    if lm > 1:
        for p in range(lm):
            ap = merge_list[p]
            buf = ap.copy()
            for q in range(lm):
                if q != p:
                    aq = merge_list[q]
                    lq = len(aq)
                    flag = False
                    for r in range(lq):
                        if aq[r] in ap: flag = True
                    if flag:
                        for r in range(lq):
                            if aq[r] not in ap: buf.append(aq[r])
                            
            merge_list_3.append(buf)

        for i in range(len(merge_list_3)):
            a = sorted(merge_list_3[i])
            if a not in merge_list_4:
                merge_list_4.append(a)
    
    return merge_list_4


#--- Input cross-sectional areas on the shear planes and their sorting,
#    output a redundent merge list using the sorted area.

def sm(area0,area1,s0,s1):

    shear = []
    i,i0,i1 = 0,0,0
    while area0[s0[i+i0]] != 0. and area1[s1[i+i1]] != 0.:
        if s0[i+i0] != 0 and s1[i+i1] != 0:
            a = [s0[i+i0],s1[i+i1]]
            
            j = i+i0
            while area0[s0[j+1]] == area0[s0[j]] and area0[s0[j+1]] != 0.:
                a += [s0[j+1]]
                i0 += 1
                j += 1
                
            j = i+i1
            while area1[s1[j+1]] == area1[s1[j]] and area1[s1[j+1]] != 0.:
                a += [s1[j+1]]
                i1 += 1
                j += 1
                
        shear.append(a)
        i += 1

    for p in range(len(shear)):
        shear[p] = sorted(shear[p])

    return shear


#---- Main script

if __name__  == '__main__':

    t_tot = time.time()
    
    #--> Load the data

    tstep = 82000

    print('Loading data file '+str(tstep)+'...\n')

    fn = 'vof'+str(tstep).zfill(7)+'_3d.bin'

    vof_raw = np.fromfile(fn)

    Lx,Ly,Lz = 20., 19., 40.  # domain size
    Nx,Ny,Nz = 160, 152, 320  # number of points
    dV = Lx*Ly*Lz/(Nx*Ny*Nz)
    
    vof = np.reshape(vof_raw,(Nx,Ny,Nz), order='F')

    #--> Count and label the bubbles in the entire domain

    print('First count in 3D...')
        
    nb    = np.ones(shape=(3,3,3))  # neighbor matrix
    nb_2D = np.ones(shape=(3,3))
    nb_1D = np.ones(3)
    thres = 0.7  # threshold value (Total volume converges as thres -> 0.)
    
    labels, num = ndimage.label(vof > thres, nb)
    
    print(num,'bubbles found.\n')

    #--> Merge bubbles across periodic planes

    merge_list = []

    # regular periodic planes
    
    print('Finding split bubbles in regular periodic planes...')
    
    vof_xz0, vof_xz1 = vof[:,0,:], vof[:,-1,:]    
    label_xz0, num_xz0 = ndimage.label(vof_xz0 > thres, nb_2D)
    label_xz1, num_xz1 = ndimage.label(vof_xz1 > thres, nb_2D)
    
    for i in range(Nx):
        for k in range(Nz):
            if label_xz0[i,k] and label_xz1[i,k]:
                a = [labels[i,0,k],labels[i,-1,k]]
                if a not in merge_list:
                    merge_list.append(a)
                    
    vof_xy0, vof_xy1 = vof[:,:,0], vof[:,:,-1]
    label_xy0, num_xy0 = ndimage.label(vof_xy0 > thres, nb_2D)
    label_xy1, num_xy1 = ndimage.label(vof_xy1 > thres, nb_2D)
    
    for i in range(Nx):
        for j in range(Ny):
            if label_xy0[i,j] and label_xy1[i,j]:
                a = [labels[i,j,0],labels[i,j,-1]]
                if a not in merge_list:
                    merge_list.append(a)
    
    # axis perpendicular to the shear plane (diagonal) [redundent]
    
    vof_x0, vof_x1 = vof[:,0,0], vof[:,-1,-1]
    label_x0, num_x0 = ndimage.label(vof_x0 > thres, nb_1D)
    label_x1, num_x1 = ndimage.label(vof_x1 > thres, nb_1D)

    for i in range(Nx):
        if label_x0[i] and label_x1[i]:
            a = [labels[i,0,0],labels[i,-1,-1]]
            if a not in merge_list:
                merge_list.append(a)
                
    print('Pairs to be merged (may contain multi-connections):\n',merge_list)

    # merge lists in the merge list

    merge_list_1 = mm1(merge_list)
    merge_list = mm2(merge_list_1)

    print('All gathered:\n',merge_list,'\n')
                
    # shear periodic planes
    
    print('Finding split bubbles in shear periodic planes...')

    vof_yz0, vof_yz1 = vof[0,:,:], vof[-1,:,:]   
    label_yz0, num_yz0 = ndimage.label(vof_yz0 > thres, nb_2D)
    label_yz1, num_yz1 = ndimage.label(vof_yz1 > thres, nb_2D)

    area0 = np.zeros(num+1) # use cross-sectional area to determine merging
    area1 = np.zeros(num+1)
    list0, list1 = [],[]
    
    for j in range(Ny):
        for k in range(Nz):
            
            if label_yz0[j,k]:
                l = labels[0,j,k] # labels in 3D
                area0[l] += vof[0,j,k]*dV
                if l not in list0: list0.append(l)
           
            if label_yz1[j,k]:
                l = labels[-1,j,k] # labels in 3D
                area1[l] += vof[-1,j,k]*dV
                if l not in list1: list1.append(l)

    list0 = sorted(list0)
    list1 = sorted(list1)
    
    print('Bubbles crossing the bottom shear plane:',list0)
    print('Bubbles crossing the top shear plane:'   ,list1)

    buf = np.zeros(num+1)
    for i in range(num+1): buf[i] = area0[i]
    for i in range(len(list0)):
        l = list0[i]
        for j in range(len(merge_list)):
            a = merge_list[j]
            if l in a: buf[l] = sum([area0[k] for k in a])
    area0 = buf
                
    buf = np.zeros(num+1)
    for i in range(num+1): buf[i] = area1[i]
    for i in range(len(list1)):
        l = list1[i]
        for j in range(len(merge_list)):
            a = merge_list[j]
            if l in a: buf[l] = sum([area1[k] for k in a])
    area1 = buf
    
    s0,s1 = [],[]
    s0,s1 = np.argsort(-area0),np.argsort(-area1) # minus to sort from max
        
    shear = sm(area0,area1,s0,s1)
    
    print('Groups to be merged:\n',shear,'\n')

    merge_list += shear
    
    # merging hopefully completed

    merge_list = mm3(merge_list)
    merge_list = mm2(merge_list)

    num_m = len(merge_list)

    print('Final merging list:\n',merge_list,'\n')    

    # Calculate the individual volumes (including splitted bubbles)

    print('Calculate the split volumes (the slowest part) ...')
    
    #t_spl = time.time()
    V_spl = np.zeros(num+1)
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                l = labels[i,j,k]
                if l: V_spl[l] += vof[i,j,k]*dV

    V1 = sum(V_spl)
    #t_spl = time.time()- t_spl
    #print('Splitting time = ',t_spl)

    print('Done.')

    # Volume correction
    
    V_corr = V_spl.copy()  
    for i in range(num_m):
        a = merge_list[i]
        for j in range(len(a)):
            V_corr[a[j]] = sum([V_spl[k] for k in a])
    V_corr = V_corr[1:] # remove label zero (background)

    # Final output (unique values)

    V_output = []
    for i in range(num):
        if V_corr[i] not in V_output: V_output.append(V_corr[i])

    V_output = sorted(V_output)
    
    # Check errors
    
    V0 = sum(vof_raw)*dV
    V_tot = sum(V_output)

    ## should be less than a few % (depending on thres)
    print('\nRelative loss of unlabeled volume = ',(V_tot-V0)/V0)
    ## should be machine zero
    print('Volume merging error = ',(V_tot-V1)/V1)

    # Save to file
    
    index = list(range(1,len(V_output)+1))
    data_array = np.column_stack((index, V_output)) # convert multiple lists to one matrix
    
    np.savetxt('ivol'+str(tstep).zfill(7)+'.txt',data_array, fmt=('%8i','%16.8f'))

    t_tot = time.time() - t_tot
    print('\nTotal wall time = ',t_tot, 'seconds.')
