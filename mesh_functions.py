from igakit.cad import *
from igakit.cad import extrude
from igakit.igalib import *
from igakit.plot import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
import math
from igakit.io import PetIGA,VTK
from numpy import linspace
import sys
np.set_printoptions(threshold=sys.maxsize)
import glob
import stretching_functions as sf

#Utilities
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


# Functions
def generate_unif_background(C1, C2, t, num_div_vec, n=1):
    #Input top and bottom line of rectangular domain as 'igakit cad line', thickness, t, a vector [numknotsx, numknotsy, numknotsh] as np.array, returns NURBS object mesh of domain, S. Generates a uniform rectangular background grid. Implementation of p and h refinement procedure.
    S=ruled(C1, C2)
    S.elevate(1,n)
    S.elevate(0,n)
    if t>0:
        S=extrude(S, displ=t, axis=2)
        S.elevate(2,n)

    S.refine(0, np.linspace(0,1,num_div_vec[0])[1:-1])
    S.refine(1, np.linspace(0,1,num_div_vec[1])[1:-1])
    if num_div_vec[2]>1: S.refine(2, np.linspace(0,1,num_div_vec[2])[1:-1])

    return(S)

def generate_Shell(C1, C2, t, num_div_vec, n=1):
    #Input top and bottom line of rectangular domain as 'igakit cad line', thickness, t, a vector [numknotsx, numknotsy, numknotsh] as np.array, returns NURBS object mesh of domain, S. Generates a uniform rectangular background grid. Implementation of p and h refinement procedure.
    S=ruled(C1, C2)
    S.elevate(1,n)
    S.elevate(0,n)
    #if t>0:
    S=extrude(S, displ=t, axis=2)
    S.elevate(2,n)

    S.refine(0, np.linspace(0,1,num_div_vec[0])[1:-1])
    S.refine(1, np.linspace(0,1,num_div_vec[1])[1:-1])
    if num_div_vec[2]>1: S.refine(2, np.linspace(0,1,num_div_vec[2])[1:-1])

    return(S)

def generate_unif_foreground(origin, xyzmax_vec, xyz_numpts_vec, material):
    # Generate volumes and foreground points from input origin, 3 orthogonal max displacements in x, y and z as [xmax, ymax, zmax], and a vector of the number of points in x, y and z to generate as np.array([numptsx, numptsy, numptsz])
    node_num=xyz_numpts_vec[0]*xyz_numpts_vec[1]*xyz_numpts_vec[2]
    x=np.zeros(xyz_numpts_vec[0])
    y=np.zeros(xyz_numpts_vec[1])
    h=np.zeros(xyz_numpts_vec[2])

    dx_vec=np.zeros(node_num)
    dy_vec=np.zeros(node_num)
    dh_vec=np.zeros(node_num)

    xr_vec=np.zeros(node_num)
    xl_vec=np.zeros(node_num)
    yr_vec=np.zeros(node_num)
    yl_vec=np.zeros(node_num)
    hr_vec=np.zeros(node_num)
    hl_vec=np.zeros(node_num)

    mat_vec=np.zeros(node_num)

    coor=np.zeros((int(node_num),4))
    print 'Generating Uniform foreground...'
    print "%s Nodes" %node_num
    vol=np.zeros(node_num)
    dx=(xyzmax_vec[0]-origin[0])/(xyz_numpts_vec[0]-1)
    dy=(xyzmax_vec[1]-origin[1])/(xyz_numpts_vec[1]-1)
    if xyz_numpts_vec[2]>1: dh=(np.double(xyzmax_vec[2])-origin[2])/(np.double(xyz_numpts_vec[2])-1)
    if xyz_numpts_vec[2]==1: dh=0
    node=0
    for i in range(xyz_numpts_vec[2]):
        if i>0: h[i]=h[i-1]+dh
        if i==0: h[i]=origin[2]
        hr=h[i]+dh/2.0
        hl=h[i]-dh/2.0
        if hr>xyzmax_vec[2]: hr=xyzmax_vec[2] #Comment these lines to switch to uniform volume
        if hl<0: hl=origin[2] #Comment these lines to switch to uniform volume
        for j in range(xyz_numpts_vec[0]):
            if j>0: x[j]=x[j-1]+dx
            if j==0: x[j]=origin[0]
            for k in range(xyz_numpts_vec[1]):
                xr=x[j]+dx/2.0
                xl=x[j]-dx/2.0
                if xr>xyzmax_vec[0]: xr=xyzmax_vec[0] #Comment these lines to switch to uniform volume
                if xl<origin[0]: xl=origin[0] #Comment these lines to switch to uniform volume
                if k>0: y[k]=y[k-1]+dy
                if k==0: y[k]=origin[1]
                yr=y[k]+dy/2.0
                yl=y[k]-dy/2.0
                if yr>xyzmax_vec[1]: yr=xyzmax_vec[1] #Comment these lines to switch to uniform volume
                if yl<origin[1]: yl=origin[1] #Comment these lines to switch to uniform volume
                vol[node]=(xr-xl)*(yr-yl)
                if hr!=hl: vol[node]=vol[node]*(hr-hl)
                coor[node,:]=[node+1,x[j],y[k],h[i]]
                xr_vec[node]=xr
                xl_vec[node]=xl
                yr_vec[node]=yr
                yl_vec[node]=yl
                hr_vec[node]=hr
                hl_vec[node]=hl
                dx_vec[node]=xr-xl
                dy_vec[node]=yr-yl
                dh_vec[node]=hr-hl
                mat_vec[node]=material
                node=node+1
                if node%1000000==0:
                    perc=np.floor(np.float(node)/np.float(node_num)*100.0)
                    str="%s %%" %perc
                    print str
    coor=np.vstack(([node,0,0,0], coor))

    ## Test statements
    #if xyzmax_vec[2]>0:
    #   assert np.absolute(np.sum(vol)-(xyzmax_vec[0]-origin[0])*(xyzmax_vec[1]-origin[1])*(xyzmax_vec[2]-origin[2]))<=10**(-12)*xyz_numpts_vec[0]*xyz_numpts_vec[1]*xyz_numpts_vec[2]
    #elif xyzmax_vec[2]==0:
    #    assert np.absolute(np.sum(vol)-(xyzmax_vec[0]-origin[0])*(xyzmax_vec[1]-origin[1]))<=10**(-12)*xyz_numpts_vec[0]*xyz_numpts_vec[1]*xyz_numpts_vec[2]
    ##

    G=sf.foreground(coor, vol, xyz_numpts_vec, dx_vec, dy_vec, dh_vec, xr_vec, xl_vec, yr_vec, yl_vec, hr_vec, hl_vec, mat_vec)
    print "Complete"
    return(G)

def generate_FGCone_cylindrical(origin, r1, r2, r_min, height, xyz_numpts_vec, material, filled, axis = 2):
    # xyz numpts vec is [layers, particles per ring, num_divZ]
    h=np.zeros(xyz_numpts_vec[2])
    dh=np.double(height-origin[2])/(np.double(xyz_numpts_vec[2])-1)

    for i in range(xyz_numpts_vec[2]):
        if i>0: h[i]=h[i-1]+dh
        if i==0: h[i]=origin[2]
        hr=h[i]+dh/2.0
        hl=h[i]-dh/2.0
        if hr>height: hr=height #Comment these lines to switch to uniform volume
        if hl<0: hl=origin[2] #Comment these lines to switch to uniform volume
        r_crit = r1-h[i]*(r1-r2)/height;
        if(i ==0):
            G = GetConeLayer(origin, h[i], hr, hl, r_crit, r_min, xyz_numpts_vec, material, filled)
        elif(i>0):
            G1 = GetConeLayer(origin, h[i], hr, hl, r_crit, r_min, xyz_numpts_vec, material, filled)
            G = sf.fg_superpose(G, G1)
    print "Complete"
    return(G)

def GetConeLayer(origin, h, hr, hl, r_crit, r_min, xyz_numpts_vec, material, filled):
    arc = 2.0*np.pi/xyz_numpts_vec[1]
    layers = xyz_numpts_vec[0]
    node_num=xyz_numpts_vec[0]*xyz_numpts_vec[1]

    r = np.zeros(xyz_numpts_vec[0])
    theta = np.zeros(xyz_numpts_vec[1])

    dr_vec=np.zeros(node_num)
    dtheta_vec=np.zeros(node_num)
    dh_vec=np.zeros(node_num)

    rr_vec=np.zeros(node_num)
    rl_vec=np.zeros(node_num)
    thetar_vec=np.zeros(node_num)
    thetal_vec=np.zeros(node_num)
    hr_vec=np.zeros(node_num)
    hl_vec=np.zeros(node_num)

    mat_vec=np.zeros(node_num)

    coor = np.zeros((int(node_num),4))
    vol = np.zeros(node_num)
    dr = (r_crit-r_min)/(layers-1)
    dtheta = arc

    node=0
    for j in range(layers):
        if j>0: r[j]=r[j-1]+dr
        if j==0: r[j]=r_min
        for i in range(xyz_numpts_vec[1]):
            if i==0: theta[i]= 0.0
            if i>0:  theta[i] = theta[i-1] + dtheta
            rr=r[j]+dr/2.0
            rl=r[j]-dr/2.0
            if rr>r_crit: rr=r_crit #Comment these lines to switch to uniform volume
            if rl<origin[0]: rl=r_min*(filled%1) #Comment these lines to switch to uniform volume
            thetar=theta[i] + dtheta/2.0
            thetal=theta[i] - dtheta/2.0
            [x, y] = pol2cart(r[j], theta[i])
            coor[node,:]=[node+1, x, y ,h]
            vol[node] = (hr-hl)*(rr**2-rl**2)*arc/2.0
            rr_vec[node]=rr
            rl_vec[node]=rl
            thetar_vec[node]=thetar
            thetal_vec[node]=thetal
            hr_vec[node]=hr
            hl_vec[node]=hl
            dr_vec[node]=rr-rl
            dtheta_vec[node]=thetar-thetal
            dh_vec[node]=hr-hl
            mat_vec[node]=material
            node=node+1
    coor=np.vstack(([node,0,0,0], coor))
    G=sf.foreground(coor, vol, xyz_numpts_vec, dr_vec, dtheta_vec, dh_vec, rr_vec, rl_vec, thetar_vec, thetal_vec, hr_vec, hl_vec, mat_vec)
    return(G)


def generate_unif_PDforeground(origin, xyzmax_vec, xyz_numpts_vec, material):
    # Generate volumes and foreground points from input origin, 3 orthogonal max displacements in x, y and z as [xmax, ymax, zmax], and a vector of the number of points in x, y and z to generate as np.array([numptsx, numptsy, numptsz])
    node_num=xyz_numpts_vec[0]*xyz_numpts_vec[1]*xyz_numpts_vec[2]
    x=np.zeros(xyz_numpts_vec[0])
    y=np.zeros(xyz_numpts_vec[1])
    h=np.zeros(xyz_numpts_vec[2])

    dx_vec=np.zeros(node_num)
    dy_vec=np.zeros(node_num)
    dh_vec=np.zeros(node_num)

    xr_vec=np.zeros(node_num)
    xl_vec=np.zeros(node_num)
    yr_vec=np.zeros(node_num)
    yl_vec=np.zeros(node_num)
    hr_vec=np.zeros(node_num)
    hl_vec=np.zeros(node_num)

    mat_vec=np.zeros(node_num)

    coor=np.zeros((int(node_num),4))
    print 'Generating Uniform foreground...'
    print "%s Nodes" %node_num
    vol=np.zeros(node_num)
    dx=(xyzmax_vec[0]-origin[0])/(xyz_numpts_vec[0]-1)
    dy=(xyzmax_vec[1]-origin[1])/(xyz_numpts_vec[1]-1)
    if xyz_numpts_vec[2]>1: dh=(np.double(xyzmax_vec[2])-origin[2])/(np.double(xyz_numpts_vec[2])-1)
    if xyz_numpts_vec[2]==1: dh=0
    node=0
    for i in range(xyz_numpts_vec[2]):
        if i>0: h[i]=h[i-1]+dh
        if i==0: h[i]=origin[2]
        hr=h[i]+dh/2.0
        hl=h[i]-dh/2.0
        if hr>xyzmax_vec[2]: hr=xyzmax_vec[2] #Comment these lines to switch to uniform volume
        if hl<0: hl=0 #Comment these lines to switch to uniform volume
        for k in range(xyz_numpts_vec[1]):
            if k>0: y[k]=y[k-1]+dy
            if k==0: y[k]=origin[1]
            yr=y[k]+dy/2.0
            yl=y[k]-dy/2.0
            if yr>xyzmax_vec[1]: yr=xyzmax_vec[1] #Comment these lines to switch to uniform volume
            if yl<origin[1]: yl=origin[1] #Comment these lines to switch to uniform volume
            for j in range(xyz_numpts_vec[0]):
                if j>0: x[j]=x[j-1]+dx
                if j==0: x[j]=origin[0]
                xr=x[j]+dx/2.0
                xl=x[j]-dx/2.0
                if xr>xyzmax_vec[0]: xr=xyzmax_vec[0] #Comment these lines to switch to uniform volume
                if xl<origin[0]: xl=origin[0] #Comment these lines to switch to uniform volume
                vol[node]=(xr-xl)*(yr-yl)
                if hr!=hl: vol[node]=vol[node]*(hr-hl)
                if hr==hl: vol[node]=vol[node]*xyzmax_vec[2] #Specifically for shells, not for PD solid in 2D
                coor[node,:]=[node+1,x[j],y[k],h[i]]
                xr_vec[node]=xr
                xl_vec[node]=xl
                yr_vec[node]=yr
                yl_vec[node]=yl
                hr_vec[node]=hr
                hl_vec[node]=hl
                dx_vec[node]=xr-xl
                dy_vec[node]=yr-yl
                dh_vec[node]=hr-hl
                mat_vec[node]=material
                node=node+1
                if node%1000000==0:
                    perc=np.floor(np.float(node)/np.float(node_num)*100.0)
                    str="%s %%" %perc
                    print str
    coor=np.vstack(([node,0,0,0], coor))

    ## Test statements
    #if xyzmax_vec[2]>0:
    #   assert np.absolute(np.sum(vol)-(xyzmax_vec[0]-origin[0])*(xyzmax_vec[1]-origin[1])*(xyzmax_vec[2]-origin[2]))<=10**(-12)*xyz_numpts_vec[0]*xyz_numpts_vec[1]*xyz_numpts_vec[2]
    #elif xyzmax_vec[2]==0:
    #    assert np.absolute(np.sum(vol)-(xyzmax_vec[0]-origin[0])*(xyzmax_vec[1]-origin[1]))<=10**(-12)*xyz_numpts_vec[0]*xyz_numpts_vec[1]*xyz_numpts_vec[2]
    ##

    G=sf.foreground(coor, vol, xyz_numpts_vec, dx_vec, dy_vec, dh_vec, xr_vec, xl_vec, yr_vec, yl_vec, hr_vec, hl_vec, mat_vec)
    print "Complete"
    return(G)
def vis_background():
    # Generate VTK from backround geometry NURBS dat file, openable in Paraview
    nrb = PetIGA().read("Geometry.dat")
    # write a function to sample the nrbs object (100 points from beginning to end)
    uniform = lambda U: linspace(U[0], U[-1], 150)
    outfile = "Geometry" + ".vtk"
    # write a binary VTK file
    VTK().write(outfile,
                nrb,
                )

def vis_foreground2d(G):
    # Generate Scatter plot of saved point cloud data
    plt.scatter(G.coor[1:,1], G.coor[1:,2])
    plt.show()

def save_geometry(G, S, processor_num):
    # outputs input_coor.dat, XVOL.dat, and Geometry.dat from input coordinates, volumes, and NURBS obj
    print "Saving input files..."
    coor=G.coor
    vol=G.vols
    dx_vec=np.transpose(np.array([G.dx_vec]))
    dy_vec=np.transpose(np.array([G.dy_vec]))
    dz_vec=np.transpose(np.array([G.dh_vec]))
    mat_vec=np.transpose(np.array([G.mat]))

    nodes=coor.shape[0]-1

    temp_var1=nodes%processor_num #number of remaining nodes added to last file
    for i in range(processor_num-1):
        num=np.transpose(np.array([[(i+1)*(nodes-temp_var1)/processor_num+1-(i*(nodes-temp_var1)/processor_num+1)],[0],[0],[0],[0],[0],[0],[0],[0]]))
        lines=np.hstack((coor[i*(nodes-temp_var1)/processor_num+1:(i+1)*(nodes-temp_var1)/processor_num+1,:], dx_vec[i*(nodes-temp_var1)/processor_num:(i+1)*(nodes-temp_var1)/processor_num], dy_vec[i*(nodes-temp_var1)/processor_num:(i+1)*(nodes-temp_var1)/processor_num], dz_vec[i*(nodes-temp_var1)/processor_num:(i+1)*(nodes-temp_var1)/processor_num], mat_vec[i*(nodes-temp_var1)/processor_num:(i+1)*(nodes-temp_var1)/processor_num], np.transpose(np.array([vol[i*(nodes-temp_var1)/processor_num:(i+1)*(nodes-temp_var1)/processor_num]]))))
        np.savetxt('foreground%d.dat'%(i+1, ), np.vstack((num,lines)), fmt='%7i %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %7i' '%20.15f')
        #np.savetxt('XVOL%d.dat'%(i+1, ), vol[i*(nodes-temp_var1)/processor_num:(i+1)*(nodes-temp_var1)/processor_num], fmt='%20.15f') #Uncomment for codes which read XVOL
        print "%s %%" %(np.floor(np.float(i)/np.float(processor_num)*np.float(100)))

    num=np.transpose(np.array([[nodes+1-((processor_num-1)*(nodes-temp_var1)/processor_num+1)],[0],[0],[0],[0],[0],[0],[0],[0]]))
    lines=np.hstack((coor[(processor_num-1)*(nodes-temp_var1)/processor_num+1:nodes+1,:], dx_vec[(processor_num-1)*(nodes-temp_var1)/processor_num:nodes], dy_vec[(processor_num-1)*(nodes-temp_var1)/processor_num:nodes], dz_vec[(processor_num-1)*(nodes-temp_var1)/processor_num:nodes], mat_vec[(processor_num-1)*(nodes-temp_var1)/processor_num:nodes], np.transpose(np.array([vol[(processor_num-1)*(nodes-temp_var1)/processor_num:nodes]]))))
    np.savetxt('foreground%d.dat'%(processor_num, ), np.vstack((num,lines)), fmt='%7i %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %7i' '%20.15f')
    #np.savetxt('XVOL%d.dat'%(processor_num, ), vol[(processor_num-1)*(nodes-temp_var1)/processor_num:nodes], fmt='%20.15f')
    PetIGA().write("./Geometry.dat",S)
    print "Save Complete, %d nodes"%nodes

def save_PDGeometry(G, processor_num, name, units):
    # outputs input_coor.dat, XVOL.dat, and Geometry.dat from input coordinates, volumes, and NURBS obj; For Peridigm, use processor_num = 1 as the points
    # are distributed internally
    if(units=="NMS"):
        print "Saving input files...(NMS)"
        coor= G.coor[:,1:]
        vol=G.vols
        vol=np.around(vol, 15)
        nodes=coor.shape[0]-1
    elif(units=="mmNS"):
        print "Saving input files...(mmNS)"
        coor= G.coor[:,1:]*1000.0
        vol=G.vols*1.0e9
        vol=np.around(vol, 15)
        nodes=coor.shape[0]-1

    temp_var1=nodes%processor_num #number of remaining nodes added to last file
    for i in range(processor_num-1):
        lines=np.hstack((coor[i*(nodes-temp_var1)/processor_num+1:(i+1)*(nodes-temp_var1)/processor_num+1,:], np.ones(((i+1)*(nodes-temp_var1)/processor_num-i*(nodes-temp_var1)/processor_num, 1)), np.transpose(np.array([vol[i*(nodes-temp_var1)/processor_num:(i+1)*(nodes-temp_var1)/processor_num]]))))
        np.savetxt( name + '%d.txt'%(i+1, ), np.vstack(lines), fmt='%15.15f %15.15f %15.15f %1i %1.15e')
        print "%s %%" %(np.floor(np.float(i)/np.float(processor_num)*np.float(100)))

    lines=np.hstack((coor[(processor_num-1)*(nodes-temp_var1)/processor_num+1:nodes+1,:], np.ones((nodes-(processor_num-1)*(nodes-temp_var1)/processor_num, 1)), np.transpose(np.array([vol[(processor_num-1)*(nodes-temp_var1)/processor_num:nodes]]))))
    np.savetxt(name + '%d.txt'%(processor_num, ), np.vstack(lines), fmt='%15.15f %15.15f %15.15f %1i %1.15e')
    print "Save Complete, %d nodes"%nodes


def subt_rect_dom_fg(G, p1, p2, p3, p4):
#Subtracts rectangle primitive domain from foreground domain given 4 corner point vectors as nparray p1 p2 p3 p4
    xyz=G.coor[1:,1:-1]
    temp=G.coor[1:,:]
    p = path.Path([p1,p2,p3,p4])
    bools=p.contains_points(xyz)
    bools=np.invert(p.contains_points(xyz))
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    G.mat = G.mat[bools]
    return(G)

def in_cube(xyz, xmin, xmax, ymin, ymax, zmin, zmax):
    out = np.logical_and(np.logical_and(np.logical_and(xyz[:, 0]>=xmin, xyz[:, 0]<=xmax), np.logical_and(xyz[:, 1]>=ymin, xyz[:, 1]<=ymax)), np.logical_and(xyz[:, 2]>=zmin, xyz[:, 2]<=zmax))
    return(out)

def subt_cubic_dom_fg(G, xmin, xmax, ymin, ymax, zmin, zmax):
    #Subtracts cubic primitive domain from foreground domain given 4 corner point vectors as nparray p1 p2 p3 p4
    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]
    bools=in_cube(xyz, xmin, xmax, ymin, ymax, zmin, zmax)
    bools=np.invert(bools)
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    G.mat = G.mat[bools]
    return(G)


def in_circle(xyz, center, radius):
    # Determine if point(s) is(are) within the 2d circular region defined by center and radius
    out = np.sqrt(np.sum(np.multiply(np.array(xyz)-center,np.array(xyz)-center), axis=1))>=radius
    return(out)

def in_cone(xyz, center, height, r_base, ax):
    # Determine if point(s) is(are) within a cone defined by the
    # height, and the radius of the base and the axis
    X = np.array(xyz[:,0])
    x0 = center[0];
    Y = np.array(xyz[:,1])
    y0 = center[1]
    Z = np.array(xyz[:,2])
    z0 = center[2]

    if(ax == 0 ):
        #cone oriented along the x-axis
        YZ = np.column_stack((Y,Z))
        cent = np.column_stack((y0,z0))
        h = xyz[:, 0];
        r_crit = r_base-h*r_base/height;
        out = np.sqrt(np.sum(np.multiply(np.array(YZ)-cent,np.array(YZ)-cent), axis=1))<=r_crit
    elif(ax == 1):
        XZ = np.column_stack((X,Z))
        cent = np.column_stack((x0, z0))
        h = xyz[:, 1];
        r_crit = r_base-h*r_base/height;
        out = np.sqrt(np.sum(np.multiply(np.array(XZ)-cent,np.array(XZ)-cent), axis=1))<=r_crit
    elif(ax == 2):
        XY = np.column_stack((X,Y))
        cent = np.column_stack((x0, y0))
        h = xyz[:, 2];
        r_crit = r_base-h*r_base/height;
        out = np.sqrt(np.sum(np.multiply(np.array(XY)-cent,np.array(XY)-cent), axis=1))<=r_crit
    else:
        print("Error, axis > 2");
    return(out)

def in_cone2(xyz, center, height, r_base, r_top, ax):
    # Determine if point(s) is(are) within a cone defined by the
    # height, and the radius of the base and the axis
    X = np.array(xyz[:,0])
    x0 = center[0];
    Y = np.array(xyz[:,1])
    y0 = center[1]
    Z = np.array(xyz[:,2])
    z0 = center[2]

    if(ax == 0 ):
        #cone oriented along the x-axis
        YZ = np.column_stack((Y,Z))
        cent = np.column_stack((y0,z0))
        h = xyz[:, 0];
        r_crit = r_base-h*(r_base-r_top)/height;
        out = np.sqrt(np.sum(np.multiply(np.array(YZ)-cent,np.array(YZ)-cent), axis=1))<=r_crit
    elif(ax == 1):
        XZ = np.column_stack((X,Z))
        cent = np.column_stack((x0, z0))
        h = xyz[:, 1];
        r_crit = r_base-h*(r_base-r_top)/height;
        out = np.sqrt(np.sum(np.multiply(np.array(XZ)-cent,np.array(XZ)-cent), axis=1))<=r_crit
    elif(ax == 2):
        XY = np.column_stack((X,Y))
        cent = np.column_stack((x0, y0))
        h = xyz[:, 2];
        r_crit = r_base-h*(r_base-r_top)/height;
        out = np.sqrt(np.sum(np.multiply(np.array(XY)-cent,np.array(XY)-cent), axis=1))<=r_crit
    else:
        print("Error, axis > 2");
    return(out)

def in_cylinder(xyz, center, radius, ax):
    # Determine if point(s) is(are) within the 2d circular region defined by center and radius
    X = np.array(xyz[:,0])
    x0 = center[0];
    Y = np.array(xyz[:,1])
    y0 = center[1]
    Z = np.array(xyz[:,2])
    z0 = center[2]
    if(ax==0):
        YZ = np.column_stack((Y,Z))
        cent = np.column_stack((y0,z0))
        out = np.sqrt(np.sum(np.multiply(np.array(YZ)-cent,np.array(YZ)-cent), axis=1))<=radius
    if(ax==1):
        XZ = np.column_stack((X,Z))
        cent = np.column_stack((x0, z0))
        out = np.sqrt(np.sum(np.multiply(np.array(XZ)-cent,np.array(XZ)-cent), axis=1))<=radius
    if(ax==2):
        XY = np.column_stack((X,Y))
        cent = np.column_stack((x0, y0))
        out = np.sqrt(np.sum(np.multiply(np.array(XY)-cent,np.array(XY)-cent), axis=1))<=radius
    return(out)

def in_sphere(xyz, center, radius):
    # Determine if point(s) is(are) within the 3d spherical region defined by center and radius
    out = np.sqrt(np.sum(np.multiply(np.array(xyz)-center,np.array(xyz)-center), axis=1))<=(radius);
    return(out)

def subt_circular_domain(G, center, radius):
    # Subtracts cirucular domain from foreground given center and radius
    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]
    bools=in_circle(xyz, center, radius)
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    G.mat = G.mat[bools]
    return(G)

def Sphere_Mask(G, center, radius):
    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]
    bools=in_sphere(xyz, center, radius)
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    G.mat = G.mat[bools]
    return(G)

def Cylinder_Mask(G, center, radius, axis = 2):
    # Make a cylindrical Lagrangian geometry. By default, make cyl over the center of a
    # cube in the z-direction

    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]

    bools = in_cylinder(xyz, center, radius, axis)
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    G.mat = G.mat[bools]
    return(G)

def Conic_Mask(G, center, height, r_base, axis = 2, keep = 1):
    # Make a conic Lagrangian geometry. By default, make cone over the center of a
    # cube in the z-direction by default

    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]

    bools = in_cone(xyz, center, height, r_base, axis)
    if keep == 0:
        bools = not bools

    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    G.mat = G.mat[bools]
    return(G)

def Conic_Mask2(G, center, height, r_base, r_top, axis = 2, keep = 1):
    # Make a conic Lagrangian geometry. By default, make cone over the center of a
    # cube in the z-direction by default

    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]

    bools = in_cone2(xyz, center, height, r_base, r_top, axis)
    if(keep == 0):
        bools = np.invert(bools)

    nodenum=np.sum(bools)

    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    G.mat = G.mat[bools]
    return(G)

def generate_FGCone_cylindrical2(origin, r1, r2, r_min, height, xyz_numpts_vec, material, filled, axis = 2, min_pts = 16):
        # xyz numpts vec is [layers, particles per ring, num_divZ]
    h=np.zeros(xyz_numpts_vec[2])
    dh=np.double(height-origin[2])/(np.double(xyz_numpts_vec[2])-1)

    for i in range(xyz_numpts_vec[2]):
        if i>0: h[i]=h[i-1]+dh
        if i==0: h[i]=origin[2]
        hr=h[i]+dh/2.0
        hl=h[i]-dh/2.0
        if hr>height: hr=height 
        if hl<0: hl=origin[2]
        r_crit = r1-h[i]*(r1-r2)/height;
        if(i ==0):
            G = GetConeLayer2(origin, h[i], hr, hl, r_crit, r_min, xyz_numpts_vec, material, filled, min_pts)
        elif(i>0):
            G1 = GetConeLayer2(origin, h[i], hr, hl, r_crit, r_min, xyz_numpts_vec, material, filled, min_pts)
            G = sf.fg_superpose(G, G1)
    print "Complete"
    return(G)

def GetConeLayer2(origin, h, hr, hl, r_crit, r_min, xyz_numpts_vec, material, filled, min_pts):
    layers = xyz_numpts_vec[0]
    r = np.zeros(xyz_numpts_vec[0])
    dr = (r_crit-r_min)/(layers-1)

    for j in range(layers):
        if j>0: r[j]=r[j-1]+dr
        if j==0: r[j]=r_min
        rr=r[j]+dr/2.0
        rl=r[j]-dr/2.0
        if rr>r_crit: rr=r_crit
        if rl<origin[0]: rl=r_min*(filled%1)
        if(j == 0):
            G = getSegmentedCircle(rl, rr, r[j], r_crit, h, hr, hl, xyz_numpts_vec, material, min_pts)
        elif(j>0):
            G1 = getSegmentedCircle(rl, rr, r[j], r_crit, h, hr, hl, xyz_numpts_vec, material, min_pts)
            G = sf.fg_superpose(G, G1)
    return(G)

def getSegmentedCircle(rl, rr, r, r_crit, h, hr, hl, xyz_numpts_vec, material, min_pts = 16):
    arc = 2.0*np.pi/xyz_numpts_vec[1]
    arcL = arc*r_crit
    node_num = int(2*np.pi*r/arcL)

    if(arcL>2*np.pi*r):
        node_num = min_pts
    arc = 2.0*np.pi/node_num
    dtheta = arc

    dr_vec=np.zeros(node_num)
    dtheta_vec=np.zeros(node_num)
    dh_vec=np.zeros(node_num)
    theta=np.zeros(node_num)
    rr_vec=np.zeros(node_num)
    rl_vec=np.zeros(node_num)
    thetar_vec=np.zeros(node_num)
    thetal_vec=np.zeros(node_num)
    hr_vec=np.zeros(node_num)
    hl_vec=np.zeros(node_num)

    mat_vec=np.zeros(node_num)

    coor = np.zeros((int(node_num),4))
    vol = np.zeros(node_num)
    node = 0;
    for i in range(node_num):
        if i==0: theta[i]= 0.0
        if i>0:  theta[i] = theta[i-1] + dtheta
        thetar=theta[i] + dtheta/2.0
        thetal=theta[i] - dtheta/2.0
        [x, y] = pol2cart(r, theta[i])
        coor[node,:]=[node+1, x, y ,h]
        vol[node] = (hr-hl)*(rr**2-rl**2)*arc/2.0
        rr_vec[node]=rr
        rl_vec[node]=rl
        thetar_vec[node]=thetar
        thetal_vec[node]=thetal
        hr_vec[node]=hr
        hl_vec[node]=hl
        dr_vec[node]=rr-rl
        dtheta_vec[node]=thetar-thetal
        dh_vec[node]=hr-hl
        mat_vec[node]=material
        node = node+1
    coor=np.vstack(([node,0,0,0], coor))
    G=sf.foreground(coor, vol, xyz_numpts_vec, dr_vec, dtheta_vec, dh_vec, rr_vec, rl_vec, thetar_vec, thetal_vec, hr_vec, hl_vec, mat_vec)
    return(G)
