from igakit.cad import *
from igakit.cad import extrude
from igakit.cad import join
from igakit.igalib import *
from igakit.plot import *
from igakit.io import PetIGA,VTK
from numpy import linspace
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import mesh_functions as ms


class foreground():
    # Contains coor, volumes, numpts, dx,dy,dz and right and left boundaries for volumes
    def __init__(self, coor, vols, xyx_numpts_vec, dx_vec, dy_vec, dh_vec, xr, xl, yr, yl, hr, hl, mat):
        self.coor = coor
        self.vols = vols
        self.xyz_numpts_vec = xyx_numpts_vec
        self.dx_vec = dx_vec
        self.dy_vec = dy_vec
        self.dh_vec = dh_vec
        self.xr = xr
        self.xl = xl
        self.yr = yr
        self.yl = yl
        self.hr = hr
        self.hl = hl
        self.mat= mat # Material flag 


def lin_stretch(axis, lin_center, poly_order, foreground_obj):
    #given axis (0,1,2), origin and polynomial order for stretching (>1 for center and <1 for edge) and foreground instance, stretches foreground mesh in 1-direction toward or away from lin_center (scalar)
    
    xtilde=foreground_obj.coor[:,axis+1]-lin_center
    mx=np.amax(xtilde)
    mn=np.absolute(np.amin(xtilde))
    x_p_t=xtilde[xtilde>0]/mx
    x_n_t=xtilde[xtilde<0]/mn
    xtilde[xtilde>0]=(np.sign(x_p_t)*np.absolute(x_p_t)**poly_order)*mx
    xtilde[xtilde<0]=(np.sign(x_n_t)*np.absolute(x_n_t)**poly_order)*mn
    xtilde=xtilde+lin_center
    foreground_obj.coor[:,axis+1]=xtilde
    if axis==0:
        vol_tilde=np.divide(foreground_obj.vols,foreground_obj.dx_vec)

        xtilder = foreground_obj.xr-lin_center
        xtildel = foreground_obj.xl-lin_center
        
        x_p_tr=xtilder[xtilder>0]/mx
        x_n_tr=xtilder[xtilder<0]/mn
        x_p_tl=xtildel[xtildel>0]/mx
        x_n_tl=xtildel[xtildel<0]/mn
        
        xtilder[xtilder>0]=(np.sign(x_p_tr)*np.absolute(x_p_tr)**poly_order)*mx
        xtilder[xtilder<0]=(np.sign(x_n_tr)*np.absolute(x_n_tr)**poly_order)*mn
        xtildel[xtildel>0]=(np.sign(x_p_tl)*np.absolute(x_p_tl)**poly_order)*mx
        xtildel[xtildel<0]=(np.sign(x_n_tl)*np.absolute(x_n_tl)**poly_order)*mn

        xtilder=xtilder+lin_center #This reassignment step is necessary to allow for stretching to edges
        xtildel=xtildel+lin_center
        
        dxhat=xtilder-xtildel
        volhat=vol_tilde*(dxhat)
        G=foreground(foreground_obj.coor, volhat, foreground_obj.xyz_numpts_vec, dxhat, foreground_obj.dy_vec, foreground_obj.dh_vec, xtilder, xtildel, foreground_obj.yr, foreground_obj.yl, foreground_obj.hr, foreground_obj.hl, foreground_obj.mat)
        
        ## Test statements:
        assert np.absolute(np.sum(volhat)-np.sum(foreground_obj.vols))<10**(-12)#, "Volume error"
        ##
    elif axis==1:
        vol_tilde=np.divide(foreground_obj.vols,foreground_obj.dy_vec)
        
        ytilder = foreground_obj.yr-lin_center
        ytildel = foreground_obj.yl-lin_center
        
        y_p_tr=ytilder[ytilder>0]/mx
        y_n_tr=ytilder[ytilder<0]/mn
        y_p_tl=ytildel[ytildel>0]/mx
        y_n_tl=ytildel[ytildel<0]/mn
        
        ytilder[ytilder>0]=(np.sign(y_p_tr)*np.absolute(y_p_tr)**poly_order)*mx
        ytilder[ytilder<0]=(np.sign(y_n_tr)*np.absolute(y_n_tr)**poly_order)*mn
        ytildel[ytildel>0]=(np.sign(y_p_tl)*np.absolute(y_p_tl)**poly_order)*mx
        ytildel[ytildel<0]=(np.sign(y_n_tl)*np.absolute(y_n_tl)**poly_order)*mn
        
        ytilder=ytilder+lin_center
        ytildel=ytildel+lin_center
        dyhat=ytilder-ytildel
        
        volhat=vol_tilde*(dyhat)
        
        G=foreground(foreground_obj.coor, volhat, foreground_obj.xyz_numpts_vec, foreground_obj.dx_vec, dyhat , foreground_obj.dh_vec, foreground_obj.xr, foreground_obj.xl, ytilder, ytildel, foreground_obj.hr, foreground_obj.hl, foreground_obj.mat)
        
        ## Test statements:
        assert np.absolute(np.sum(volhat)-np.sum(foreground_obj.vols))<10**(-12)#, "Volume error"
        ##
    elif axis==2:
        vol_tilde=np.divide(foreground_obj.vols,foreground_obj.dh_vec)
        
        htilder = foreground_obj.hr-lin_center
        htildel = foreground_obj.hl-lin_center
        
        h_p_tr=htilder[htilder>0]/mx
        h_n_tr=htilder[htilder<0]/mn
        h_p_tl=htildel[htildel>0]/mx
        h_n_tl=htildel[htildel<0]/mn
        
        htilder[htilder>0]=(np.sign(h_p_tr)*np.absolute(h_p_tr)**poly_order)*mx
        htilder[htilder<0]=(np.sign(h_n_tr)*np.absolute(h_n_tr)**poly_order)*mn
        htildel[htildel>0]=(np.sign(h_p_tl)*np.absolute(h_p_tl)**poly_order)*mx
        htildel[htildel<0]=(np.sign(h_n_tl)*np.absolute(h_n_tl)**poly_order)*mn
        
        htilder=htilder+lin_center
        htildel=htildel+lin_center
        dhhat=htilder-htildel
        volhat=vol_tilde*(dhhat)
        
        G=foreground(foreground_obj.coor, volhat, foreground_obj.xyz_numpts_vec, foreground_obj.dx_vec, foreground_obj.dy_vec , dhhat, foreground_obj.xr, foreground_obj.xl, foreground_obj.yr, foreground_obj.yl, htilder, htildel, foreground_obj.mat)
        
        ## Test statements:
        assert np.absolute(np.sum(volhat)-np.sum(foreground_obj.vols))<10**(-12)#, "Volume error"
        ##
    return(G)

def cent_stretch_fg(xyz, poly_order, foreground_obj):
    #Foreground stretching function toward a coordinate point xyz given np.array([x0,y0,z0]), polynomial order as np.array and foreground instance
    G=lin_stretch(0, xyz[0], poly_order[0], foreground_obj)
    G=lin_stretch(1, xyz[1], poly_order[1],  G)
    if xyz[2]>0: G=lin_stretch(2, xyz[2], poly_order[2], G)
    return(G)

def cent_stretch_bg(C1, C2, t, num_div_vec, xyz, poly_order, n=1):
    #Stretched background generation function, given line1 line2 thickness numknots vector xyz origin of stretch and polynomial order as np.array. Returns NURBS object stretched mesh.
    S=ruled(C1, C2)
    S.elevate(1,n)
    S.elevate(0,n)
    O=C1.points[0,:]
    if t>0:
        S=extrude(S, t, 2)
        S.elevate(2,n)

    x=np.linspace(0,1,num_div_vec[0])
    y=np.linspace(0,1,num_div_vec[1])
    if num_div_vec[2]>1: z=np.linspace(0,1,num_div_vec[2])

    xtilde=x-(xyz[0]-O[0])/(C1.points[1,0]-O[0])
    x_p_t=xtilde[xtilde>0]/np.amax(xtilde)
    x_n_t=xtilde[xtilde<0]/np.amin(xtilde)
    x_p_t=(np.sign(x_p_t)*np.absolute(x_p_t)**poly_order[0])*np.amax(xtilde)
    x_n_t=(np.sign(x_n_t)*np.absolute(x_n_t)**poly_order[0])*np.amin(xtilde)
    xtilde[xtilde>0]=x_p_t
    xtilde[xtilde<0]=x_n_t
    xhat=xtilde+(xyz[0]-O[0])/(C1.points[1,0]-O[0])

    ytilde=y-(xyz[1]-O[1])/(C2.points[1,1]-O[1])
    y_p_t=ytilde[ytilde>0]/np.amax(ytilde)
    y_n_t=ytilde[ytilde<0]/np.amin(ytilde)
    y_p_t=(np.sign(y_p_t)*np.absolute(y_p_t)**poly_order[1])*np.amax(ytilde)
    y_n_t=(np.sign(y_n_t)*np.absolute(y_n_t)**poly_order[1])*np.amin(ytilde)
    ytilde[ytilde>0]=y_p_t
    ytilde[ytilde<0]=y_n_t
    yhat=ytilde+(xyz[1]-O[1])/(C2.points[1,1]-O[1])



    if num_div_vec[2]>1:
        ztilde=z-(xyz[2]-O[2])/(t)
        z_p_t=ztilde[ztilde>0]/np.amax(ztilde)
        z_n_t=ztilde[ztilde<0]/np.amin(ztilde)
        z_p_t=(np.sign(z_p_t)*np.absolute(z_p_t)**poly_order[2])*np.amax(ztilde)
        z_n_t=(np.sign(z_n_t)*np.absolute(z_n_t)**poly_order[2])*np.amin(ztilde)
        ztilde[ztilde>0]=z_p_t
        ztilde[ztilde<0]=z_n_t
        zhat=ztilde+xyz[2]/t
        S.refine(2,zhat[1:-1])
    S.refine(0,xhat[1:-1])
    S.refine(1,yhat[1:-1])
    
    #When we make stretched grids, information about the element size is required to properly compute Lagrangian point ranks:
    
    #lines = np.column_stack((np.transpose(xhat), np.transpose(yhat), np.transpose(zhat)))
    #lines = np.vstack((np.array([xhat.shape[0], yhat.shape[0], zhat.shape[0]]), lines))
    #print(lines)
    #np.savetxt('ElementInfo.dat', lines, fmt='%15.10f %15.10f %15.10f')
    return(S)


def translate(G,x,y,z):
    G.coor[1:,1]=G.coor[1:,1]+x*np.ones(G.coor[1:,1].shape)
    G.coor[1:,2]=G.coor[1:,2]+y*np.ones(G.coor[1:,2].shape)
    G.coor[1:,3]=G.coor[1:,3]+z*np.ones(G.coor[1:,3].shape)
    G.xr = G.xr+x*np.ones(G.coor[1:,1].shape)
    G.xl = G.xl+x*np.ones(G.coor[1:,1].shape)
    G.yr = G.yr+y*np.ones(G.coor[1:,2].shape)
    G.yl = G.yl+y*np.ones(G.coor[1:,2].shape)
    G.hr = G.hr+z*np.ones(G.coor[1:,3].shape)
    G.hl = G.hl+z*np.ones(G.coor[1:,3].shape)
    fg_output=foreground(G.coor, G.vols, G.xyz_numpts_vec, G.dx_vec, G.dy_vec, G.dh_vec, G.xr, G.xl, G.yr, G.yl, G.hr, G.hl, G.mat)
    return(fg_output)

def merge(G1,G2, axis, lin_center):
    assert axis<=2
    inputvol=np.sum(G1.vols)+np.sum(G2.vols)
    coor1=G1.coor[1:,:]
    coor2=G2.coor[1:,:]
    ind1=np.invert(np.absolute(coor1[:,axis+1]-lin_center)<=10**(-12))
    ind2=np.invert(np.absolute(coor2[:,axis+1]-lin_center)<=10**(-12))
    assert ind1[ind1==0].size==ind2[ind2==0].size, "All particle discretizations other than the discritization along merge axis must match"
    coor1=coor1[ind1, :]
    coor=np.vstack([coor1,coor2])
    
    for i in range(coor.shape[0]):
        coor[i,0]=i+1
    coor=np.vstack(([coor.shape[0],0,0,0], coor))
    if axis==0:
        xl1=G1.xl
        xl2=G2.xl
        xr1=G1.xr
        xr2=G2.xr
        vol1=np.divide(G1.vols,G1.dx_vec)
        vol2=np.divide(G2.vols,G2.dx_vec)
        
        xl2[np.invert(ind2)]=xl1[np.invert(ind1)]
        xl1=xl1[ind1]
        
        xl=np.hstack([xl1,xl2])
        xr1=xr1[ind1]
        xr=np.hstack([xr1,xr2])
        vol1=vol1[ind1]
        vol=np.hstack([vol1, vol2])
        
        # Update attributes and construct output object:
        
        dxhat=xr-xl
        vol=np.multiply(vol,dxhat)
        xyz_numpts_vec=G1.xyz_numpts_vec+G2.xyz_numpts_vec-np.array([1,0,0])
        dx_vec=dxhat
        dy_vec=np.hstack([G1.dy_vec[ind1],G2.dy_vec])
        dh_vec=np.hstack([G1.dh_vec[ind1],G2.dh_vec])
        yr=np.hstack([G1.yr[ind1],G2.yr])
        yl=np.hstack([G1.yl[ind1],G2.yl])
        hr=np.hstack([G1.hr[ind1],G2.hr])
        hl=np.hstack([G1.hl[ind1],G2.hl])
        mat=np.hstack([G1.mat[ind1],G2.mat])
        fg_output=foreground(coor, vol, xyz_numpts_vec, dx_vec, dy_vec, dh_vec, xr, xl, yr, yl, hr, hl, mat)
        #Volume preservation test
        print np.absolute(np.sum(vol)-inputvol)<=10**(-12)#, "Volume error"
        #
        return(fg_output)
    
    elif axis==1:
        yl1=G1.yl
        yl2=G2.yl
        yr1=G1.yr
        yr2=G2.yr
        vol1=np.divide(G1.vols,G1.dy_vec)
        vol2=np.divide(G2.vols,G2.dy_vec)
        
        yl2[np.invert(ind2)]=yl1[np.invert(ind1)]
        yl1=yl1[ind1]
        
        yl=np.hstack([yl1,yl2])
        yr1=yr1[ind1]
        yr=np.hstack([yr1,yr2])
        vol1=vol1[ind1]
        vol=np.hstack([vol1, vol2])
        
        # Update attributes and construct output object:
        
        dyhat=yr-yl
        vol=np.multiply(vol,dyhat)
        xyz_numpts_vec=G1.xyz_numpts_vec+G2.xyz_numpts_vec-np.array([0,1,0])
        dy_vec=dyhat
        dx_vec=np.hstack([G1.dx_vec[ind1],G2.dx_vec])
        dh_vec=np.hstack([G1.dh_vec[ind1],G2.dh_vec])
        xr=np.hstack([G1.xr[ind1],G2.xr])
        xl=np.hstack([G1.xl[ind1],G2.xl])
        hr=np.hstack([G1.hr[ind1],G2.hr])
        hl=np.hstack([G1.hl[ind1],G2.hl])
        mat=np.hstack([G1.mat[ind1],G2.mat])
        fg_output=foreground(coor, vol, xyz_numpts_vec, dx_vec, dy_vec, dh_vec, xr, xl, yr, yl, hr, hl, mat)
        #Volume preservation test
        print np.absolute(np.sum(vol)-inputvol)<=10**(-12)
        #
        return(fg_output)
    
    elif axis==2:
        hl1=G1.hl
        hl2=G2.hl
        hr1=G1.hr
        hr2=G2.hr
        vol1=np.divide(G1.vols,G1.dh_vec)
        vol2=np.divide(G2.vols,G2.dh_vec)
        
        hl2[np.invert(ind2)]=hl1[np.invert(ind1)]
        hl1=hl1[ind1]
        
        hl=np.hstack([hl1,hl2])
        hr1=hr1[ind1]
        hr=np.hstack([hr1,hr2])
        vol1=vol1[ind1]
        vol=np.hstack([vol1, vol2])
        
        # Update attributes and construct output object:
        
        dhhat=hr-hl
        vol=np.multiply(vol,dhhat)
        xyz_numpts_vec=G1.xyz_numpts_vec+G2.xyz_numpts_vec-np.array([0,0,1])
        dh_vec=dhhat
        dx_vec=np.hstack([G1.dx_vec[ind1],G2.dx_vec])
        dy_vec=np.hstack([G1.dy_vec[ind1],G2.dy_vec])
        xr=np.hstack([G1.xr[ind1],G2.xr])
        xl=np.hstack([G1.xl[ind1],G2.xl])
        yr=np.hstack([G1.yr[ind1],G2.yr])
        yl=np.hstack([G1.yl[ind1],G2.yl])
        mat=np.hstack([G1.mat[ind1],G2.mat])
        fg_output=foreground(coor, vol, xyz_numpts_vec, dx_vec, dy_vec, dh_vec, xr, xl, yr, yl, hr, hl, mat)
        #Volume preservation test
        assert np.absolute(np.sum(vol)-inputvol)<=10**(-12)
        #
        return(fg_output)

def merge_bg(S1, S2, axis):
    # Wraps igakit.cad function join to merge two background domains
    S1.unclamp(axis)
    S2.unclamp(axis)
    out=join(S1, S2, axis)
    out.remap(axis,0,1)
    out.clamp(axis)
    return(out)
    
    
def fg_superpose(G1,G2):
    #Superimposes two separate foreground domains for computation. Stacks and re-indexes only
    node=np.int64(G1.coor[0,0]+G2.coor[0,0])
    coor=np.vstack((G1.coor[1:,:],G2.coor[1:,:]))
    for k in range(node):
        coor[k, 0]=k+1
    coor=np.vstack(([node,0,0,0], coor))
    
    vols=np.array(np.hstack((G1.vols,G2.vols)))
    xyz_numpts_vec=np.hstack((G1.xyz_numpts_vec,G2.xyz_numpts_vec))
    dx_vec=np.hstack((G1.dx_vec,G2.dx_vec))
    dy_vec=np.hstack((G1.dy_vec,G2.dy_vec))
    dh_vec=np.hstack((G1.dh_vec,G2.dh_vec))
    xr=np.array(np.hstack((G1.xr,G2.xr)))
    xl=np.array(np.hstack((G1.xl,G2.xl)))
    yr=np.array(np.hstack((G1.yr,G2.yr)))
    yl=np.array(np.hstack((G1.yl,G2.yl)))
    hr=np.array(np.hstack((G1.hr,G2.hr)))
    hl=np.array(np.hstack((G1.hl,G2.hl)))
    mat=np.array(np.hstack((G1.mat,G2.mat)))
    G=foreground(coor, vols, xyz_numpts_vec, dx_vec, dy_vec, dh_vec, xr, xl, yr, yl, hr, hl, mat)
    return(G)

def nonlinear_parameterization(S):
    # Takes as input a NURBS surface object, S, and uniformly spaces the control points to yield a nonlinearly parametrized NURBS surface (see IGA_vibrations Ch. 5). n is order of continuity of the nurbs object
    controlpts=S.control
    #print(controlpts.shape)
    xcpts=np.linspace(controlpts[0,0, 0], controlpts[-1,0, 0], num=controlpts.shape[0]);
    #print(xcpts)
    ycpts=np.linspace(controlpts[0,0, 1], controlpts[0,-1, 1], num=controlpts.shape[1]);
    #print(ycpts)
    zcpts=np.linspace(controlpts[1, 1, 0], controlpts[1, 1,-1], num=controlpts.shape[2]);
    #print(zcpts)
    #for i in range(controlpts.shape[1]):
    #    for j in range(controlpts.shape[2]):
    #        if j==0:
    #            controlpts[:,i, j]=xcpts
    #        if j==1:
    #            controlpts[:,i,j]=ycpts[i]*np.ones(controlpts[:,i,j].size)
    #        if j==2:
    #            controlpts[:,i,j]=np.zeros(controlpts[:,i,j].size)
    #        if j==3:
    #            controlpts[:,i,j]=np.ones(controlpts[:,i,j].size)
            #print(controlpts[:,i, j])
    #for i in range(controlpts.shape[0]):
    #    for j in range(controlpts.shape[2]):
    #            if j==0:
    #                controlpts[i, :, j]=xcpts[i]*np.ones(controlpts[i,:,j].size)
    #            if j==1:
    #                controlpts[i, :, j]=ycpts
    #            if j==2:
    #                controlpts[i, :, j]=np.zeros(controlpts[i,:,j].size)
    #             if j==3:
    #                controlpts[i, :, j]=np.ones(controlpts[i,:,j].size)

    for i in range(controlpts.shape[0]):
       for j in range(controlpts.shape[1]):
            #print(controlpts[i, j, :])
            controlpts[i, j, :]=np.array([xcpts[i], ycpts[j], 0.0, 1.0])
    
    S1=NURBS(S.knots, control=controlpts, fields=S.fields)
    print(S1.control)
    return(S1)
