    def nonlinear_parameterization(S):
    # Takes as input a NURBS surface object, S, and uniformly spaces the control points to yield a nonlinearly parametrized NURBS surface (see IGA_vibrations Ch. 5)
        controlpts=S.control

        for i in range(controlpts.shape[1]):
            for j in range(controlpts.shape[2]):
                controlpts[:, i, j]=np.linspace(controlpts[0, i, j], controlpts[-1 , i , j], num=controlpts.shape[0])
        
        for i in range(controlpts.shape[0]):
            for j in range(controlpts.shape[2]):
                controlpts[i, :, j]=np.linspace(controlpts[i, 0, j], controlpts[i, -1, j], num=controlpts.shape[1])
        
#for i in range(controlpts.shape[0]):
   # for j in range(controlpts.shape[1]):
       # controlpts[i, j, :]=np.linspace(controlpts[i, j, 0], controlpts[i, j,-1], num=controlpts.shape[2])

    S1=NURBS(S.knots, control=controlpts, fields=S.fields, weights=None)
    print(S1.control)
    return(S1)

