#------------------------------------------------------------------------------------------------------------
#
#              _____ ______ _____             __     _____  ______  
#             / ____|  ____|  __ \           / _|   |  __ \|  ____| 
#            | |    | |__  | |  | |     ___ | |_    | |__) | |__    
#            | |    |  __| | |  | |    / _ \|  _|   |  _  /|  __|   
#            | |____| |    | |__| |   | (_) | |     | | \ \| |      
#             \_____|_|    |_____/     \___/|_|     |_|  \_\_|      
#                                                        
#                                                        
#------------------------------------------------------------------------------------------------------------
#
#   Technische Universiteit Delft - TUDelft (2022)
#
#   Master of Science in Chemical Engineering
#
#   This code was developed and tested by Elia Ferretti
#
#   You can redistribute the code and/or modify it
#   Whenever the code is used to produce any publication or document,
#   reference to this work and author should be reported
#   No warranty of fitness for a particular purpose is offered
#   The user must assume the entire risk of using this code
#
#------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------
#                                          --- FUNCTIONS ---
#------------------------------------------------------------------------------------------------------------


def preAllocation(nx,ny):
    p = np.zeros((nx+2,ny+2))
    u = np.zeros((nx+1,ny+2))
    v = np.zeros((nx+2,ny+1))
    ut = np.zeros((nx+1,ny+2))
    vt = np.zeros((nx+2,ny+1))
    pp = np.zeros((nx+1,ny+1))
    uu = np.zeros((nx+1,ny+1))
    vv = np.zeros((nx+1,ny+1))
    CA = np.zeros((nx+2,ny+2))
    CB = np.zeros((nx+2,ny+2))
    CC = np.zeros((nx+2,ny+2))
    CD = np.zeros((nx+2,ny+2))
    CE = np.zeros((nx+2,ny+2))
    CAstar = np.zeros((nx+2,ny+2))
    CBstar = np.zeros((nx+2,ny+2))
    CCstar = np.zeros((nx+2,ny+2))
    CDstar = np.zeros((nx+2,ny+2))
    CEstar = np.zeros((nx+2,ny+2))
    ca = np.zeros((nx+1,ny+1))
    cb = np.zeros((nx+1,ny+1))
    cc = np.zeros((nx+1,ny+1))
    cd = np.zeros((nx+1,ny+1))
    ce = np.zeros((nx+1,ny+1))
    den = np.zeros((nx+2,ny+2))
    position = np.zeros((4,4))
    return p,u,v,ut,vt,pp,uu,vv,CA,CB,CC,CD,CE,CAstar,CBstar,CCstar,CDstar,CEstar,ca,cb,cc,cd,ce,den,position

def initializeDenominator(nx,ny,hx,hy):
    #dimension of den are: (nx+2)*(ny+2)
    den[:,:] = 2*(hx**2+hy**2);
    #horizontal sides
    den[1:int(nx*2/5+1),1] = hx**2+2*hy**2;                         #1
    den[1:int(nx*2/5+1),ny] = hx**2+2*hy**2;                        #2
    den[int(nx*2/5+1):int(nx+1),int(ny/4+1)] = hx**2+2*hy**2;       #3
    den[int(nx*2/5+1):int(nx+1),int(ny*3/4)] = hx**2+2*hy**2;       #4
    den[1:int(nx/5+1),int(ny/10)] = hx**2+2*hy**2;                  #5
    den[1:int(nx/5+1),int(ny/4+1)] = hx**2+2*hy**2;                 #6
    den[1:int(nx/5+1),int(ny*3/4)] = hx**2+2*hy**2;                 #7
    den[1:int(nx/5+1),int(ny*9/10+1)] = hx**2+2*hy**2;              #8
    #vertical sides
    den[int(nx/5+1),int(ny/10+1):int(ny/4+1)] = 2*hx**2+hy**2;      #9
    den[int(nx/5+1),int(ny*3/4+1):int(ny*9/10+1)] = 2*hx**2+hy**2;  #10
    den[int(nx*2/5),1:int(ny/4+1)] = 2*hx**2+hy**2;                 #11
    den[int(nx*2/5),int(ny*3/4+1):ny+1] = 2*hx**2+hy**2;            #12
    #corner
    den[int(nx*2/5),ny] = hx**2+hy**2;
    den[int(nx*2/5),1] = hx**2+hy**2;
    return den

def solvingDomain(nx,ny):
    solve = np.zeros((nx+2,ny+2));
    solve_u = np.zeros((nx+1,ny+2));
    solve_v = np.zeros((nx+2,ny+1));
    solve_r = np.zeros((nx+1,ny+1));
    
    #1-solve_u staggered grid
    for i in range(1,int(nx/5+1)):
        for j in range(1,int(ny/10+1)):
            solve_u[i,j] = 1;

    #2-solve_u staggered grid
    for i in range(1,int(nx/5+1)):
        for j in range(int(ny*9/10+1),ny+1):
            solve_u[i,j] = 1;

    #3-solve_u staggered grid
    for i in range(1,int(nx/5+1)):
        for j in range(int(ny/4+1),int(ny*3/4+1)):
            solve_u[i,j] = 1;

    #4-solve_u staggered grid
    for i in range(int(nx/5+1),int(nx*2/5)):
        for j in range(1,ny+1):
            solve_u[i,j] = 1;

    #5-solve_u staggered grid
    for i in range(int(nx*2/5),nx):
        for j in range(int(ny/4+1),int(ny*3/4+1)):
            solve_u[i,j] = 1;
            

    #1-solve_v staggered grid
    for i in range(1,int(nx/5+1)):
        for j in range(1,int(ny/10)):
            solve_v[i,j] = 1;

    #2-solve_v staggered grid
    for i in range(1,int(nx/5+1)):
        for j in range(int(ny*9/10+1),ny):
            solve_v[i,j] = 1;

    #3-solve_v staggered grid
    for i in range(1,int(nx/5+1)):
        for j in range(int(ny/4+1),int(ny*3/4)):
            solve_v[i,j] = 1;

    #4-solve_v staggered grid
    for i in range(int(nx/5+1),int(nx*2/5+1)):
        for j in range(1,ny):
            solve_v[i,j] = 1;

    #5-solve_v staggered grid
    for i in range(int(nx*2/5+1),int(nx+1)):
        for j in range(int(ny/4+1),int(ny*3/4)):
            solve_v[i,j] = 1;

    
    #1-solve
    for i in range(1,int(nx/5+1)):
        for j in range(1,int(ny/10+1)):
            solve[i,j] = 1;

    #2-solve
    for i in range(1,int(nx/5+1)):
        for j in range(int(ny*9/10+1),ny+1):
            solve[i,j] = 1;

    #3-solve
    for i in range(1,int(nx/5+1)):
        for j in range(int(ny/4+1),int(ny*3/4+1)):
            solve[i,j] = 1;

    #4-solve
    for i in range(int(nx/5+1),int(nx*2/5+1)):
        for j in range(1,ny+1):
            solve[i,j] = 1;
 
    #5-solve
    for i in range(int(nx*2/5+1),nx+1):
        for j in range(int(ny/4+1),int(ny*3/4+1)):
            solve[i,j] = 1;


    #1-reconstruct fields
    for i in range(0,int(nx/5+1)):
        for j in range(0,int(ny/10+1)):          
            solve_r[i,j] = 1;

    #2-reconstruct fields
    for i in range(0,int(nx/5+1)):
        for j in range(int(ny*9/10),ny+1):           
            solve_r[i,j] = 1;
  
    #3-reconstruct fields
    for i in range(0,int(nx/5+1)):
        for j in range(int(ny/4),int(ny*3/4+1)):           
            solve_r[i,j] = 1;

    #4-reconstruct fields
    for i in range(int(nx/5+1),int(nx*2/5+1)):
        for j in range(0,ny+1):          
            solve_r[i,j] = 1;

    #5-reconstruct fields
    for i in range(int(nx*2/5+1),nx+1):
        for j in range(int(ny/4),int(ny*3/4+1)):            
            solve_r[i,j] = 1;

    return solve,solve_u,solve_v,solve_r

def velocityBoundaryConditions(u,v,nx,ny,u_top,u_mid,u_bot): 
    #horizontal sides u-velocity
    u[0:int(nx*2/5),0] = -u[0:int(nx*2/5),1]                                                #1
    u[0:int(nx*2/5),ny+1] = -u[0:int(nx*2/5),ny]                                            #2                              
    u[int(nx*2/5):nx+1,int(ny/4)] = -u[int(nx*2/5):nx+1,int(ny/4+1)]                        #3
    u[int(nx*2/5):nx+1,int(ny*3/4+1)] = -u[int(nx*2/5):nx+1,int(ny*3/4)]                    #4                                           
    u[0:int(nx/5+1),int(ny/10+1)] = -u[0:int(nx/5+1),int(ny/10)]                            #5
    u[0:int(nx/5+1),int(ny/4)] = -u[0:int(nx/5+1),int(ny/4+1)]                              #6                           
    u[0:int(nx/5+1),int(ny*3/4+1)] = -u[0:int(nx/5+1),int(ny*3/4)]                          #7
    u[0:int(nx/5+1),int(ny*9/10)] = -u[0:int(nx/5+1),int(ny*9/10+1)]                        #8

    #horizontal sides v-velocity
    v[1:int(nx*2/5+1),0] = 0                                                                #1
    v[1:int(nx*2/5+1),int(ny)] = 0                                                          #2
    v[int(nx*2/5+1):nx+1,int(ny/4)] = 0                                                     #3
    v[int(nx*2/5+1):nx+1,int(ny*3/4)] = 0                                                   #4
    v[1:int(nx/5+1),int(ny/10)] = 0                                                         #5
    v[1:int(nx/5+1),int(ny/4)] = 0                                                          #6
    v[1:int(nx/5+1),int(ny*3/4)] = 0                                                        #7
    v[1:int(nx/5+1),int(ny*9/10)] = 0                                                       #8

    #vertical sides u-velocity
    u[int(nx/5),int(ny/10+1):int(ny/4+1)] = 0                                               #9
    u[int(nx/5),int(ny*3/4+1):int(ny*9/10+1)] = 0                                           #10
    u[int(nx*2/5),1:int(ny/4+1)] = 0                                                        #11
    u[int(nx*2/5),int(ny*3/4+1):ny+1] = 0                                                   #12

    #vertical sides v-velocity
    v[int(nx/5),int(ny/10+1):int(ny/4)] = -v[int(nx/5+1),int(ny/10+1):int(ny/4)]            #9
    v[int(nx/5),int(ny*3/4+1):int(ny*9/10)] = -v[int(nx/5+1),int(ny*3/4+1):int(ny*9/10)]    #10
    v[int(nx*2/5+1),0:int(ny/4)] = -v[int(nx*2/5),0:int(ny/4)]                              #11
    v[int(nx*2/5+1),int(ny*3/4+1):ny+1] = -v[int(nx*2/5),int(ny*3/4+1):ny+1]                #12

    #inlet u-velocity
    u[0,1:int(ny/10+1)] = u_bot                                                             #13
    u[0,int(ny/4+1):int(ny*3/4+1)] = u_mid                                                  #14
    u[0,int(ny*9/10+1):ny+1] = u_top                                                        #15    

    #inlet v-velocity
    v[0,0:int(ny/10+1)] = -v[1,0:int(ny/10+1)]                                              #13
    v[0,int(ny/4):int(ny*3/4+1)] = -v[1,int(ny/4):int(ny*3/4+1)]                            #14
    v[0,int(ny*9/10):ny+1] = -v[1,int(ny*9/10):ny+1]                                        #15

    #outlet u-velocity
    u[nx,int(ny/4+1):int(ny*3/4+1)] = u[nx-1,int(ny/4+1):int(ny*3/4+1)]                     #16   

    #outlet v-velocity
    v[nx+1,int(ny/4):int(ny*3/4+1)] = v[nx,int(ny/4):int(ny*3/4+1)]                         #16  
    return u,v

def temporaryVelocity(ut,vt,u,v,nx,ny,hx,hy,dt,nu,solve_u,solve_v):
    
    u = np.append(u,np.zeros((1,np.shape(u)[1])),axis=0)
    v = np.append(v,np.zeros((np.shape(v)[0],1)),axis=1)
    
    # u-velocity
    q_x = 0.5*hx*u/nu
    q_y = 0.5*hx*( v+np.roll(v,(-1,0),axis=(0,1))+np.roll(v,(0,+1),axis=(0,1))+np.roll(v,(-1,+1),axis=(0,1)) )/4/nu
    
    noTaylor_x = np.where(abs(q_x)>0.1)
    noTaylor_y = np.where(abs(q_y)>0.1)
    
    alpha_x = q_x/3 - q_x**3/45 + 2*q_x**5/945 - q_x**7/4725
    alpha_x[noTaylor_x] = 1/np.tanh(q_x[noTaylor_x]) - 1/q_x[noTaylor_x]
    
    alpha_y = q_y/3 - q_y**3/45 + 2*q_y**5/945 - q_y**7/4725
    alpha_y[noTaylor_y] = 1/np.tanh(q_y[noTaylor_y])-1/q_y[noTaylor_y]
    
    alpha = 0.25*( ( np.roll(u,(-1,0),axis=(0,1))+u )*( (1-alpha_x)*np.roll(u,(-1,0),axis=(0,1))+(1+alpha_x)*u ) - \
                   ( np.roll(u,(+1,0),axis=(0,1))+u )*( (1-alpha_x)*u + (1+alpha_x)*np.roll(u,(+1,0),axis=(0,1)) ) )/hx + \
            0.25*( ( np.roll(v,(-1,0),axis=(0,1))+v )*( (1-alpha_y)*np.roll(u,(0,-1),axis=(0,1))+(1+alpha_y)*u ) - \
                   ( np.roll(v,(-1,+1),axis=(0,1))+np.roll(v,(0,+1),axis=(0,1)) )*( (1-alpha_y)*u + (1+alpha_y)*np.roll(u,(0,+1),axis=(0,1)) ) )/hy 
                  
    delta = nu*( (np.roll(u,(-1,0),axis=(0,1))-2*u+np.roll(u,(+1,0),axis=(0,1)))/hx**2 + (np.roll(u,(0,-1),axis=(0,1))-2*u+np.roll(u,(0,+1),axis=(0,1)))/hy**2 )
    
    ut = u + dt*(-alpha+delta)
    
    # v_velocity
    q_x = 0.5*hx*( u+np.roll(u,(+1,0),axis=(0,1))+np.roll(u,(0,-1),axis=(0,1))+np.roll(u,(+1,-1),axis=(0,1)) )/4/nu
    q_y = 0.5*hy*v/nu
    
    noTaylor_x = np.where(abs(q_x)>0.1)
    noTaylor_y = np.where(abs(q_y)>0.1)
    
    alpha_x = q_x/3 - q_x**3/45 + 2*q_x**5/945 - q_x**7/4725
    alpha_x[noTaylor_x] = 1/np.tanh(q_x[noTaylor_x]) - 1/q_x[noTaylor_x]
    
    alpha_y = q_y/3 - q_y**3/45 + 2*q_y**5/945 - q_y**7/4725
    alpha_y[noTaylor_y] = 1/np.tanh(q_y[noTaylor_y])-1/q_y[noTaylor_y]
    
    alpha = 0.25*( ( np.roll(v,(0,-1),axis=(0,1))+v )*( (1-alpha_y)*np.roll(v,(0,-1),axis=(0,1))+(1+alpha_y)*v ) - \
                   ( np.roll(v,(0,+1),axis=(0,1))+v )*( (1-alpha_y)*v + (1+alpha_y)*np.roll(v,(0,+1),axis=(0,1)) ) )/hx + \
            0.25*( ( np.roll(u,(0,-1),axis=(0,1))+u )*( (1-alpha_x)*np.roll(v,(-1,0),axis=(0,1))+(1+alpha_x)*v ) - \
                   ( np.roll(u,(+1,-1),axis=(0,1))+np.roll(u,(+1,0),axis=(0,1)) )*( (1-alpha_x)*v + (1+alpha_x)*np.roll(v,(+1,0),axis=(0,1)) ) )/hy 
                  
    delta = nu*( (np.roll(v,(-1,0),axis=(0,1))-2*v+np.roll(v,(+1,0),axis=(0,1)))/hx**2 + (np.roll(v,(0,-1),axis=(0,1))-2*v+np.roll(v,(0,+1),axis=(0,1)))/hy**2 )
    
    vt = v + dt*(-alpha+delta)
    
    ut = solve_u*np.delete(ut,-1,0)
    vt = solve_v*np.delete(vt,-1,1)
    
    return ut,vt

def temporaryVelocityBoundaryConditions(ut,vt,u,v,nx,ny):
    #inlet ut-velocity
    ut[0,1:int(ny/10+1)] = u[0,1:int(ny/10+1)]                                      #13
    ut[0,int(ny/4+1):int(ny*3/4+1)] = u[0,int(ny/4+1):int(ny*3/4+1)]                #14
    ut[0,int(ny*9/10+1):ny+1] = u[0,int(ny*9/10+1):ny+1]                            #15    

    #outlet ut-velocity
    ut[nx,int(ny/4+1):int(ny*3/4+1)] = u[nx,int(ny/4+1):int(ny*3/4+1)]              #16   

    #outlet v-velocity
    vt[nx+1,int(ny/4):int(ny*3/4+1)] = v[nx+1,int(ny/4):int(ny*3/4+1)]              #16  
    return ut,vt

def poissonEquation(p,ut,vt,nx,ny,hx,hy,dt,beta,maxIter,maxError,den,solve,solver_type):
    #SOR algorithm
    if solver_type == 1:
        #ut = np.append(ut,np.zeros((1,np.shape(ut)[1])),axis=0)
        #vt = np.append(vt,np.zeros((np.shape(vt)[0],1)),axis=1)
        
        poissonIteration = 0
        for k in range(0,maxIter):
            p_old = p.copy();
            poissonIteration = poissonIteration + 1
    
            #SOR algorithm
            for i in range(0,nx+2):
                for j in range(0,ny+2):
                    if solve[i,j] == 1:
                        first = hx**2/den[i,j]*(p[i,j+1]+p[i,j-1])
                        second = hy**2/den[i,j]*(p[i+1,j]+p[i-1,j])
                        S = hx*hy/den[i,j]/dt*(hy*(ut[i,j]-ut[i-1,j])+hx*(vt[i,j]-vt[i,j-1]))
                        gaussSeidler = first+second-S
                        p[i,j] = beta*gaussSeidler+(1-beta)*p[i,j]
    
            #residual estimation
            R = 0
            for i in range(0,nx+2):
                for j in range(0,ny+2):
                    if solve[i,j] == 1:
                        R = R + abs(p_old[i,j]-p[i,j]);
    
            R = R/(0.64*nx*ny)
            #check for convergence
            if R<=maxError:
                break       
    else:
    
        ne = (nx+2)*(ny+2)
        b = np.zeros(ne)
        I = np.zeros(ne)
        J = np.zeros(ne)
        V = np.zeros(ne)

        # Fill main diagonal
        counter = 0;
        for i in range(0,ne):
            I[counter] = i
            J[counter] = i
            V[counter] = 1
            counter = counter+1

        # Fill equations
        for i in range(1,nx+1):
            for j in range(1,ny+1):
                if solve[i,j] == 1:
                    k = (nx+2)*j + i;
    
                    #position and value for parameter of p(i+1,j)
                    I = np.append(I,k)
                    J = np.append(J,k+1)
                    V = np.append(V,-hy**2/den[i,j])
                    
                    #position and value for parameter of p(i-1,j)
                    I = np.append(I,k)
                    J = np.append(J,k-1)
                    V = np.append(V,-hy**2/den[i,j])
                    
                    #position and value for parameter of p(i,j+1)
                    I = np.append(I,k)
                    J = np.append(J,k+(nx+2))
                    V = np.append(V,-hx**2/den[i,j])
                    
                    #position and value for parameter of p(i,j-1)
                    I = np.append(I,k)
                    J = np.append(J,k-(nx+2))
                    V = np.append(V,-hx**2/den[i,j])
                    
                    b[k] = -hx*hy/den[i,j]/dt*(hy*(ut[i,j]-ut[i-1,j])+hx*(vt[i,j]-vt[i,j-1]))
                    
        #M = sparse(I,J,V,ne,ne)
        
        M = csc_matrix((V,(I,J)),shape = (ne,ne))
        
        if solver_type == 3:
            #Direct algorithm
            #p = scipy.sparse.linalg.spsolve(M,b);
            p = pypardiso.spsolve(M,b);
            poissonIteration = 1;

        #elif solver_type == 2:
            #GMRES algorithm
            #tol = 1e-5;
            #maxit = 10;
            #[p,~,~,iter] = gmres(M,b,[],tol,maxit,[],[],p(:));
            #poissonIteration=iter(2) 

        p = p.reshape(ny+2,nx+2);
        p = p.transpose();
        
    return p,poissonIteration

def velocityCorrection(u,v,ut,vt,p,nx,ny,hx,hy,dt,solve_u,solve_v):
    
    u = np.append(u,np.zeros((1,np.shape(u)[1])),axis=0)
    v = np.append(v,np.zeros((np.shape(v)[0],1)),axis=1)
    ut = np.append(ut,np.zeros((1,np.shape(u)[1])),axis=0)
    vt = np.append(vt,np.zeros((np.shape(v)[0],1)),axis=1)
    
    # u-velocity
    u = ( ut - dt/hx*(np.roll(p,(-1,0),axis=(0,1))-p) )
    
    # v-velocity
    v = ( vt - dt/hy*(np.roll(p,(0,-1),axis=(0,1))-p) )
    
    ut = np.delete(ut,-1,0)
    vt = np.delete(vt,-1,1)
    u = solve_u*np.delete(u,-1,0)
    v = solve_v*np.delete(v,-1,1)

    return u,v

def reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r):
    '''
    uu = 0.50*( u  + np.roll(u ,(0,-1),axis=(0,1)) )
    vv = 0.50*( v  + np.roll(v ,(-1,0),axis=(0,1)) )
    pp = 0.25*( p  + np.roll(p ,(-1,0),axis=(0,1)) + np.roll(p ,(0,-1),axis=(0,1)) + np.roll(p ,(-1,-1),axis=(0,1)) )
    ca = 0.25*( ca + np.roll(ca,(-1,0),axis=(0,1)) + np.roll(ca,(0,-1),axis=(0,1)) + np.roll(ca,(-1,-1),axis=(0,1)) )
    cb = 0.25*( cb + np.roll(cb,(-1,0),axis=(0,1)) + np.roll(cb,(0,-1),axis=(0,1)) + np.roll(cb,(-1,-1),axis=(0,1)) )
    cc = 0.25*( cc + np.roll(cc,(-1,0),axis=(0,1)) + np.roll(cc,(0,-1),axis=(0,1)) + np.roll(cc,(-1,-1),axis=(0,1)) )
    cd = 0.25*( cd + np.roll(cd,(-1,0),axis=(0,1)) + np.roll(cd,(0,-1),axis=(0,1)) + np.roll(cd,(-1,-1),axis=(0,1)) )
    ce = 0.25*( ce + np.roll(ce,(-1,0),axis=(0,1)) + np.roll(ce,(0,-1),axis=(0,1)) + np.roll(ce,(-1,-1),axis=(0,1)) )
    
    '''
    
    for i in range(0,nx+1):
        for j in range(0,ny+1):
            if solve_r[i,j]==1:
                uu[i,j] = 0.5*(u[i,j]+u[i,j+1]);
                vv[i,j] = 0.5*(v[i,j]+v[i+1,j]);
                pp[i,j] = 0.25*(p[i,j]+p[i+1,j]+p[i,j+1]+p[i+1,j+1]);
                ca[i,j] = 0.25*(CA[i,j]+CA[i+1,j]+CA[i,j+1]+CA[i+1,j+1]);
                cb[i,j] = 0.25*(CB[i,j]+CB[i+1,j]+CB[i,j+1]+CB[i+1,j+1]);
                cc[i,j] = 0.25*(CC[i,j]+CC[i+1,j]+CC[i,j+1]+CC[i+1,j+1]);
                cd[i,j] = 0.25*(CD[i,j]+CD[i+1,j]+CD[i,j+1]+CD[i+1,j+1]);
                ce[i,j] = 0.25*(CE[i,j]+CE[i+1,j]+CE[i,j+1]+CE[i+1,j+1]);

    return uu,vv,pp,ca,cb,cc,cd,ce

def graphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color,typeOfPlot):

    #reconstruct main field
    [uu,vv,pp,ca,cb,cc,cd,ce] = reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r);
    
    if "quiver" in typeOfPlot:
        #velocity field vector (quiver)
        plt.figure(1,figsize=(30,20),dpi=500);
        plt.quiver(X,Y,np.transpose(uu),np.transpose(vv),color = 'blue',width=0.0005);
        plt.title("Velocity Field",fontsize=15); 
        plt.ylabel('y [m]',fontsize=10); plt.xlabel('x [m]',fontsize=10);
        plt.gca().add_patch(Rectangle((position[0,0],position[0,1]),position[0,2],position[0,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[1,0],position[1,1]),position[1,2],position[1,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[2,0],position[2,1]),position[2,2],position[2,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[3,0],position[3,1]),position[3,2],position[3,3],facecolor = color));        
        plt.savefig("quiver.png",dpi = 1000)
        plt.show();
    
    if "streamline" in typeOfPlot:
        #velocity field vector (streamlines)
        plt.figure(2,figsize=(30,20),dpi=500);
        plt.streamplot(X,Y,np.transpose(uu),np.transpose(vv),density=2,arrowsize=1.5);
        plt.title("Streamlines",fontsize=15);
        plt.gca().add_patch(Rectangle((position[0,0],position[0,1]),position[0,2],position[0,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[1,0],position[1,1]),position[1,2],position[1,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[2,0],position[2,1]),position[2,2],position[2,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[3,0],position[3,1]),position[3,2],position[3,3],facecolor = color));
        plt.ylabel('y [m]',fontsize=10); plt.xlabel('x [m]',fontsize=10);
        plt.savefig("streamlines.png",dpi = 1000)
        plt.show()
    
    if "surface" in typeOfPlot:
        #velocity field vector (streamlines)
        plt.figure(3,figsize=(30,20),dpi=500);
        colorVector = ["#FFFFFF","#FFF1F1","#FFD6D6","#FFC5C5","#FFABAB","#FF8F8F","#FF7474","#FF4B4B","#FF0000"]
        plt.contourf(X, Y, np.transpose(uu),colors=colorVector)
        plt.colorbar()
        plt.title("u velocity surface",fontsize=15);
        plt.gca().add_patch(Rectangle((position[0,0],position[0,1]),position[0,2],position[0,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[1,0],position[1,1]),position[1,2],position[1,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[2,0],position[2,1]),position[2,2],position[2,3],facecolor = color));
        plt.gca().add_patch(Rectangle((position[3,0],position[3,1]),position[3,2],position[3,3],facecolor = color));
        plt.ylabel('y [m]',fontsize=10); plt.xlabel('x [m]',fontsize=10);
        plt.savefig("surface.png",dpi = 1000)
        plt.show()

#------------------------------------------------------------------------------------------------------------
#                                          --- IMPORT LIBRARIES ---
#------------------------------------------------------------------------------------------------------------

import numpy as np
import pypardiso
import matplotlib.pyplot as plt
import math
import time
from matplotlib.patches import Rectangle
from scipy.sparse import csc_matrix

#------------------------------------------------------------------------------------------------------------
#                                         --- SETTING UP SIMULATION ---
#------------------------------------------------------------------------------------------------------------

#grid
nx = 100                        #[-]
ny = 80                         #[-]
Lx = 1                          #[m]
Ly = 0.1                        #[m]
hx = Lx/nx                      #[m]
hy = Ly/ny                      #[m]

#grid construction
x = np.linspace(0,Lx,nx+1)      # grid coordinates (x axis)
y = np.linspace(0,Ly,ny+1)      # grid coordinates (y axis)

[X,Y] = np.meshgrid(x,y);       # PYTHON grid

#onTheFly graphical post-processing
onTheFlyPP = False;
graph_step = 100;

#Type of plot available: "streamline" "quiver"
typeOfPlot = ["streamline","quiver","surface"]

#save the results in txt file
saveResults = False

#data
u_top = 1;                             #[m/s]
u_mid = 0.5;                           #[m/s]
u_bot = 1;                             #[m/s]
u_out = (0.2*(u_top+u_bot)+u_mid);     #[m/s]
MW = 0.046                             #[kg/mol]
mu = 1e-4                              #[Pa*s]
pressure = 101325                      #[Pa]
R = 8.31446261815324                   #[J/mol/K]
T = 273.15+70                          #[K]
rho = MW*pressure/R/T                  #[kg/m^3]
nu = mu/rho                            #[m^2/s]
p_star = 0
CA_in = pressure/R/T                   #[mol/^3]
CD_in = pressure/R/T                   #[mol/m^3] 

gamma = 4e-5;
kappa = np.array([3.5,1.2,95e-6]);

#time step
sigma = 0.5
dt_diff=min(hx,hy)**2/4/nu                               # time step (diffusion stability) [s]
dt_conv=4*nu/(max(u_top,max(u_mid,u_bot)))**2            # time step (convection stability) [s]
dt_diff_sp=min(hx,hy)**2/4/gamma                         # time step (species diffusion stability) [s]
dt_conv_sp=4*gamma/(max(u_top,max(u_mid,u_bot)))**2      # time step (species convection stability) [s]

dt_sp = min(dt_conv_sp,dt_diff_sp)                       # time step NAVIER-STOKES equation [s]
dt_ns = min(dt_conv,dt_diff)                             # time step SPECIES EQUATIONS [s]
dt = sigma*min(dt_ns,dt_sp)                              # actual time-step [s]

peclet = max(u_top,max(u_mid,u_bot))*max(hx,hy)/gamma;
print("Peclet number: ",round(peclet))

#simulation time
tau = 10;                                                #[s]
t_steps = tau/dt;

#solver_type for POISSON EQUATION
# 1 = SOR
# 2 = GMRES
# 3 = direct
solver_type = 3;
poissonSolverName = ["SOR","GMRES","DIRECT"]
solver_type_change = 2e3;  #switch poisson solver after "solver_type_change" time iters
solver_change_to = 1;      #solver switch to "solver_change_to"

#SOR 
beta = 1.3;
maxIter = int(1e6);
maxError = 1e-5;

#pre allocation
p,u,v,ut,vt,pp,uu,vv,CA,CB,CC,CD,CE,CAstar,CBstar,CCstar,CDstar,CEstar,ca,cb,cc,cd,ce,den,position = preAllocation(nx,ny);

#pressure inlet/outlet
for i in range(ny+2):
    p[0,i] = p_star
    p[nx+1,i] = p_star

#initialize denominator for poissonEquation
den = initializeDenominator(nx,ny,hx,hy);

#initialize solving domain
[solve,solve_u,solve_v,solve_r] = solvingDomain(nx,ny);

#rectangles position (not solved domain)
position[0,:] = [0,Ly/10,Lx/5,3*Ly/20];
position[1,:] = [0,Ly*3/4,Lx/5,3*Ly/20];
position[2,:] = [Lx*2/5,0,Lx*3/5,Ly/4];
position[3,:] = [Lx*2/5,Ly*3/4,Lx*3/5,Ly/4];
color = '#2D3033';

t = 0;
startTime = time.time()
oldTime = startTime

#------------------------------------------------------------------------------------------------------------
#                                           --- SIMULATION ---
#------------------------------------------------------------------------------------------------------------


for n in range(0,int(round(t_steps+1))):
    
    if n==solver_type_change:
        solver_type = solver_change_to

    #update boundary conditiods
    u,v = velocityBoundaryConditions(u,v,nx,ny,u_top,u_mid,u_bot)
    
    #solve for temporary velocity
    ut,vt = temporaryVelocity(ut,vt,u,v,nx,ny,hx,hy,dt,nu,solve_u,solve_v)
    
    #temporary velocity boundary conditions
    ut,vt = temporaryVelocityBoundaryConditions(ut,vt,u,v,nx,ny)
    
    #solve poisson equation
    p,iter = poissonEquation(p,ut,vt,nx,ny,hx,hy,dt,beta,maxIter,maxError,den,solve,solver_type)
    
    #velocity field correction
    u,v = velocityCorrection(u,v,ut,vt,p,nx,ny,hx,hy,dt,solve_u,solve_v)
    
    #forcing mass conservation (if not guaranteed by the soultion)
    mass_cons = "nop";
    avg_u_velocity_outlet = np.mean((u[nx,int(ny/4+1):int(ny*3/4+2)]+u[nx,int(ny/4):int(ny*3/4+1)])/2)
    if abs(avg_u_velocity_outlet-u_out)>1e-10:
        mass_cons = "yep"
        if abs(avg_u_velocity_outlet)<u_out*0.7:
            u[nx,int(ny/4+1):int(ny*3/4+1)] = u_out
        else:     
            u_correction = u_out/avg_u_velocity_outlet
            u[nx,int(ny/4+1):int(ny*3/4+1)] = u_correction*u[nx,int(ny/4+1):int(ny*3/4+1)]

#------------------------------------------------------------------------------------------------------------
#                                          --- REACTION PART ---
#------------------------------------------------------------------------------------------------------------
    
    #MISSING ...
    
#------------------------------------------------------------------------------------------------------------
#                                          --- FINAL OPERATION ---
#------------------------------------------------------------------------------------------------------------
    
    #advance in time
    t = t+dt;

    if time.time()-oldTime>2:
        if round(time.time()-startTime)<60:
            print("elapsedTime =\t",round(time.time()-startTime)," [s] - (#",n+1, \
                  ") - PS =\t",iter," - t =\t",np.round_(t,3)," [s] - poissonSolver: ",poissonSolverName[solver_type-1])
            oldTime = time.time()
        elif round(time.time()-startTime)>=60 and round(time.time()-startTime)<3600:
            timeMinutes = math.floor((time.time()-startTime)/60)
            timeSeconds = round(time.time()-startTime-60*timeMinutes)
            if timeSeconds>=60:
                timeMinutes = timeMinutes + 1
                timeSeconds = timeSeconds - 60
            print("elapsedTime =\t",timeMinutes,"[min] ",timeSeconds," [s] - (#",n+1, \
                  ") - PS =\t",iter," - t =\t",np.round_(t,3)," [s] - poissonSolver: ",poissonSolverName[solver_type-1])
            oldTime = time.time()
        else:
            if time.time()-oldTime>61:
                timeHours = math.floor((time.time()-startTime)/3600)
                timeMinutes = math.floor((time.time()-startTime-3600*timeHours)/60)
                timeSeconds = round(time.time()-startTime-60*timeMinutes-300*timeHours)
                timeMinutes = timeMinutes + 1
                print("elapsedTime =\t",timeHours," [h] ",timeMinutes," [min] - (#",n+1, \
                      ") - PS =\t",iter," - t =\t",np.round_(t,3)," [s] - poissonSolver: ",poissonSolverName[solver_type-1])
                oldTime = time.time()

    #OnTheFly graphical post-processing
    
    if onTheFlyPP:
        if n%graph_step==0 and n!=0:
            graphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color,typeOfPlot)

#final graphical post-processing
graphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color,typeOfPlot)
            
#SAVE THE RESULTS
if saveResults:
    [uu,vv,pp,ca,cb,cc,cd,ce] = reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r);

    file = open("results.txt","w+")
    uu_str = repr(uu)
    vv_str = repr(vv)
    pp_str = repr(pp)
      
    file.write("uu = " + uu_str + "\n")    
    file.write("vv = " + vv_str + "\n")    
    file.write("pp = " + pp_str + "\n")  

    file.close()      
    
  
    


