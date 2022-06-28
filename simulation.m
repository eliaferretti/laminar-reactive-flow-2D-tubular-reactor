clear variables
close all
clc

tic
cmpt_time0 = cputime;

%grid
nx = 25;
ny = 20;
Lx = 1;                         %[m]
Ly = 0.1;                       %[m]
hx = Lx/nx;                     %[m]
hy = Ly/ny;                     %[m]

%grid construction
x=0:hx:Lx;                      % grid coordinates (x axis)
y=0:hy:Ly;                      % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);          % MATLAB grid

%onTheFly graphical post-processing
onTheFlyPP="false";
graph_step = 5e3;

%data
u_top = 1;                             %[m/s]
u_mid = 0.5;                           %[m/s]
u_bot = 1;                             %[m/s]
u_out = (0.2*(u_top+u_bot)+u_mid);     %[m/s]
MW = 0.046;                            %[kg/mol]
mu = 1e-4;                             %[Pa*s]
pressure = 101325;                     %[Pa]
R = 8.31446261815324;                  %[J/mol/K]
T = 273.15+70;                         %[K]
rho = MW*pressure/R/T;                 %[kg/m^3]
nu = mu/rho;                           %[m^2/s]
p_star = 0;
CA_in = pressure/R/T;                  %[mol/m^3]
CD_in = pressure/R/T;                  %[mol/m^3] 

gamma = 4e-5;
kappa = [3.5 1.2 95e-6];

%time step
sigma = 0.5;
dt_diff=min(hx,hy)^2/4/nu;                               % time step (diffusion stability) [s]
dt_conv=4*nu/(max(u_top,max(u_mid,u_bot)))^2;            % time step (convection stability) [s]
dt_diff_sp=min(hx,hy)^2/4/gamma;                         % time step (species diffusion stability) [s]
dt_conv_sp=4*gamma/(max(u_top,max(u_mid,u_bot)))^2;      % time step (species convection stability) [s]

dt_sp = min(dt_conv_sp,dt_diff_sp);                      % time step NAVIER-STOKES equation [s]
dt_ns = min(dt_conv,dt_diff);                            % time step SPECIES EQUATIONS [s]
dt = sigma*min(dt_ns,dt_sp);                             % actual time-step [s]

peclet = max(u_top,max(u_mid,u_bot))*max(hx,hy)/gamma;
fprintf("Peclet number: %f \n",peclet)

%simulation time
tau = 10;                                                %[s]
t_steps = tau/dt;

%solver_type for POISSON EQUATION
% 1 = SOR
% 2 = GMRES
% 3 = direct
solver_type = 3;
solver_type_change = 1e4;  %switch poisson solver after "solver_type_change" time iters
solver_change_to = 1;      %solver switch to "solver_change_to"

%SOR 
beta = 1.5;
maxIter = 1e6;
maxError = 1e-5;

%pre allocation
[p,u,v,ut,vt,pp,uu,vv,CA,CB,CC,CD,CE,CAstar,CBstar,CCstar,CDstar,CEstar,ca,cb,cc,cd,ce] = preAllocation(nx,ny);

%pressure inlet/outlet
p(1,:) = p_star;
p(nx+2,:) = p_star;

%initialize denominator for poissonEquation
den = zeros(nx+2,ny+2);
[den] = initializeDenominator(den,nx,ny,hx,hy);

%initialize solving domain
[solve,solve_u,solve_v,solve_r] = solvingDomain(nx,ny);

% Rectangle position (not solved domain)
position(1,:) = [0,Ly/10,Lx/5,Ly*3/20];
position(2,:) = [0,Ly*3/4,Lx/5,Ly*3/20];
position(3,:) = [Lx*2/5,0,Lx*3/5,Ly/4];
position(4,:) = [Lx*2/5,Ly*3/4,Lx*3/5,Ly/4];
color = [0.85 0.85 0.85];

t = 0;

%Simulation
for n=1:t_steps+1
    if n>solver_type_change
        solver_type = solver_change_to;
    end
    %update boundary conditiods
    [u,v] = velocityBoundaryConditions(u,v,nx,ny,u_top,u_mid,u_bot);
  
    %solve for temporary velocity
    [ut,vt] = temporaryVelocity(ut,vt,u,v,nx,ny,hx,hy,dt,nu,solve_u,solve_v);
    
    %temporary velocity buondary conditions
    [ut,vt] = temporaryVelocityBoundaryConditions(ut,vt,u,v,nx,ny);

    %solve poisson equation
    [p,iter] = poissonEquation(p,ut,vt,nx,ny,hx,hy,dt,beta,maxIter,maxError,den,solve,solver_type);
    
    %velocity field correction
    [u,v] = velocityCorrection(u,v,ut,vt,p,nx,ny,hx,hy,dt,solve_u,solve_v);
    
    %forcing mass conservation (if not guaranteed by the soultion)
    mass_cons = "nop";
    avg_u_velocity_outlet = mean((u(nx+1,ny/4+2:ny*3/4+2)+u(nx+1,ny/4+1:ny*3/4+1))/2);
    if ( abs(avg_u_velocity_outlet-u_out)>1e-10 )
        mass_cons = "yep";
        if abs(avg_u_velocity_outlet)<u_out*0.7
            u(nx+1,ny/4+2:ny*3/4+1) = u_out;
        else     
            u_correction = u_out/avg_u_velocity_outlet;
            u(nx+1,ny/4+2:ny*3/4+1) = u_correction*u(nx+1,ny/4+2:ny*3/4+1);
        end
    end
    nTS = 2;
    for TS=1:nTS
        dt_TS = dt/nTS;
        %concentration boundary conditions
        [CA,CB,CC,CD,CE] = concentrationBuondaryConditions(CA,CB,CC,CD,CE,nx,ny,CA_in,CD_in);
    
        %species OPERATOR-SPLITTING
    
        [CAstar] = advectionDiffusionEquation(CA,nx,ny,hx,hy,gamma,u,v,dt_TS,solve);
        [CBstar] = advectionDiffusionEquation(CB,nx,ny,hx,hy,gamma,u,v,dt_TS,solve);
        [CCstar] = advectionDiffusionEquation(CC,nx,ny,hx,hy,gamma,u,v,dt_TS,solve);
        [CDstar] = advectionDiffusionEquation(CD,nx,ny,hx,hy,gamma,u,v,dt_TS,solve);
        [CEstar] = advectionDiffusionEquation(CE,nx,ny,hx,hy,gamma,u,v,dt_TS,solve);
        
        for i=1:nx+2
            for j=1:ny+2
                if solve(i,j) == 1               
                    %LOCAL TIME STEP
                    eigen(1) = kappa(2)+kappa(3)*CDstar(i,j)+kappa(3)*CBstar(i,j);
                    eigen(2) = 4*kappa(2)*kappa(3)*CBstar(i,j);
                    dt_chem(1) = abs( ( -kappa(1) )^(-1) );
                    dt_chem(2) = abs( ( 0.5*(-(eigen(1)^2-eigen(2))^0.5 - eigen(1)) + 1e-12)^(-1) );
                    dt_chem(3) = abs( ( 0.5*((eigen(1)^2-eigen(2))^0.5 - eigen(1))  + 1e-12)^(-1) );
                    dt_chem = min(dt_chem);
        
                    n_chemStep = ceil(dt_TS/dt_chem);
                    dt_chem = dt_TS/n_chemStep;
        
                    for k = 1:n_chemStep
                        r(1) = kappa(1)*CAstar(i,j);
                        r(2) = kappa(2)*CBstar(i,j);
                        r(3) = kappa(3)*CBstar(i,j)*CDstar(i,j);
    
                        CAstar(i,j) = CAstar(i,j)-r(1)*dt_chem;
                        CBstar(i,j) = CBstar(i,j)+(r(1)-r(2)-r(3))*dt_chem;
                        CCstar(i,j) = CCstar(i,j)+r(2)*dt_chem;
                        CDstar(i,j) = CDstar(i,j)-r(3)*dt_chem;
                        CEstar(i,j) = CEstar(i,j)+2*r(3)*dt_chem;
                    end
                    
                    CA(i,j) = CAstar(i,j);
                    CB(i,j) = CBstar(i,j);
                    CC(i,j) = CCstar(i,j);
                    CD(i,j) = CDstar(i,j);
                    CE(i,j) = CEstar(i,j);
                end
            end
        end
    end
  
    %advance in time
    t = t+dt;

    %reconstruct main field
    [uu,vv,pp,ca,cb,cc,cd,ce] = reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r);
    
    %molar flux of species E
    n_e = mean( uu(nx+1,ny/4+1:ny*3/4+1).*ce(nx+1,ny/4+1:ny*3/4+1) );
    fprintf("Ne [mol/m^2/s]: %f - PS: %u time [s]: %f (#%u) TOT_time [s]: %u simulationTime [s]: %f \n",n_e,iter,t,n,tau,cputime-cmpt_time0)

    %OnTheFly graphical post-processing
    if onTheFlyPP=="true"
        if(mod(n,graph_step)==0)
            graphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color);
            pause(0.001)
        end
    end
end
close all

%reconstruct main field
[uu,vv,pp,ca,cb,cc,cd,ce] = reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r);

%molar flux of species E
n_e = mean( uu(nx+1,ny/4+1:ny*3/4+1) )*mean( ce(nx+1,ny/4+1:ny*3/4+1) );
fprintf("Molar flux of species E @ the outlet [mol/m^2/s]: %f \n",n_e);

%post-processing
graphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color)
finalGraphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color)

%Total time evaluation
fprintf("Total simulation time [s]: %f \n",cputime-cmpt_time0);
toc

function [solve,solve_u,solve_v,solve_r] = solvingDomain(nx,ny)
    solve = zeros(nx+2,ny+2);
    solve_u = zeros(nx+1,ny+2);
    solve_v = zeros(nx+2,ny+1);
    solve_r = zeros(nx+1,ny+1);
    %1-solve_u staggered grid
    for i=2:nx/5+1
        for j=2:ny/10+1
            solve_u(i,j) = 1;
        end
    end
    %2-solve_u staggered grid
    for i=2:nx/5+1
        for j=ny*9/10+2:ny+1
            solve_u(i,j) = 1;
        end
    end
    %3-solve_u staggered grid
    for i=2:nx/5+1
        for j=ny/4+2:ny*3/4+1
            solve_u(i,j) = 1;
        end
    end
    %4-solve_u staggered grid
    for i=nx/5+2:nx*2/5
        for j=2:ny+1
            solve_u(i,j) = 1;
        end
    end
    %5-solve_u staggered grid
    for i=nx*2/5+1:nx
        for j=ny/4+2:ny*3/4+1
            solve_u(i,j) = 1;
        end
    end

    %1-solve_v staggered grid
    for i=2:nx/5+1
        for j=2:ny/10
            solve_v(i,j) = 1;
        end
    end
    %2-solve_v staggered grid
    for i=2:nx/5+1
        for j=ny*9/10+2:ny
            solve_v(i,j) = 1;
        end
    end
    %3-solve_v staggered grid
    for i=2:nx/5+1
        for j=ny/4+2:ny*3/4
            solve_v(i,j) = 1;
        end
    end
    %4-solve_v staggered grid
    for i=nx/5+2:nx*2/5+1
        for j=2:ny
            solve_v(i,j) = 1;
        end
    end
    %5-solve_v staggered grid
    for i=nx*2/5+2:nx+1
        for j=ny/4+2:ny*3/4
            solve_v(i,j) = 1;
        end
    end
    
    %1-solve
    for i=2:nx/5+1
        for j=2:ny/10+1
            solve(i,j) = 1;
        end
    end
    %2-solve
    for i=2:nx/5+1
        for j=ny*9/10+2:ny+1
            solve(i,j) = 1;
        end
    end
    %3-solve
    for i=2:nx/5+1
        for j=ny/4+2:ny*3/4+1
            solve(i,j) = 1;
        end
    end
    %4-solve
    for i=nx/5+2:nx*2/5+1
        for j=2:ny+1
            solve(i,j) = 1;
        end
    end
    %5-solve
    for i=nx*2/5+2:nx+1
        for j=ny/4+2:ny*3/4+1
            solve(i,j) = 1;
        end
    end

    %1-reconstruct fields
    for i=1:nx/5+1
        for j=1:ny/10+1            
            solve_r(i,j) = 1;
        end
    end
    %2-reconstruct fields
    for i=1:nx/5+1
        for j=ny*9/10+1:ny+1            
            solve_r(i,j) = 1;
        end
    end
    %3-reconstruct fields
    for i=1:nx/5+1
        for j=ny/4+1:ny*3/4+1            
            solve_r(i,j) = 1;
        end
    end
    %4-reconstruct fields
    for i=nx/5+2:nx*2/5+1
        for j=1:ny+1            
            solve_r(i,j) = 1;
        end
    end
    %5-reconstruct fields
    for i=nx*2/5+2:nx+1
        for j=ny/4+1:ny*3/4+1            
            solve_r(i,j) = 1;
        end
    end
end


function [p,u,v,ut,vt,pp,uu,vv,CA,CB,CC,CD,CE,CAstar,CBstar,CCstar,CDstar,CEstar,ca,cb,cc,cd,ce] = preAllocation(nx,ny)
    p = zeros(nx+2,ny+2);
    u = zeros(nx+1,ny+2);
    v = zeros(nx+2,ny+1);
    ut = zeros(nx+1,ny+2);
    vt = zeros(nx+2,ny+1);
    pp = zeros(nx+1,ny+1);
    uu = zeros(nx+1,ny+1);
    vv = zeros(nx+1,ny+1);
    CA = zeros(nx+2,ny+2);
    CB = zeros(nx+2,ny+2);
    CC = zeros(nx+2,ny+2);
    CD = zeros(nx+2,ny+2);
    CE = zeros(nx+2,ny+2);
    CAstar = zeros(nx+2,ny+2);
    CBstar = zeros(nx+2,ny+2);
    CCstar = zeros(nx+2,ny+2);
    CDstar = zeros(nx+2,ny+2);
    CEstar = zeros(nx+2,ny+2);
    ca = zeros(nx+1,ny+1);
    cb = zeros(nx+1,ny+1);
    cc = zeros(nx+1,ny+1);
    cd = zeros(nx+1,ny+1);
    ce = zeros(nx+1,ny+1);
end


function [den] = initializeDenominator(den,nx,ny,hx,hy)
    den(:,:) = 2*(hx^2+hy^2);
    %horizontal sides
    den(2:nx*2/5+1,2) = hx^2+2*hy^2;                %1
    den(2:nx*2/5+1,ny+1) = hx^2+2*hy^2;             %2
    den(nx*2/5+2:nx+1,ny/4+2) = hx^2+2*hy^2;        %3
    den(nx*2/5+2:nx+1,ny*3/4+1) = hx^2+2*hy^2;      %4
    den(2:nx/5+1,ny/10+1) = hx^2+2*hy^2;            %5
    den(2:nx/5+1,ny/4+2) = hx^2+2*hy^2;             %6
    den(2:nx/5+1,ny*3/4+1) = hx^2+2*hy^2;           %7
    den(2:nx/5+1,ny*9/10+2) = hx^2+2*hy^2;          %8
    %vertical sides
    den(nx/5+2,ny/10+2:ny/4+1) = 2*hx^2+hy^2;       %9
    den(nx/5+2,ny*3/4+2:ny*9/10+1) = 2*hx^2+hy^2;   %10
    den(nx*2/5+1,2:ny/4+1) = 2*hx^2+hy^2;           %11
    den(nx*2/5+1,ny*3/4+2:ny+1) = 2*hx^2+hy^2;      %12
    %corner
    den(nx*2/5+1,ny+1) = hx^2+hy^2;
    den(nx*2/5+1,2) = hx^2+hy^2;
end


function [uu,vv,pp,ca,cb,cc,cd,ce] = reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r)
    
    for i=1:nx+1
        for j=1:ny+1
            if ( solve_r(i,j)==1 )
                uu(i,j) = 0.5*(u(i,j)+u(i,j+1));
                vv(i,j) = 0.5*(v(i,j)+v(i+1,j));
                pp(i,j) = 0.25*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1));
                ca(i,j) = 0.25*(CA(i,j)+CA(i+1,j)+CA(i,j+1)+CA(i+1,j+1));
                cb(i,j) = 0.25*(CB(i,j)+CB(i+1,j)+CB(i,j+1)+CB(i+1,j+1));
                cc(i,j) = 0.25*(CC(i,j)+CC(i+1,j)+CC(i,j+1)+CC(i+1,j+1));
                cd(i,j) = 0.25*(CD(i,j)+CD(i+1,j)+CD(i,j+1)+CD(i+1,j+1));
                ce(i,j) = 0.25*(CE(i,j)+CE(i+1,j)+CE(i,j+1)+CE(i+1,j+1));
            end
        end
    end
end


function graphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color)  

    %reconstruct main field
    [uu,vv,pp,ca,cb,cc,cd,ce] = reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r);

    % Surface map: u-velocity
    subplot(251);   
    title('u'); xlabel('x'); ylabel('y');
    surface(X,Y,uu','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    % Surface map: v-velocity
    subplot(252);
    surface(X,Y,vv','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('v'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);    

    % Pressure
    subplot(253);
    surface(X,Y,pp','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('pressure'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);
    
    % Surface map: velocity vectors
    subplot(254);
    quiver(X,Y,uu',vv');
    axis([0 Lx 0 Ly]);
    title('velocity vector field'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    % Streamlines
    subplot(255);
    sx = [Lx/100:-Lx/4000:0 0:Lx/4000:Lx/100];
    sy = [0:Ly/80:Ly/2 Ly/2:Ly/80:Ly];
    streamline(X,Y,uu',vv',sx,sy);
    axis([0 Lx 0 Ly]);
    title('streamlines'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);
    
    %Surface map: species concentration
    subplot(256);
    surface(X,Y,ca','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species A'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    subplot(257);
    surface(X,Y,cb','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species B'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    subplot(258);
    surface(X,Y,cc','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species C'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    subplot(259);
    surface(X,Y,cd','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species D'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    subplot(2,5,10);
    surface(X,Y,ce','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species E'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);  
end


function finalGraphicalPostProcessing(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,Lx,Ly,hx,hy,X,Y,x,y,solve_r,position,color)  

    %reconstruct main field
    [uu,vv,pp,ca,cb,cc,cd,ce] = reconstructFields(u,v,p,CA,CB,CC,CD,CE,uu,vv,pp,ca,cb,cc,cd,ce,nx,ny,solve_r);

    % Surface map: u-velocity
    figure  
    title('u'); xlabel('x'); ylabel('y');
    surface(X,Y,uu','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    % Surface map: v-velocity
    figure
    surface(X,Y,vv','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('v'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);    

    % Pressure
    figure
    surface(X,Y,pp','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('pressure'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);
    
    % Surface map: velocity vectors
    figure
    quiver(X,Y,uu',vv');
    axis([0 Lx 0 Ly]);
    title('velocity vector field'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    % Streamlines
    figure
    sx = [Lx/100:-Lx/16000:0 0:Lx/16000:Lx/100];
    sy = [0:Ly/320:Ly/2 Ly/2:Ly/320:Ly];
    quiver(X,Y,uu',vv');
    streamline(X,Y,uu',vv',sx,sy);
    axis([0 Lx 0 Ly]);
    title('streamlines'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);
    
    %Surface map: species concentration
    figure
    surface(X,Y,ca','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species A'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    figure
    surface(X,Y,cb','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species B'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    figure
    surface(X,Y,cc','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species C'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    figure
    surface(X,Y,cd','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species D'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    figure
    surface(X,Y,ce','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('Concentration of species E'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);

    % total concentration
    ctot = ca + cb + cc + cd + ce;
    figure
    surface(X,Y,ctot','EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colorbar; shading interp;
    title('total concentration'); xlabel('x'); ylabel('y');
    rectangle( 'Position', position(1,:), 'FaceColor', color);
    rectangle( 'Position', position(2,:), 'FaceColor', color);
    rectangle( 'Position', position(3,:), 'FaceColor', color);
    rectangle( 'Position', position(4,:), 'FaceColor', color);
end


function [ut,vt] = temporaryVelocityBoundaryConditions(ut,vt,u,v,nx,ny)
    %inlet ut-velocity
    ut(1,2:ny/10+1) = u(1,2:ny/10+1);                          %13
    ut(1,ny/4+2:ny*3/4+1) = u(1,ny/4+2:ny*3/4+1);              %14
    ut(1,ny*9/10+2:ny+1) = u(1,ny*9/10+2:ny+1);                %15    

    %outlet ut-velocity
    ut(nx+1,ny/4+2:ny*3/4+1) = u(nx+1,ny/4+2:ny*3/4+1);        %16   

    %outlet v-velocity
    vt(nx+2,ny/4+1:ny*3/4+1) = v(nx+2,ny/4+1:ny*3/4+1);         %16  
end


function [u,v] = velocityBoundaryConditions(u,v,nx,ny,u_top,u_mid,u_bot) 
    %horizontal sides u-velcity
    u(1:nx*2/5,1) = -u(1:nx*2/5,2);                             %1
    u(1:nx*2/5,ny+2) = -u(1:nx*2/5,ny+1);                       %2
    u(nx*2/5+1:nx+1,ny/4+1) = -u(nx*2/5+1:nx+1,ny/4+2);         %3
    u(nx*2/5+1:nx+1,ny*3/4+2) = -u(nx*2/5+1:nx+1,ny*3/4+1);     %4
    u(1:nx/5+1,ny/10+2) = -u(1:nx/5+1,ny/10+1);                 %5
    u(1:nx/5+1,ny/4+1) = -u(1:nx/5+1,ny/4+2);                   %6
    u(1:nx/5+1,ny*3/4+2) = -u(1:nx/5+1,ny*3/4+1);               %7
    u(1:nx/5+1,ny*9/10+1) = -u(1:nx/5+1,ny*9/10+2);             %8

    %horizontal sides v-velcity
    v(2:nx*2/5+1,1) = 0;                                        %1
    v(2:nx*2/5+1,ny+1) = 0;                                     %2
    v(nx*2/5+2:nx+1,ny/4+1) = 0;                                %3
    v(nx*2/5+2:nx+1,ny*3/4+1) = 0;                              %4
    v(2:nx/5+1,ny/10+1) = 0;                                    %5
    v(2:nx/5+1,ny/4+1) = 0;                                     %6
    v(2:nx/5+1,ny*3/4+1) = 0;                                   %7
    v(2:nx/5+1,ny*9/10+1) = 0;                                  %8

    %vertical sides u-velocity
    u(nx/5+1,ny/10+2:ny/4+1) = 0;                               %9
    u(nx/5+1,ny*3/4+2:ny*9/10+1) = 0;                           %10
    u(nx*2/5+1,2:ny/4+1) = 0;                                   %11
    u(nx*2/5+1,ny*3/4+2:ny+1) = 0;                              %12

    %vertical sides v-velocity
    v(nx/5+1,ny/10+2:ny/4) = -v(nx/5+2,ny/10+2:ny/4);           %9
    v(nx/5+1,ny*3/4+2:ny*9/10) = -v(nx/5+2,ny*3/4+2:ny*9/10);   %10
    v(nx*2/5+2,1:ny/4) = -v(nx*2/5+1,1:ny/4);                   %11
    v(nx*2/5+2,ny*3/4+2:ny+1) = -v(nx*2/5+1,ny*3/4+2:ny+1);     %12

    %inlet u-velocity
    u(1,2:ny/10+1) = u_bot;                                     %13
    u(1,ny/4+2:ny*3/4+1) = u_mid;                               %14
    u(1,ny*9/10+2:ny+1) = u_top;                                %15    

    %inlet v-velocity
    v(1,1:ny/10+1) = -v(2,1:ny/10+1);                           %13
    v(1,ny/4+1:ny*3/4+1) = -v(2,ny/4+1:ny*3/4+1);               %14
    v(1,ny*9/10+1:ny+1) = -v(1,ny*9/10+1:ny+1);                 %15

    %outlet u-velocity
    u(nx+1,ny/4+2:ny*3/4+1) = u(nx,ny/4+2:ny*3/4+1);            %16   

    %outlet v-velocity
    v(nx+2,ny/4+1:ny*3/4+1) = v(nx+1,ny/4+1:ny*3/4+1);          %16  
end


function [ut,vt] = temporaryVelocity(ut,vt,u,v,nx,ny,hx,hy,dt,nu,solve_u,solve_v)
    % u-velocity
    for i=1:nx+1
        for j=1:ny+2
            if (solve_u(i,j) == 1)
                q_x = 0.5*hx*u(i,j)/nu;
                vy = ( v(i,j)+v(i+1,j)+v(i,j-1)+v(i+1,j-1) )/4;
                q_y = 0.5*hy*vy/nu;
                if abs(q_x)<0.1
                    %taylor expansion
                    alfa_x = q_x/3 - q_x^3/45 + 2*q_x^5/945 - q_x^7/4725;
                else
                    alfa_x = coth(q_x)-1/q_x;
                end
                if abs(q_y)<0.1
                    %taylor expansion
                    alfa_y = q_y/3 - q_y^3/45 + 2*q_y^5/945 - q_y^7/4725;
                else
                    alfa_y = coth(q_y)-1/q_y;
                end
                alfa = 0.25*( ( u(i+1,j)+u(i,j) )*( (1-alfa_x)*u(i+1,j)+(1+alfa_x)*u(i,j) ) - ( u(i,j)+u(i-1,j) )*( (1-alfa_x)*u(i,j)+(1+alfa_x)*u(i-1,j) ) )/hx + ...
                       0.25*( ( v(i+1,j)+v(i,j) )*( (1-alfa_y)*u(i,j+1)+(1+alfa_y)*u(i,j) ) - ( v(i+1,j-1)+v(i,j-1) )*( (1-alfa_y)*u(i,j)+(1+alfa_y)*u(i,j-1) ) )/hy;
                delta = nu*( (u(i+1,j)-2*u(i,j)+u(i-1,j) )/hx^2 + ( u(i,j+1)-2*u(i,j)+u(i,j-1) )/hy^2);
                ut(i,j) = u(i,j) + dt*(-alfa+delta);
            end
        end
    end
    
    % v_velocity
    for i=1:nx+2
        for j=1:ny+1
            if (solve_v(i,j) == 1)
                vx = ( u(i,j)+u(i-1,j)+u(i,j+1)+u(i-1,j+1) )/4;
                q_x = 0.5*hx*vx/nu;
                q_y = 0.5*hy*v(i,j)/nu;
                if abs(q_x)<0.1
                    %taylor expansion
                    alfa_x = q_x/3 - q_x^3/45 + 2*q_x^5/945 - q_x^7/4725;
                else
                    alfa_x = coth(q_x)-1/q_x;
                end
                if abs(q_y)<0.1
                    %taylor expansion
                    alfa_y = q_y/3 - q_y^3/45 + 2*q_y^5/945 - q_y^7/4725;
                else
                    alfa_y = coth(q_y)-1/q_y;
                end
                alfa = 0.25*(( v(i,j+1)+v(i,j) )*( (1-alfa_y)*v(i,j+1)+(1+alfa_y)*v(i,j) ) - ( v(i,j)+v(i,j-1) )*( (1-alfa_y)*v(i,j)+(1+alfa_y)*v(i,j-1) ) )/hy + ...
                       0.25*(( u(i,j+1)+u(i,j) )*( (1-alfa_x)*v(i+1,j)+(1+alfa_x)*v(i,j) ) - ( u(i-1,j+1)+u(i-1,j) )*( (1-alfa_x)*v(i,j)+(1+alfa_x)*v(i-1,j) ))/hx;
                delta = nu*((v(i+1,j)-2*v(i,j)+v(i-1,j))/hx^2+(v(i,j+1)-2*v(i,j)+v(i,j-1))/hy^2);
                vt(i,j) = v(i,j) + dt*(-alfa+delta);
            end
        end
    end
end


function [p,k] = poissonEquation(p,ut,vt,nx,ny,hx,hy,dt,beta,maxIter,maxError,den,solve,solver_type)
    %SOR algorithm
    if solver_type == 1
        for k=1:maxIter
            p_old = p;
            
            %SOR algorithm
            for i=1:nx+2
                for j=1:ny+2
                    if (solve(i,j) == 1)
                        first = hx^2/den(i,j)*(p(i,j+1)+p(i,j-1));
                        second = hy^2/den(i,j)*(p(i+1,j)+p(i-1,j));
                        S = hx*hy/den(i,j)/dt*(hy*(ut(i,j)-ut(i-1,j))+hx*(vt(i,j)-vt(i,j-1)));
                        gaussSeidler = first+second-S;
                        p(i,j) = beta*gaussSeidler+(1-beta)*p(i,j);
                    end
                end
            end
    
            %residuale estimation
            R = 0;
            for i=1:nx+2
                for j=1:ny+2
                    if (solve(i,j) == 1)
                        R = R + abs(p_old(i,j)-p(i,j));
                    end
                end
            end
    
            R = R/(0.64*nx*ny);
            %check for convergence
            if R<=maxError
                break;
            end
        end
    else
    
        ne = (nx+2)*(ny+2);
        b = zeros(ne,1);
        I = zeros(ne,1);
        J = zeros(ne,1);
        V = zeros(ne,1);

        % Fill main diagonal
        counter = 1;
        for i=1:ne
            I(counter) = i; J(counter) = i; V(counter) = 1.; counter = counter+1;
        end

        % Fill equations
        for i=2:nx+1
            for j=2:ny+1
                if (solve(i,j) == 1)
                    k = (nx+2)*(j-1) + i;
    
                    I(counter) = k; J(counter) = k+1; V(counter) = -hy^2/den(i,j); counter = counter+1;
                    I(counter) = k; J(counter) = k-1; V(counter) = -hy^2/den(i,j); counter = counter+1;
                    I(counter) = k; J(counter) = k+(nx+2); V(counter) = -hx^2/den(i,j); counter = counter+1;
                    I(counter) = k; J(counter) = k-(nx+2); V(counter) = -hx^2/den(i,j); counter = counter+1;
    
                    b(k) = -hx*hy/den(i,j)/dt*(hy*(ut(i,j)-ut(i-1,j))+hx*(vt(i,j)-vt(i,j-1)));
                end
            end
        end

        M = sparse(I,J,V,ne,ne);

        if solver_type == 3
            %Direct algorithm
            p = M\b;
            k = 0;

        elseif solver_type == 2
            %GMRES algorithm
            tol = 1e-5;
            maxit = 10;
            [p,~,~,iter] = gmres(M,b,[],tol,maxit,[],[],p(:));
            k=iter(2);
        end
        p = reshape(p,[nx+2 ny+2]);
    end    
end


function [u,v] = velocityCorrection(u,v,ut,vt,p,nx,ny,hx,hy,dt,solve_u,solve_v)
    % u-velocity
    for i=1:nx+1
        for j=1:ny+2
            if (solve_u(i,j) == 1)
                u(i,j) = ut(i,j)-dt/hx*(p(i+1,j)-p(i,j));
            end
        end
    end

    % v-velocity
    for i=1:nx+2
        for j=1:ny+1
            if (solve_v(i,j) == 1)
                v(i,j) = vt(i,j)-dt/hy*(p(i,j+1)-p(i,j));
            end
        end
    end
end


function [CA,CB,CC,CD,CE] = concentrationBuondaryConditions(CA,CB,CC,CD,CE,nx,ny,CA_in,CD_in)

    %inlet
    CA(1,2:ny/10+2) = 2*CA_in - CA(2,2:ny/10+2);
    CA(1,ny*9/10+2:ny+1) = 2*CA_in - CA(2,ny*9/10+2:ny+1);
    CA(1,ny/4+2:ny*3/4+1) = -CA(2,ny/4+2:ny*3/4+1);
    
    CB(1,ny/4+2:ny*3/4+1) = -CB(2,ny/4+2:ny*3/4+1);
    CB(1,2:ny/10+2) = -CB(2,2:ny/10+2);
    CB(1,ny*9/10+2:ny+1) = -CB(2,ny*9/10+2:ny+1);

    CC(1,ny/4+2:ny*3/4+1) = -CC(2,ny/4+2:ny*3/4+1);
    CC(1,2:ny/10+2) = -CC(2,2:ny/10+2);
    CC(1,ny*9/10+2:ny+1) = -CC(2,ny*9/10+2:ny+1);

    CD(1,ny/4+2:ny*3/4+1) = 2*CD_in - CD(2,ny/4+2:ny*3/4+1);
    CD(1,2:ny/10+2) = -CD(2,2:ny/10+2);
    CD(1,ny*9/10+2:ny+1) = -CD(2,ny*9/10+2:ny+1);

    CE(1,ny/4+2:ny*3/4+1) = -CE(2,ny/4+2:ny*3/4+1);
    CE(1,2:ny/10+2) = -CE(2,2:ny/10+2);
    CE(1,ny*9/10+2:ny+1) = -CE(2,ny*9/10+2:ny+1);
    
    %CA-horizontal sides
    CA(2:nx*2/5+1,1) = CA(2:nx*2/5+1,2);                          %1
    CA(2:nx*2/5+1,ny+2) = CA(2:nx*2/5+1,ny+1);                    %2
    CA(nx*2/5+2:nx+1,ny/4+1) = CA(nx*2/5+2:nx+1,ny/4+2);          %3
    CA(nx*2/5+2:nx+1,ny*3/4+2) = CA(nx*2/5+2:nx+1,ny*3/4+1);      %4
    CA(2:nx/5+1,ny/10+2) = CA(2:nx/5+1,ny/10+1);                  %5
    CA(2:nx/5+1,ny/4+1) = CA(2:nx/5+1,ny/4+2);                    %6
    CA(2:nx/5+1,ny*3/4+2) = CA(2:nx/5+1,ny*3/4+1);                %7
    CA(2:nx/5+1,ny*9/10+1) = CA(2:nx/5+1,ny*9/10+2);              %8
    %CA-vertical sides
    CA(nx/5+1,ny/10+2:ny/4+1) = CA(nx/5+2,ny/10+2:ny/4+1);        %9
    CA(nx/5+1,ny*3/4+2:ny*9/10+1) = CA(nx/5+2,ny*3/4+2:ny*9/10+1);%10
    CA(nx*2/5+2,2:ny/4+1) = CA(nx*2/5+1,2:ny/4+1);                %11
    CA(nx*2/5+2,ny*3/4+2:ny+1) = CA(nx*2/5+1,ny*3/4+2:ny+1);      %12
    CA(nx+2,ny/4+1:ny*3/4+2) = CA(nx+1,ny/4+1:ny*3/4+2);          %16

    %CB-horizontal sides
    CB(2:nx*2/5+1,1) = CB(2:nx*2/5+1,2);                          %1
    CB(2:nx*2/5+1,ny+2) = CB(2:nx*2/5+1,ny+1);                    %2
    CB(nx*2/5+2:nx+1,ny/4+1) = CB(nx*2/5+2:nx+1,ny/4+2);          %3
    CB(nx*2/5+2:nx+1,ny*3/4+2) = CB(nx*2/5+2:nx+1,ny*3/4+1);      %4
    CB(2:nx/5+1,ny/10+2) = CB(2:nx/5+1,ny/10+1);                  %5
    CB(2:nx/5+1,ny/4+1) = CB(2:nx/5+1,ny/4+2);                    %6
    CB(2:nx/5+1,ny*3/4+2) = CB(2:nx/5+1,ny*3/4+1);                %7
    CB(2:nx/5+1,ny*9/10+1) = CB(2:nx/5+1,ny*9/10+2);              %8
    %CB-vertical sides
    CB(nx/5+1,ny/10+2:ny/4+1) = CB(nx/5+2,ny/10+2:ny/4+1);        %9
    CB(nx/5+1,ny*3/4+2:ny*9/10+1) = CB(nx/5+2,ny*3/4+2:ny*9/10+1);%10
    CB(nx*2/5+2,2:ny/4+1) = CB(nx*2/5+1,2:ny/4+1);                %11
    CB(nx*2/5+2,ny*3/4+2:ny+1) = CB(nx*2/5+1,ny*3/4+2:ny+1);      %12
    CB(nx+2,ny/4+1:ny*3/4+2) = CB(nx+1,ny/4+1:ny*3/4+2);          %16

    %CC-horizontal sides
    CC(2:nx*2/5+1,1) = CC(2:nx*2/5+1,2);                          %1
    CC(2:nx*2/5+1,ny+2) = CC(2:nx*2/5+1,ny+1);                    %2
    CC(nx*2/5+2:nx+1,ny/4+1) = CC(nx*2/5+2:nx+1,ny/4+2);          %3
    CC(nx*2/5+2:nx+1,ny*3/4+2) = CC(nx*2/5+2:nx+1,ny*3/4+1);      %4
    CC(2:nx/5+1,ny/10+2) = CC(2:nx/5+1,ny/10+1);                  %5
    CC(2:nx/5+1,ny/4+1) = CC(2:nx/5+1,ny/4+2);                    %6
    CC(2:nx/5+1,ny*3/4+2) = CC(2:nx/5+1,ny*3/4+1);                %7
    CC(2:nx/5+1,ny*9/10+1) = CC(2:nx/5+1,ny*9/10+2);              %8
    %CC-vertical sides
    CC(nx/5+1,ny/10+2:ny/4+1) = CC(nx/5+2,ny/10+2:ny/4+1);        %9
    CC(nx/5+1,ny*3/4+2:ny*9/10+1) = CC(nx/5+2,ny*3/4+2:ny*9/10+1);%10
    CC(nx*2/5+2,2:ny/4+1) = CC(nx*2/5+1,2:ny/4+1);                %11
    CC(nx*2/5+2,ny*3/4+2:ny+1) = CC(nx*2/5+1,ny*3/4+2:ny+1);      %12
    CC(nx+2,ny/4+1:ny*3/4+2) = CC(nx+1,ny/4+1:ny*3/4+2);          %16

    %CD-horizontal sides
    CD(2:nx*2/5+1,1) = CD(2:nx*2/5+1,2);                          %1
    CD(2:nx*2/5+1,ny+2) = CD(2:nx*2/5+1,ny+1);                    %2
    CD(nx*2/5+2:nx+1,ny/4+1) = CD(nx*2/5+2:nx+1,ny/4+2);          %3
    CD(nx*2/5+2:nx+1,ny*3/4+2) = CD(nx*2/5+2:nx+1,ny*3/4+1);      %4
    CD(2:nx/5+1,ny/10+2) = CD(2:nx/5+1,ny/10+1);                  %5
    CD(2:nx/5+1,ny/4+1) = CD(2:nx/5+1,ny/4+2);                    %6
    CD(2:nx/5+1,ny*3/4+2) = CD(2:nx/5+1,ny*3/4+1);                %7
    CD(2:nx/5+1,ny*9/10+1) = CD(2:nx/5+1,ny*9/10+2);              %8
    %CD-vertical sides
    CD(nx/5+1,ny/10+2:ny/4+1) = CD(nx/5+2,ny/10+2:ny/4+1);        %9
    CD(nx/5+1,ny*3/4+2:ny*9/10+1) = CD(nx/5+2,ny*3/4+2:ny*9/10+1);%10
    CD(nx*2/5+2,2:ny/4+1) = CD(nx*2/5+1,2:ny/4+1);                %11
    CD(nx*2/5+2,ny*3/4+2:ny+1) = CD(nx*2/5+1,ny*3/4+2:ny+1);      %12
    CD(nx+2,ny/4+1:ny*3/4+2) = CD(nx+1,ny/4+1:ny*3/4+2);          %16

    %CE-horizontal sides
    CE(2:nx*2/5+1,1) = CE(2:nx*2/5+1,2);                          %1
    CE(2:nx*2/5+1,ny+2) = CE(2:nx*2/5+1,ny+1);                    %2
    CE(nx*2/5+2:nx+1,ny/4+1) = CE(nx*2/5+2:nx+1,ny/4+2);          %3
    CE(nx*2/5+2:nx+1,ny*3/4+2) = CE(nx*2/5+2:nx+1,ny*3/4+1);      %4
    CE(2:nx/5+1,ny/10+2) = CE(2:nx/5+1,ny/10+1);                  %5
    CE(2:nx/5+1,ny/4+1) = CE(2:nx/5+1,ny/4+2);                    %6
    CE(2:nx/5+1,ny*3/4+2) = CE(2:nx/5+1,ny*3/4+1);                %7
    CE(2:nx/5+1,ny*9/10+1) = CE(2:nx/5+1,ny*9/10+2);              %8
    %CE-vertical sides
    CE(nx/5+1,ny/10+2:ny/4+1) = CE(nx/5+2,ny/10+2:ny/4+1);        %9
    CE(nx/5+1,ny*3/4+2:ny*9/10+1) = CE(nx/5+2,ny*3/4+2:ny*9/10+1);%10
    CE(nx*2/5+2,2:ny/4+1) = CE(nx*2/5+1,2:ny/4+1);                %11
    CE(nx*2/5+2,ny*3/4+2:ny+1) = CE(nx*2/5+1,ny*3/4+2:ny+1);      %12
    CE(nx+2,ny/4+1:ny*3/4+2) = CE(nx+1,ny/4+1:ny*3/4+2);          %16

    %outlet
    CA(end,ny/4+2:ny*3/4+1) = CA(end-1,ny/4+2:ny*3/4+1);
    CB(end,ny/4+2:ny*3/4+1) = CB(end-1,ny/4+2:ny*3/4+1);
    CC(end,ny/4+2:ny*3/4+1) = CC(end-1,ny/4+2:ny*3/4+1);
    CD(end,ny/4+2:ny*3/4+1) = CD(end-1,ny/4+2:ny*3/4+1);
    CE(end,ny/4+2:ny*3/4+1) = CE(end-1,ny/4+2:ny*3/4+1);
end


function [f] = advectionDiffusionEquation(f,nx,ny,hx,hy,gamma,u,v,dt,solve)
        fo = f;

        for i=1:nx+2
            for j=1:ny+2
                if (solve(i,j) == 1)
                    vx = (u(i,j)+u(i-1,j))/2;
                    vy = (v(i,j)+v(i,j-1))/2;
                    q_x = 0.5*hx*vx/gamma;
                    q_y = 0.5*hy*vy/gamma;
                    if abs(q_x)<0.1
                        %taylor expansion
                        alfa_x = q_x/3 - q_x^3/45 + 2*q_x^5/945 - q_x^7/4725;
                    else
                        alfa_x = coth(q_x)-1/q_x;
                    end
                    if abs(q_y)<0.1
                        %taylor expansion
                        alfa_y = q_y/3 - q_y^3/45 + 2*q_y^5/945 - q_y^7/4725;
                    else
                        alfa_y = coth(q_y)-1/q_y;
                    end
                    alfa = 0.5/hx*( u(i,j)*((1-alfa_x)*fo(i+1,j)+(1+alfa_x)*fo(i,j)) - u(i-1,j)*((1-alfa_x)*fo(i,j)+(1+alfa_x)*fo(i-1,j)) ) + ...
                           0.5/hy*( v(i,j)*((1-alfa_y)*fo(i,j+1)+(1+alfa_y)*fo(i,j)) - v(i,j-1)*((1-alfa_y)*fo(i,j)+(1+alfa_y)*fo(i,j-1)) );
                    beta = gamma*( (fo(i+1,j)-2*fo(i,j)+fo(i-1,j))/hx^2 + (fo(i,j+1)-2*fo(i,j)+fo(i,j-1))/hy^2 );
                    f(i,j) = fo(i,j) + dt*(-alfa+beta);
                end
            end
        end
end


