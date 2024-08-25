% Hua-sheng XIE, FDTD 2D refraction
% Update from Ananth Krishnan's version
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "2D FDTD solution refraction at the interface of two media"
% 
% Objective of the program is to solve for the Maxwell's equation for a TM 
% wave containing the xy-plane polarized magnetic field having components Hy
% and Hx and z-polarized electric field Ez in a domain having two media which 
% are glass in the upper half and air in the lower half. The field update is 
% done using standard update equations obtained from the difference form of 
% Maxwell's curl equations. The field points are defined in a grid described 
% by Yee's algorithm. The H fields are defined at every half coordinate of 
% spacesteps. More precisely, the Hx part is defined at every half y-coordinate 
% and full x-coordinate and the Hy part is defined at every half x-coordinate 
% and full y-coordinate and E fields i.e the Ez part is defined at every full 
% x and full y-coordinate points.Also here, the space-step length is taken 
% as 1 micron. 
%
% Also, the time update is done using Leapfrog time-stepping. Here, H-fields
% i.e. Hx and Hy are updated every half time-step and E fields i.e Ez are 
% updated every full time-step. This is shown by two alternating matrix updates 
% spanning the entire spatial grid at every instant of time. These spatial 
% updates are inside the main for-loop for time update, spanning the entire 
% time grid. Also, here, the matrices used as multiplication factors for update 
% equations are initialized before the loop starts to avoid repeated calculation 
% of the same in every loop iteration, a minor attempt at optimization. One of 
% the boundary condition options here is Mur's Absorbing Boundary Condition (ABC)
% where the fields at the grid points have electric field values formulated 
% using Engquist Majda one way wave equations [1] where the boundaries give 
% a sense of absorbing the total field incident on them and reflecting none 
% back to the domain. Another option here is the Perfectly Matched Layer (PML) 
% boundary condition where the fields near the boundary are attenuated over 
% a predetermined length of boundary width before they reach the boudary to 
% a zero value at the boundary using a polynomially increasing electrical 
% conductivity value over the boundary width with maximum at the boundary and 
% also chosing a magnetic conductivity value at every point in the boundary 
% width to avoid reflection at that point given in [2]. Also, here, Berenger's 
% PML condition is used where in the field Ez is split into two components Ezx 
% and Ezy and the components are attenuated using separate electric and magnetic 
% conductivities in the two directions (sigmax and sigma_starx in x direction 
% and sigmay and sigma_stary in the y direction).
%
% A source of electric field is defined at the top edge of the spatial domain,(in 
% glass) on the left half, which is a hard sinusoidal source, in that it does not 
% change its value due to interference from external fields i.e in other words, 
% the source is a perfect electric conductor. The free space wavelength and the 
% angle of the plane of propagation of the wave from the source with the x-axis 
% can be specified using the variables free_space_wavelength and alpha. The color
% scaled plot of Ez field over the entire spatial domain is shown at every time 
% step as the wave propagates and hits the interface of the two media where it 
% refracts (bends) into the other medium (viz. air) as angle alpha would be lesser 
% than critical angle (~42 degrees) between the two interfaces. The simulation 
% can be ended by closing this plot window or by waiting till all the time step 
% updates are completed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

%Boundary condition selection, any one can be selected by making it 1
pml=1;
abc=0;

%Courant stability factor
S=1/(2^0.5);

% Parameters of free space (permittivity and permeability and speed of
% light) are given corresponding values
epsilon0=(1/(36*pi))*1e-9;
mu0=4*pi*1e-7;
c=3e+8;

% Spatial grid step length (spatial grid step= 1 micron and can be changed)
delta=1e-6;
% Temporal grid step obtained using Courant condition
deltat=S*delta/c;

%Total no of time steps
time_tot=1000;

% Grid Dimension in x (xdim) and y (ydim) directions
ydim=240;%The domain is 240 space steps or 240*1=240 microns long
xdim=400;%The domain is 400 space steps or 400*1=400 microns wide

%Free-space wavelength 
free_space_wavelength=10e-6;

%Index of the top half (glass)
index=1.5;

% Initialization of permittivity and permeability matrices
epsilon=epsilon0*ones(xdim,ydim);
mu=mu0*ones(xdim,ydim);


% Defining of the permittivity profile of the region:-

% Specifying the top half of the spatial domain which is made of glass
% epsilon(:,1:ydim/2)=index*index*epsilon0;
epsilon(:,ydim/2:end)=index*index*epsilon0;

% 2D FDTD update for PML boundary as used in previous program
if pml==1
    
    % Initialization of field matrices
    Ez=zeros(xdim,ydim);
    Ezx=zeros(xdim,ydim);
    Ezy=zeros(xdim,ydim);
    Hy=zeros(xdim,ydim);
    Hx=zeros(xdim,ydim);
    
    % Initializing electric conductivity matrices in x and y directions
    sigmax=zeros(xdim,ydim);
    sigmay=zeros(xdim,ydim);
    
    
    %Perfectly matched layer boundary design
    %[2]-Reference:-http://dougneubauer.com/wp-content/uploads/wdata/yee2dpml1/yee2d_c.txt
    %(An adaptation of 2-D FDTD TE code of Dr. Susan Hagness)
    
    %Boundary width of PML in all directions
    bound_width=20;
    
    %Order of polynomial on which sigma is modeled
    gradingorder=6;
    
    %Required reflection co-efficient
    refl_coeff=1e-6;
    
    %Polynomial model for sigma
    sigmamax=(-log10(refl_coeff)*(gradingorder+1)*epsilon0*c)/(2*bound_width*delta);
    boundfact1=((epsilon(xdim/2,bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact2=((epsilon(xdim/2,ydim-bound_width)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact3=((epsilon(bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    boundfact4=((epsilon(xdim-bound_width,ydim/2)/epsilon0)*sigmamax)/((bound_width^gradingorder)*(gradingorder+1));
    x=0:1:bound_width;
    for i=1:1:xdim
        sigmax(i,bound_width+1:-1:1)=boundfact1*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
        sigmax(i,ydim-bound_width:1:ydim)=boundfact2*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1));
    end
    for i=1:1:ydim
        sigmay(bound_width+1:-1:1,i)=boundfact3*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
        sigmay(xdim-bound_width:1:xdim,i)=boundfact4*((x+0.5*ones(1,bound_width+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width)]).^(gradingorder+1))';
    end
    
    %Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
    %This is also split into x and y directions in Berenger's model
    sigma_starx=(sigmax.*mu)./epsilon;
    sigma_stary=(sigmay.*mu)./epsilon;
    
    %Multiplication factor matrices for H matrix update to avoid being calculated many times
    %in the time update loop so as to increase computation speed
    G=((mu-0.5*deltat*sigma_starx)./(mu+0.5*deltat*sigma_starx));
    H=(deltat/delta)./(mu+0.5*deltat*sigma_starx);
    A=((mu-0.5*deltat*sigma_stary)./(mu+0.5*deltat*sigma_stary));
    B=(deltat/delta)./(mu+0.5*deltat*sigma_stary);
    
    %Multiplication factor matrices for E matrix update to avoid being calculated many times
    %in the time update loop so as to increase computation speed
    C=((epsilon-0.5*deltat*sigmax)./(epsilon+0.5*deltat*sigmax));
    D=(deltat/delta)./(epsilon+0.5*deltat*sigmax);
    E=((epsilon-0.5*deltat*sigmay)./(epsilon+0.5*deltat*sigmay));
    F=(deltat/delta)./(epsilon+0.5*deltat*sigmay);
    
    % Update loop begins
    for n=1:1:time_tot
        
        %matrix update instead of for-loop for Hy and Hx fields
        Hy(1:xdim-1,1:ydim-1)=A(1:xdim-1,1:ydim-1).*Hy(1:xdim-1,1:ydim-1)+B(1:xdim-1,1:ydim-1).*(Ezx(2:xdim,1:ydim-1)-Ezx(1:xdim-1,1:ydim-1)+Ezy(2:xdim,1:ydim-1)-Ezy(1:xdim-1,1:ydim-1));
        Hx(1:xdim-1,1:ydim-1)=G(1:xdim-1,1:ydim-1).*Hx(1:xdim-1,1:ydim-1)-H(1:xdim-1,1:ydim-1).*(Ezx(1:xdim-1,2:ydim)-Ezx(1:xdim-1,1:ydim-1)+Ezy(1:xdim-1,2:ydim)-Ezy(1:xdim-1,1:ydim-1));
        
        %matrix update instead of for-loop for Ez field
        Ezx(2:xdim,2:ydim)=C(2:xdim,2:ydim).*Ezx(2:xdim,2:ydim)+D(2:xdim,2:ydim).*(-Hx(2:xdim,2:ydim)+Hx(2:xdim,1:ydim-1));
        Ezy(2:xdim,2:ydim)=E(2:xdim,2:ydim).*Ezy(2:xdim,2:ydim)+F(2:xdim,2:ydim).*(Hy(2:xdim,2:ydim)-Hy(1:xdim-1,2:ydim));
        
        % Source condition incorporating given free space wavelength 'free_space_wavelength'
        % and having a location at the left half of the top edge (in glass) of the domain just
        % after the PML boundary with the plane of propagation of the emanating wave inclined 
        % at an angle 'alpha' with x-axis
        N_lambda=free_space_wavelength/delta;
        t_start=1;
        alpha=35; %Angle 'alpha' should be < 42 degrees for refraction
        i=75:1:125;
        Ezx(i,bound_width+1)=0.5*sin(((2*pi*c*S*index)/(N_lambda*delta))*(n-t_start)*deltat*(i./i)-((2*pi*(i-76)*index*sin(alpha*pi/180))/(N_lambda)));
        Ezy(i,bound_width+1)=0.5*sin(((2*pi*c*S*index)/(N_lambda*delta))*(n-t_start)*deltat*(i./i)-((2*pi*(i-76)*index*sin(alpha*pi/180))/(N_lambda)));
                
        Ez=Ezx+Ezy;
        
        %Movie type colour scaled image plot of Ez
        h=imagesc(delta*1e+6*(1:1:xdim),(delta*1e+6*(1:1:ydim))',Ez',[-1,1]);colorbar;
        set(h,'AlphaData',10*epsilon'/epsilon0);
        title(['\fontsize{20}Color-scaled image plot of Ez to see refraction at glass/air interface with PML boundary at time = ',num2str(n*deltat*1e+15),' fs']);
        xlabel('x in microns','FontSize',20);
        ylabel('y in microns','FontSize',20);
        set(gca,'FontSize',20);
        getframe;
    end
end

%Fdtd update for Mur's absorbing boundary conditions
if abc==1
    
    % Initialization of field matrices
    Ez=zeros(xdim,ydim);
    Hy=zeros(xdim,ydim);
    Hx=zeros(xdim,ydim);
    
    % Initializing electric and magnetic conductivity matrices
    sigma=4e-4*ones(xdim,ydim);
    sigma_star=4e-4*ones(xdim,ydim);
    
    %Multiplication factor matrices for H matrix update to avoid being calculated many times
    %in the time update loop so as to increase computation speed
    A=((mu-0.5*deltat*sigma_star)./(mu+0.5*deltat*sigma_star));
    B=(deltat/delta)./(mu+0.5*deltat*sigma_star);
    
    %Multiplication factor matrices for E matrix update to avoid being calculated many times
    %in the time update loop so as to increase computation speed
    C=((epsilon-0.5*deltat*sigma)./(epsilon+0.5*deltat*sigma));
    D=(deltat/delta)./(epsilon+0.5*deltat*sigma);
    
    %Mur's absorbing boundary condition parameters
    p0=1;
    p2=-0.5;
    
    %Co-efficients of present and previous (regarding time-step) boundary Ez values
    %in boundary update equation for forward/up and backward/down boundaries in
    %the domain (as given in [1])
    c0=(c/(2*S))*(1-(p0/S));
    c1=-(c/(2*S))*(1+(p0/S));
    c2=(c/(S^2))*(p0+(p2*S*S));
    c3=-(p2*c)/2;
    c0efffor=-(c0/c1);
    c2efffor=-(c2/c1);
    c3efffor=-(c3/c1);
    c0=(c/(2*S))*(1+(p0/S));
    c1=-(c/(2*S))*(1-(p0/S));
    c2=-(c/(S^2))*(p0+(p2*S*S));
    c3=(p2*c)/2;
    c1effrev=-(c1/c0);
    c2effrev=-(c2/c0);
    c3effrev=-(c3/c0);
    
    %Storage vectors for Ez boundary and boundary-1 values of previous and its
    %previous timesteps
    prev_xfor=zeros(1,ydim);
    prev_x_minus_1for=zeros(1,ydim);
    prev_yfor=zeros(xdim,1);
    prev_y_minus_1for=zeros(xdim,1);
    prev_xrev=zeros(1,ydim);
    prev_x_minus_1rev=zeros(1,ydim);
    prev_yrev=zeros(xdim,1);
    prev_y_minus_1rev=zeros(xdim,1);
    
    % Update loop begins
    for n=1:1:time_tot
        
        %Vector update instead of for-loop for Hy and Hx fields
        Hx(2:xdim-3,2:ydim-3)=A(2:xdim-3,2:ydim-3).*Hx(2:xdim-3,2:ydim-3)-B(2:xdim-3,2:ydim-3).*(Ez(2:xdim-3,3:ydim-2)-Ez(2:xdim-3,2:ydim-3));
        Hy(2:xdim-3,2:ydim-3)=A(2:xdim-3,2:ydim-3).*Hy(2:xdim-3,2:ydim-3)+B(2:xdim-3,2:ydim-3).*(Ez(3:xdim-2,2:ydim-3)-Ez(2:xdim-3,2:ydim-3));
        
        %Vector update instead of for-loop for Ez field
        Ez(3:xdim-3,3:ydim-3)=C(3:xdim-3,3:ydim-3).*Ez(3:xdim-3,3:ydim-3)+(Hy(3:xdim-3,3:ydim-3)-Hy(2:xdim-4,3:ydim-3)-Hx(3:xdim-3,3:ydim-3)+Hx(3:xdim-3,2:ydim-4)).*D(3:xdim-3,3:ydim-3);
        
        %Mur's abc conditions obtained from Mur's difference equation for
        %forward boundary
        if n>1
            Ez(xdim-2,3:1:ydim-3)=c0efffor*(Ez(xdim-3,3:1:ydim-3)+prev_prev_xfor(1,3:1:ydim-3))-prev_prev_x_minus_1for(1,3:1:ydim-3)+c2efffor*(prev_xfor(1,3:1:ydim-3)+prev_x_minus_1for(1,3:1:ydim-3))+c3efffor*(prev_x_minus_1for(1,2:1:ydim-4)+prev_x_minus_1for(1,4:1:ydim-2)+prev_xfor(1,2:1:ydim-4)+prev_xfor(1,4:1:ydim-2));
        end
        
        %Storage vectors for boundary and boundary-1 values of previous and its
        %previous time steps updated at forward boundary
        prev_prev_xfor=prev_xfor;
        prev_prev_x_minus_1for=prev_x_minus_1for;
        prev_xfor(1,1:1:ydim)=Ez(xdim-2,1:1:ydim);
        prev_x_minus_1for(1,1:1:ydim)=Ez(xdim-3,1:1:ydim);
        
        %Mur's abc conditions obtained from Mur's difference equation for
        %backward boundary
        if n>1
            Ez(2,3:1:ydim-3)=-prev_prev_xrev(1,3:1:ydim-3)+c1effrev*(Ez(3,3:1:ydim-3)+prev_prev_x_minus_1rev(1,3:1:ydim-3))+c2effrev*(prev_xrev(1,3:1:ydim-3)+prev_x_minus_1rev(1,3:1:ydim-3))+c3effrev*(prev_x_minus_1rev(1,2:1:ydim-4)+prev_x_minus_1rev(1,4:1:ydim-2)+prev_xrev(1,2:1:ydim-4)+prev_xrev(1,4:1:ydim-2));
        end
        
        %Storage vectors for boundary and boundary-1 values of previous and its
        %previous time steps updated at backward boundary
        prev_prev_xrev=prev_xrev;
        prev_prev_x_minus_1rev=prev_x_minus_1rev;
        prev_xrev(1,1:1:ydim)=Ez(3,1:1:ydim);
        prev_x_minus_1rev(1,1:1:ydim)=Ez(2,1:1:ydim);
        
        %Mur's abc conditions obtained from Mur's difference equation for
        %upward boundary
        if n>1
            Ez(3:1:xdim-3,ydim-2)=c0efffor*(Ez(3:1:xdim-3,ydim-3)+prev_prev_yfor(3:1:xdim-3,1))-prev_prev_y_minus_1for(3:1:xdim-3,1)+c2efffor*(prev_yfor(3:1:xdim-3,1)+prev_y_minus_1for(3:1:xdim-3,1))+c3efffor*(prev_y_minus_1for(2:1:xdim-4,1)+prev_y_minus_1for(4:1:xdim-2,1)+prev_yfor(2:1:xdim-4,1)+prev_yfor(4:1:xdim-2,1));
        end
        
        %Storage vectors for boundary and boundary-1 values of previous and its
        %previous time steps updated at upward boundary
        prev_prev_yfor=prev_yfor;
        prev_prev_y_minus_1for=prev_y_minus_1for;
        prev_yfor(1:1:xdim,1)=Ez(1:1:xdim,ydim-2);
        prev_y_minus_1for(1:1:xdim,1)=Ez(1:1:xdim,ydim-3);
        
        %Mur's abc conditions obtained from Mur's difference equation for
        %downward boundary
        if n>1
            Ez(3:1:xdim-3,2)=-prev_prev_yrev(3:1:xdim-3,1)+c1effrev*(Ez(3:1:xdim-3,3)+prev_prev_y_minus_1rev(3:1:xdim-3,1))+c2effrev*(prev_yrev(3:1:xdim-3,1)+prev_y_minus_1rev(3:1:xdim-3,1))+c3effrev*(prev_y_minus_1rev(2:1:xdim-4,1)+prev_y_minus_1rev(4:1:xdim-2,1)+prev_yrev(2:1:xdim-4,1)+prev_yrev(4:1:xdim-2,1));
        end
        
        %Storage vectors for boundary and boundary-1 values of previous and its
        %previous time steps updated at downward boundary
        prev_prev_yrev=prev_yrev;
        prev_prev_y_minus_1rev=prev_y_minus_1rev;
        prev_yrev(1:1:xdim,1)=Ez(1:1:xdim,3);
        prev_y_minus_1rev(1:1:xdim,1)=Ez(1:1:xdim,2);
        
        %Mirroring of corner values taking the fact that corners are reached by the fields from the previous corners
        %in two time steps as S=1/sqrt(2) viz. sqrt(2)*delta(distance between two corners) is reached in 2 time steps
        Ez(2,2)=prev_prev_xrev(3);
        Ez(2,ydim-2)=prev_prev_xrev(ydim-3);
        Ez(xdim-2,2)=prev_prev_x_minus_1for(3);
        Ez(xdim-2,ydim-2)=prev_prev_x_minus_1for(ydim-3);
        
        % Source condition incorporating given free space wavelength 'free_space_wavelength'
        % and having a location at the left half of the top edge (in glass) of the domain just
        % after the absorbing boundary with the plane of propagation of the emanating wave 
        % inclined at an angle 'alpha' with x-axis
        N_lambda=free_space_wavelength/delta;
        t_start=1;
        alpha=55; %Angle 'alpha' should be < 42 degrees for refraction
        i=75:1:125;
        Ez(i,2)=sin(((2*pi*c*S*index)/(N_lambda*delta))*(n-t_start)*deltat*(i./i)-((2*pi*(i-76)*index*sin(alpha*pi/180))/(N_lambda)));
        
        %Movie type colour scaled image plot of Ez
        h=imagesc(1e+6*delta*(1:1:xdim),1e+6*(delta*(1:1:ydim))',Ez',[-1,1]);colorbar;
        set(h,'AlphaData',10*epsilon'/epsilon0);
        title(['\fontsize{20}Color-scaled image plot of Ez to see refraction at glass/air interface with ABC boundary at time = ',num2str(n*deltat*1e+15),' fs']);
        xlabel('x in microns','FontSize',20);
        ylabel('y in microns','FontSize',20);
        set(gca,'FontSize',20);
        getframe;
    end
end

% [1] "Computational Electrodynamics - The Finite Difference Time Domain
%      Method" - Allen Taflove, Susan.c.Hagness, Third Edition, Artech House