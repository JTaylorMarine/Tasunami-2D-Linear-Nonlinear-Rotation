clear all; clc; close all;
               %% Block 1: Input parametrs  %%%%%%%%%%%
               %  Block 1: Input parametrs  %%%%%%%%%%%
g=9.8;                                      % acceleration due to gravity                
fi=0;                                       % Latitude (deg)
Omega=0.0000729;                            % Angular velocity of the Earth rotation (1/rad)
f=2*Omega*sind(fi);                         % Coriolis parameter
               % Parameters of the earthquake 
A=1;                                        % Amplitude of surface displacement (m)
La=200000;                                  % Horizontal scale of the earthquake ()
               % Parameters of the basin
Lx=1500000;                                 % Length of model domain (m)
Lx1=400000;                                 % Beginning of the slope (m)  
Lx2=1000000;                                % End of the slope (m)
H1=500;                                     % Depth in the open part (m) 
H2=50;                                      % Depth on the shelf (m) 
               % Parameters for grid
dx=5000;                                    % Horizontal step (m)
x=[0:dx:Lx];                                % Grid vector 
Jx=length(x);                               % Number of grid points
dt=dx/sqrt(2*g*H2)/5;                       % Temporal step (sec)
Nt=120;                                     % Number of time steps    
              %%%%%%%%%%^%%%%% End of Block 1

             %% Block2: Definition of the bottom 
             %  Block2: Definition of the bottom 
for j=1:Jx                              % Loop for bottom definition
    H(j)=H1;                            % Deep part
    if x(j)>Lx1 & x(j)<=Lx2             % Beginning of the continental slope 
        H(j)=H1-(H1-H2)*(x(j)-Lx1)/(Lx2-Lx1);   % Linear topography from Lx1 to Lx2
    end                                 % End of the continental slope 
    if x(j)>Lx2                         % Start of the shelf
        H(j)=H2;                        % Shelf
    end                                 % End of the shelf
end             % End of the loop for bottom definition 
                %%%%%%%%%%%%%%% End of Block 2

             %% Block3: Definition of initial conditions for 
for j=1:Jx                              % Loop for initial conditions 
    if (j-1)*dx < La                    % Initial free surface displacement
        z(1,j)=A*(cos(pi*(j-1)*dx/2/La))^2;  % in the centre of the earthquake 
    else
        z(1,j)=0;                       % Free surface displacement beyond 
    end                                 % the epicentre
    u(1,j)=0;                           % Eastward velocity 
end
               %%%%%%%%%%%%%% End of Block 3 
  
              %% Block4:  First step FTCS: n=1
              %  Block4:  First step FTCS: n=1
    for n=1:Nt-1
       
    for j=2:Jx-1                      % Loop for u,v,z at first time dtep
       u(n+1,j)=u(n,j)-g*dt*(z(n,j+1)-z(n,j-1))/2/dx;
       z(n+1,j)=z(n,j)-dt*(H(j+1)*u(n,j+1)-H(j-1)*u(n,j-1))/2/dx;
    end                          % end of the loop
%%%%%%%%%%%%%%%%% Boundary conditions at first temporal step    
    u(n+1,1)=0;                          % In the epicentre j=1
    z(n+1,1)=z(n+1,2);                   % In the epicentre j=1 
    u(n+1,Jx)=0;                         % Right boundary on the shelf j=Jx
    z(n+1,Jx)=0;                         % Right boundary on the shelf j=Jx  
  
    end
              %%%%%%%%%%%%%%%%  End of Block 4

                       %% Block 6: Visualization
                       % Block 6: Visualization
fig1=figure(1);
numframes=Nt;                % Number of frames equals to number of time steps
B=moviein(numframes,fig1);   % Initialize movie frame memory
set(fig1,'NextPlot','replacechildren')
for i=1:5:numframes
subplot (2,1,1), plot(x/1000,z(i,:),'LineWidth',2)  % plot elevation
set(gca,'XTick',[0:100:Lx/1000],'YTick',[-A:A/5:A]); 
set(gca,'XMinorTick','on')
grid on                                  % Activation of a grid
ylim([-A A]);                            % Limits for y-axis    
xlabel('Distance (km)');                 % Create x-label
ylabel('Depth (m)');                     % Create y-label
%%%%%%%%%%%%%%%%%%%%%%%%%                    
subplot (2,1,2), plot(x/1000,-H,'LineWidth',6,'Color',[0.5 0.5 0.5]) %plot bottom 
set(gca,'XTick',[0:100:Lx/1000]);
ylim([-H1 0]);                          % Limits for y-axis   
xlabel('Distance (km)');                % Create x-label
ylabel('Depth (m)');                    % Create y-label
B(:,i)=getframe(fig1);                  % Get movie frame
end
                     %%%%%%%% End of Block 6

%% Block7:  Nonlinear model CTCS: n>1
              
    for n=2:Nt-1
       
    for j=2:Jx-1                      % Loop for u,v,z at first time dtep
       u(n+1,j)=u(n-1,j)-2*g*dt*(z(n,j+1)-z(n,j-1))/2/dx-u(n,j)*(u(n,j+1)-u(n,j-1))*dt/dx;
       z(n+1,j)=z(n-1,j)-2*dt*((H(j+1)+z(n,j+1))*u(n,j+1)-(H(j-1)+z(n,j-1))*u(n,j-1))/2/dx;
    end                          % end of the loop
%%%%%%%%%%%%%%%%% Boundary conditions at first temporal step    
    u(n+1,1)=0;                          % In the epicentre j=1
    z(n+1,1)=z(n+1,2);                   % In the epicentre j=1 
    u(n+1,Jx)=0;                         % Right boundary on the shelf j=Jx
    z(n+1,Jx)=0;                         % Right boundary on the shelf j=Jx  
  
    end
              %%%%%%%%%%%%%%%%  End of Block 7
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
