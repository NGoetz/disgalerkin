% Driver script for solving the 1D Euler equations
Globals1D;

% Polynomial order used for approximation 
N = 2;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0, 1, 100);

% Initialize solver and construct grid and metric
StartUp1D;
gamma = 1.4;

% Set up initial conditions -- Sedov problem
MassMatrix = inv(V')/V;
cx = ones(Np,1)*sum(MassMatrix*x,1)/2; 

rho=ones(Np,K);

rhou=zeros(Np,K);

Ener=10^-12*ones(Np,K);
Ener(:,K/2)=3200000*1/(N+1);

FinalTime = 0.001;

% Solve Problem
[rho,rhou,Ener] = Euler1D(rho,rhou,Ener,FinalTime);


figure;plot(x,rho);drawnow;
