function [rho,rhou,Ener] = Euler1D(rho, rhou, Ener, FinalTime)

% function [rho, rhou, Ener] = Euler1D(rho, rhou, Ener, FinalTime)
% Purpose  : Integrate 1D Euler equations until FinalTime starting with
%            initial conditions [rho, rhou, Ener]

Globals1D;



% Parameters
gamma = 1.4; CFL = 0.5;%! war 1
 time = 0;

% Prepare for adaptive time stepping
mindx = min(x(2,:)-x(1,:));

% Limit initial solution
[rho,rhou,Ener]=SlopeLimitNc(rho,rhou,Ener);

[rho,rhou,Ener]=PLimit2(rho,rhou,Ener);

% outer time step loop 
while(time<FinalTime)

rhoh = invV*rho; rhoh(2:Np,:)=0; rhoavg = V*rhoh; vrho = rhoavg(1,:);
rhouh = invV*rhou; rhouh(2:Np,:)=0; rhouavg = V*rhouh; vrhou = rhouavg(1,:);
Enerh = invV*Ener; Enerh(2:Np,:)=0; Eneravg = V*Enerh; vEner = Eneravg(1,:);
dt=10^20;
for j=1:K
  v=vrhou(j)/vrho(j);
  p=(gamma-1)*(vEner(j)-0.5*vrhou(j)*v);
  s=sqrt(gamma*p/vrho(j));
  sp=abs(v)+s;
  dt=min(dt,mindx/sp);
end
dt=(1/((N+1)*N))*dt;

  
  % 3rd order SSP Runge-Kutta
  
  % SSP RK Stage 1.
 [rhsrho,rhsrhou,rhsEner]  = EulerRHS1D(rho, rhou, Ener,time, FinalTime);
  rho1  = rho  + dt*rhsrho;
  rhou1 = rhou + dt*rhsrhou;
  Ener1 = Ener + dt*rhsEner;

  % Limit fields
  [rho1,rhou1,Ener1]=SlopeLimitNc(rho1,rhou1,Ener1);

 
[rho1,rhou1,Ener1]=PLimit2(rho1,rhou1,Ener1);

  % SSP RK Stage 2.
  [rhsrho,rhsrhou,rhsEner]  = EulerRHS1D(rho1, rhou1, Ener1, time+dt, FinalTime);
  rho2   = (3*rho  + rho1  + dt*rhsrho )/4;
  rhou2  = (3*rhou + rhou1 + dt*rhsrhou)/4;
  Ener2  = (3*Ener + Ener1 + dt*rhsEner)/4;

  % Limit fields
   [rho2,rhou2,Ener2]=SlopeLimitNc(rho2,rhou2,Ener2);

 [rho2,rhou2,Ener2]=PLimit2(rho2,rhou2,Ener2);
  % SSP RK Stage 3.
  [rhsrho,rhsrhou,rhsEner]  = EulerRHS1D(rho2, rhou2, Ener2, time+0.5*dt, FinalTime);
  rho  = (rho  + 2*rho2  + 2*dt*rhsrho )/3;
  rhou = (rhou + 2*rhou2 + 2*dt*rhsrhou)/3;
  Ener = (Ener + 2*Ener2 + 2*dt*rhsEner)/3;
  
  % Limit solution
  [rho,rhou,Ener]=SlopeLimitNc(rho,rhou,Ener);
 
  [rho,rhou,Ener]=PLimit2(rho,rhou,Ener);
  % Increment time and adapt timestep
  time=time+dt;
 
  end
end
return