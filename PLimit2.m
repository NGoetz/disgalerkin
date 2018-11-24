function [rhol,rhoul,Enerl]=PLimit2(rho,rhou,Ener)

% function ulimit = PLimit2(rho,rhou,Ener);
% Purpose: Apply positivity preserving limiter

Globals1D;
gamma=1.4;

rhoh = invV*rho; rhoh(2:Np,:)=0; rhoavg = V*rhoh; vrho = rhoavg(1,:);
rhouh = invV*rhou; rhouh(2:Np,:)=0; rhouavg = V*rhouh; vrhou = rhouavg(1,:);
Enerh = invV*Ener; Enerh(2:Np,:)=0; Eneravg = V*Enerh; vEner = Eneravg(1,:);
rhol1=zeros(Np,K);


eps=min([10^-13,min(vrho),min((gamma-1.0)*(vEner-0.5*(vrhou.^2)./vrho))]);

for j=1:K
 
 minrho=10^20;
 for alpha=1:Np
    if rho(alpha,j)<minrho
      minrho=rho(alpha,j);
    end
  end
    theta1=min(vrho(j)-eps,vrho(j)-minrho);
    rhol1(:,j)=(theta1*(rho(:,j)-vrho(j)))./(vrho(j)-minrho+10^-13)+vrho(j);
    for alpha=1:Np
    if(rhol1(alpha,j)<0) 
      p=(gamma-1.0)*(Ener - 0.5*((rhou).^2)./rho);
    
      theta1
      minrho;
      eps;
      vrho(j);
    end
   end
end

p=(gamma-1.0)*(Ener - 0.5*((rhou).^2)./rhol1);

tf=1;
for j=1:K
  
  for alpha=1:Np
    if(p(alpha,j)<eps)
            %rhol1;
            diffrho = rhol1(alpha,j)-vrho(j);
            diffrhou=rhou(alpha,j)-vrhou(j);
            diffEner = Ener(alpha,j)-vEner(j);
            a1 = 2.0*diffrho*diffEner - diffrhou^2;
            b1 = 2.0*diffrho*(vEner(j) - eps/(gamma-1.0))+ 2.0*vrho(j)*diffEner- 2.0*vrhou(j)*diffrhou;
            c1 = 2.0*vrho(j)*vEner(j)- vrhou(j)^2- 2.0*eps*vrho(j)/(gamma-1.0);
            b1=b1/a1 ;
            c1=c1/a1;
            sr = sqrt(abs(b1*b1 - 4.0*c1) );
            t1 = 0.5*(-b1 - sr);
            t2 = 0.5*(-b1 + sr);
    
       if(1*10^-14<t1&&t1<1+1*10^-14)
         ts=t1;
       else
         ts=t2; 
      endif
      
 
    
    else
      ts=1;
    endif
  
    if (ts<tf)
      tf=ts;
    endif
 
 end
 
 for alpha=1:Np
    rhol(alpha,j)=tf*(rhol1(alpha,j)-vrho(j))+vrho(j);
    rhoul(alpha,j)=tf*(rhou(alpha,j)-vrhou(j))+vrhou(j);
    Enerl(alpha,j)=tf*(Ener(alpha,j)-vEner(j))+vEner(j);
    
 end
end



return;