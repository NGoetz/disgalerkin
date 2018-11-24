function [rhos,rhous,Eners] = SlopeLimitNc(rho, rhou, Ener);

% function ulimit = SlopeLimitN(u);
% Purpose: Apply slopelimiter (Pi^N) to u assuming u an N'th orher polynomial            
gamma=1.4;
Globals1D;

rhoh = invV*rho; rhoh(2:Np,:)=0; rhoavg = V*rhoh; vrho = rhoavg(1,:);
rhouh = invV*rhou; rhouh(2:Np,:)=0; rhouavg = V*rhouh; vrhou = rhouavg(1,:);
Enerh = invV*Ener; Enerh(2:Np,:)=0; Eneravg = V*Enerh; vEner = Eneravg(1,:);


ph=(gamma-1.0)*(vEner - 0.5*((vrhou).^2)./vrho);
ch=sqrt(gamma*ph./vrho);
Hh=(vEner+ph)./vrho;
uvh=vrhou./vrho;

aa=((gamma-1)/4)*((uvh.^2)./(ch.^2))+(uvh./(2*ch)); 
ab=-((gamma-1)./(2*ch.^2)).*uvh-(1)./(2*ch); 
ac=(gamma-1)./(2*(ch.^2)); 
ba=1-((gamma-1)./(2*ch.^2)).*(uvh.^2);
bb=(gamma-1)*uvh./(ch.^2); 
bc=-2*ac; 
ca=((gamma-1)/4)*((uvh.^2)./ch.^2)-uvh./(2*ch);
cb=-((gamma-1)./(2*ch.^2)).*(uvh)+(1)./(2*ch); 
cc=ac;

maa=mab=mac=ones(1,K);
mba=uvh-ch;
mbb=uvh;
mbc=uvh+ch;
mca=Hh-uvh.*ch;
mcb=0.5*uvh.^2;
mcc=Hh+uvh.*ch;


rhon=zeros(Np,K);
rhoun=zeros(Np,K);
Enern=zeros(Np,K);
for i=1:Np
  
  rhon(i,:)=aa.*rho(i,:)+ab.*rhou(i,:)+ac.*Ener(i,:);
  rhoun(i,:)=ba.*rho(i,:)+bb.*rhou(i,:)+bc.*Ener(i,:);
  Enern(i,:)=ca.*rho(i,:)+cb.*rhou(i,:)+cc.*Ener(i,:);
end

h = x(Np,:)-x(1,:); 

for i=1:3
  if i==1
    u=rhon;
  elseif i==2
    u=rhoun;
  elseif i==3
    u=Enern;
  end
  % Compute cell averages
  uh = invV*u; uh(2:Np,:)=0; uavg = V*uh; v = uavg(1,:);

  % Apply slope limiter as needed.
  ulimit = u; eps0=1.0e-8;

  % find end values of each element
  ue1 = u(1,:); ue2 = u(end,:);

  % find cell averages
  vk = v; vkm1 = [v(1),v(1:K-1)]; vkp1 = [v(2:K),v(K)]; 

  % Apply reconstruction to find elements in need of limiting

  ve1 = vk - minmodB([(vk-ue1);vk-vkm1;vkp1-vk],2000,h);
  ve2 = vk + minmodB([(ue2-vk);vk-vkm1;vkp1-vk],2000,h);

  ids = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);

  % Check to see if any elements require limiting
  if(~isempty(ids))
    % create piecewise linear solution for limiting on specified elements
    uhl = invV*u(:,ids); uhl(3:Np,:)=0; ul = V*uhl;
    
    % apply slope limiter to selected elements
    ulimit(:,ids) = SlopeLimitLin(ul,x(:,ids),vkm1(ids),vk(ids),vkp1(ids));
  end
  if i==1
    rho1=ulimit;
  elseif i==2
    rhou1=ulimit;
  elseif i==3
    Ener1=ulimit;
  end
end
for i=1:Np
  
   rhos(i,:)=maa.*rho1(i,:)+mab.*rhou1(i,:)+mac.*Ener1(i,:);
  rhous(i,:)=mba.*rho1(i,:)+mbb.*rhou1(i,:)+mbc.*Ener1(i,:);
  Eners(i,:)=mca.*rho1(i,:)+mcb.*rhou1(i,:)+mcc.*Ener1(i,:);
end

return;
