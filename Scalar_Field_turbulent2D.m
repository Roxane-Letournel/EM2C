function [x,y,Ynew] = Scalar_Field_turbulent2D(N,T_end)

% N = 1000;
lll=6; %6:8
    
    
scheme=2; %0=upwind, 1=LW, 2=RK4-DF4, 3=Upwind second order

%gas field
L=3;

load Ugfull.dat
load Vgfull.dat
 Nxin=length(Ugfull);
%Choose mesh

Nx=2^lll+1;
ratio=(Nxin-1)/(Nx-1);
Ug=zeros(Nx,Nx);
Vg=zeros(Nx,Nx);


Ug(1:Nx,1:Nx)=Ugfull(1:ratio:Nxin,1:ratio:Nxin);
Vg(1:Nx,1:Nx)=Vgfull(1:ratio:Nxin,1:ratio:Nxin);

 dx=L/(Nx-1);
% 
 for i=1:Nx
     for j=1:Nx
         x(i,j)=(i-1)*dx;
         y(i,j)=(j-1)*dx;
     end
 end


% T_end=15;



Umax=max(max([abs(Ug); abs(Vg)]));


CFL=0.4;
Fo=0.2;

%molecular diffusivity
D=1.d-3;
dt=min([Fo*dx^2/D,CFL*dx/Umax]);
ReCell=Umax*dx/D;


Nt=floor(T_end/dt);

Ynew=zeros(Nx,Nx);
Yold=zeros(Nx,Nx);

Yint1=zeros(Nx,Nx);
Yint2=zeros(Nx,Nx);
Yint3=zeros(Nx,Nx);

% Position initiale: particules en cercle
% posinitx=1.5;%rand(1,1)*3;
% posinity=1.5;%rand(1,1)*3;
% r=sqrt((x-posinitx).^2+(y-posinity).^2);
% rlim=0.5;
% rsave(1,l)=rlim;
% for i=1:Nx
%     for j=1:Nx
%         if (r(i,j)<rlim);
%             Yold(i,j)=1.0d0;
%         end
%     end
% end

% Position initiale: particules aléatoires
for part=1:N
    i = floor(rand*Nx)+1;
    j = floor(rand*Nx)+1;
    Yold(i,j) = 1.d0;
end

mass(1)=sum(sum(Yold))/Nx^2;
RMS(1)=std(std(Yold))/mean(mean(Yold));
Ynew=Yold;

t(1)=0;


for k=1:Nt
    %%RK4-DF4

    %%step 1 
    for i=1:Nx
        for j=1:Nx
            ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
            jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
            im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
            jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
            ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
            jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
            im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
            jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


            dYdx=primederivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(im1,j),Yold(im2,j),dx);
            dYdy=primederivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,jm1),Yold(i,jm2),dx);
            d2Ydx2=secondderivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(i,j),Yold(im1,j),Yold(im2,j),dx);
            d2Ydy2=secondderivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,j),Yold(i,jm1),Yold(i,jm2),dx);

            %LAX WENDROFF + CENTERED
            Yint1(i,j) = Yold(i,j)   +0.5d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                        + D*dt*(d2Ydx2+d2Ydy2 )  );      
        end
    end  
    %%step 2
    for i=1:Nx
        for j=1:Nx
            ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
            jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
            im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
            jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
            ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
            jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
            im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
            jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


            dYdx=primederivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(im1,j),Yint1(im2,j),dx);
            dYdy=primederivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,jm1),Yint1(i,jm2),dx);
            d2Ydx2=secondderivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(i,j),Yint1(im1,j),Yint1(im2,j),dx);
            d2Ydy2=secondderivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,j),Yint1(i,jm1),Yint1(i,jm2),dx);

            %LAX WENDROFF + CENTERED
            Yint2(i,j) = Yold(i,j)   +0.5d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                        + D*dt*(d2Ydx2+d2Ydy2 )  );     
        end
    end 
    %%step 3
    for i=1:Nx
        for j=1:Nx
            ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
            jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
            im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
            jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
            ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
            jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
            im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
            jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


            dYdx=primederivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(im1,j),Yint2(im2,j),dx);
            dYdy=primederivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,jm1),Yint2(i,jm2),dx);
            d2Ydx2=secondderivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(i,j),Yint2(im1,j),Yint2(im2,j),dx);
            d2Ydy2=secondderivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,j),Yint2(i,jm1),Yint2(i,jm2),dx);

            %LAX WENDROFF + CENTERED
            Yint3(i,j) = Yold(i,j)   +(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                        + D*dt*(d2Ydx2+d2Ydy2 )  );     
        end
    end 
    %%step 4 (final)
    for i=1:Nx
        for j=1:Nx
            ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
            jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
            im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
            jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
            ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
            jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
            im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
            jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


            dYdx=primederivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(im1,j),Yold(im2,j),dx);
            dYdy=primederivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,jm1),Yold(i,jm2),dx);
            d2Ydx2=secondderivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(i,j),Yold(im1,j),Yold(im2,j),dx);
            d2Ydy2=secondderivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,j),Yold(i,jm1),Yold(i,jm2),dx);

            %LAX WENDROFF + CENTERED
            Ynew(i,j) = Yold(i,j)   +1.d0/6.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                        + D*dt*(d2Ydx2+d2Ydy2 )  ); 


            dYdx=primederivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(im1,j),Yint1(im2,j),dx);
            dYdy=primederivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,jm1),Yint1(i,jm2),dx);
            d2Ydx2=secondderivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(i,j),Yint1(im1,j),Yint1(im2,j),dx);
            d2Ydy2=secondderivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,j),Yint1(i,jm1),Yint1(i,jm2),dx);


            Ynew(i,j) = Ynew(i,j)   +1.d0/3.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                        + D*dt*(d2Ydx2+d2Ydy2 )  ); 

            dYdx=primederivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(im1,j),Yint2(im2,j),dx);
            dYdy=primederivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,jm1),Yint2(i,jm2),dx);
            d2Ydx2=secondderivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(i,j),Yint2(im1,j),Yint2(im2,j),dx);
            d2Ydy2=secondderivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,j),Yint2(i,jm1),Yint2(i,jm2),dx);

            Ynew(i,j) = Ynew(i,j)   +1.d0/3.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                        + D*dt*(d2Ydx2+d2Ydy2 )  );

            dYdx=primederivative_4th(Yint3(ip2,j),Yint3(ip1,j),Yint3(im1,j),Yint3(im2,j),dx);
            dYdy=primederivative_4th(Yint3(i,jp2),Yint3(i,jp1),Yint3(i,jm1),Yint3(i,jm2),dx);
            d2Ydx2=secondderivative_4th(Yint3(ip2,j),Yint3(ip1,j),Yint3(i,j),Yint3(im1,j),Yint3(im2,j),dx);
            d2Ydy2=secondderivative_4th(Yint3(i,jp2),Yint3(i,jp1),Yint3(i,j),Yint3(i,jm1),Yint3(i,jm2),dx);

            Ynew(i,j) = Ynew(i,j)   +1.d0/6.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                        + D*dt*(d2Ydx2+d2Ydy2 )  );                               
        end
    end
    
    t(k+1)=t(k)+dt;
    t(k+1);
    Yold=Ynew;
%     mass(k+1)=sum(sum(Ynew))/Nx^2;
%     RMS(k+1)=std(Ynew(:))/mean(mean(Ynew));
%     if RMS(k+1)<0.05
%         fprintf('mixing reached \n')
%         fprintf('time=%f s\n',t(k+1))
%         break
%     end

end
return 
end







