clear all
close all

%PARAMETRES
taup=10;

%Parametres Champ scalaire
D=1d-2;
alpha=0.1;
beta=1;

%Taille du maillage pour la projection
resolution=100; 
x=linspace(0,18,resolution+1); Nx = resolution; dx = 1/Nx;

%Pas de temps
deltat=0.01;
Fo=0.25;
CFL=1.0;
deltat=min([deltat,Fo*dx^2/D,CFL*dx/1 ]);

Tend=100;
Npas=Tend/deltat;
N = 100;


%Initialisation Particules Taylor Green    
X=zeros(N,Npas);
U=zeros(N,Npas);
Ug=zeros(N,Npas);


%initialisation particules
X(:,1)=random('uniform',0,1,N,1);
U(:,1)=random('uniform',0,0,N,1);

%initialisation champ scalaire
C = zeros(Nx,1);
C_old = zeros(Nx,1); 

for i=1:Npas-1

    for k=1:N
        Ug(k,i)=sin(X(k,i));
    end
    
    X(:,i+1) = X(:,i)+deltat*U(:,i);
    U(:,i+1) = U(:,i)+deltat/taup*(Ug(:,i)-U(:,i));
    
    
    %Densité de particules
    A = histogram(X(:,i),x);
    density= A.Values / N;
    
    % champ scalaire
    C_old(:)=C(:);
    
%     %  First node.
%       C(1) = C_old(1);
%     for j=2:Nx-1
%     %  Interior nodes.
%           C(j) =                C_old(j) + deltat * ( ...
%           D    * (   C_old(j+1) - 2.0 * C_old(j) + C_old(j-1)  ) / dx^2 ...
%           - 0.5 * ( f(C_old(j+1))                - f(C_old(j-1)) ) / dx ...
%           + alpha*density(j).^beta);
%     end
% 
% %  Last node.
%        Cr = 2.0 * C_old(Nx) - C_old(Nx-1);
%        C(Nx) =               C_old(Nx) + deltat * ( ...
%         D    * (   Cr - 2.0 * C_old(Nx) + C_old(Nx-1)  ) / dx^2 ...
%         -       ( f(Cr)            - f(C_old(Nx-1)) ) / dx ...
%         +alpha*density(Nx).^beta);
            
        for k=1:Nx           
            
            kp1 = (k+1<Nx+1)*(k+1)+(k+1==Nx+1)*1;
            km1 = (k-1>0)*(k-1)+(k-1==0)*Nx;

           %advection
              C(k) = C_old(k) - deltat*(f(C(kp1)) - f(C(k))) /dx;


           %diffusion
              C(k) = C(k) + deltat*D*(C_old(kp1) + C_old(km1) - 2*C_old(k)) /dx^2;
            %source
              C(k) = C(k) + deltat*alpha*density(k).^beta;  
            
        end
    
end
Ug(:,Npas)=sin(X(:,Npas));

figure
subplot(1,2,1)
plot(C)
title(['Scalar Field for N = ',num2str(N),' particles'])
subplot(1,2,2)
scatter(X(1,:),U(1,:))
hold on
scatter(X(1,:),Ug(1,:))