clearvars
close all

%Paramètres Taylor Green
tau=5/8/3.14;
N=100;

%Parametres Champ scalaire
D=1d-2;
alpha=0.1;
beta=2;
% D = 0;
% alpha = 0;
%Taille du maillage pour la projection
resolution=200; 
x=linspace(0,1,resolution+1); Nx = resolution; dx = 1/Nx;
y=linspace(0,1,resolution+1); Ny = resolution; dy = 1/Ny;

%Pas de temps
deltat=0.01;
Fo=0.25;
CFL=0.4;
deltat=min([deltat,Fo*dx^2/D,CFL*dx/1 ]);

Tend=0.5;
Npas=Tend/deltat;


%Initialisation Particules Taylor Green    
X=zeros(N,2);
U=zeros(N,2);
Y=zeros(N,2);
V=zeros(N,2);
Ug=zeros(N,2);
Vg=zeros(N,2);


%initialisation particules
X(:,1)=random('uniform',0,1,N,1);
U(:,1)=random('uniform',0,0,N,1);
Y(:,1)=random('uniform',0,1,N,1);
V(:,1)=random('uniform',0,0,N,1);

%initialisation champ scalaire
C = zeros(Nx,Ny);
C_old = zeros(Nx,Ny); 

% Center = 0.5*Nx;
% Radius = 0.3*Nx;
% for i = Center - Radius : Center + Radius 
%     j1 = floor(Center - sqrt(Radius^2 - (i-Center)^2));
%     j2 = floor(Center + sqrt(Radius^2 - (i-Center)^2));
%     C(j1:j2,i) = 1;
% end

figure
pcolor(C)
title('Initial scalar field')


%Boucle
for i=1:Npas-1
    for k=1:N
        Ug(k,1)=-sin(2*pi*X(k,1))*cos(2*pi*Y(k,1));
        Vg(k,1)=cos(2*pi*X(k,1))*sin(2*pi*Y(k,1));
    end
    
    X(:,2)=mod(X(:,1)+deltat*U(:,1),1);
    Y(:,2)=mod(Y(:,1)+deltat*V(:,1),1);
    
    
    U(:,2)=U(:,1)+deltat/tau*(Ug(:,1)-U(:,1));
    V(:,2)=V(:,1)+deltat/tau*(Vg(:,1)-V(:,1));
   
    X(:,1) = X(:,2);
    Y(:,1) = Y(:,2);
    U(:,1) = U(:,2);
    V(:,1) = V(:,2);
    
    %Densité de particules
    mX=zeros(N,2);
    mX(:,1)=X(:,1);
    mX(:,2)=Y(:,1);
    A=hist2d(mX,y,x);
    density=transpose(A/N);
    
    % champ scalaire
    C_old(:,:)=C(:,:);
        for k=1:Nx
            for l=1:Ny
                kp1=(k+1<Nx+1)*(k+1)+(k+1==Nx+1)*1;
                km1=(k-1>0)*(k-1)+(k-1==0)*Nx;
                lp1=(l+1<Ny+1)*(l+1)+(l+1==Ny+1)*1;
                lm1=(l-1>0)*(l-1)+(l-1==0)*Ny;
%                 
                Ug_kp1_l = -sin(2*pi*kp1/Nx)*cos(2*pi*l/Ny);
                Ug_km1_l = -sin(2*pi*km1/Nx)*cos(2*pi*l/Ny);
                Vg_k_lp1 = cos(2*pi*k/Nx)*sin(2*pi*lp1/Ny);
                Vg_k_lm1 = cos(2*pi*k/Nx)*sin(2*pi*lm1/Ny);
                
                Ug_k_l = -sin(2*pi*k/Nx)*cos(2*pi*l/Ny);
                Vg_k_l = cos(2*pi*k/Nx)*sin(2*pi*l/Ny);
                
                %C(k,l) = C(yi,xj) ??
%                 %advection

                   %Lax-Windroff
%                    C(k,l)=C_old(k,l)...
%                        -deltat*Ug_k_l*(C_old(kp1,l)-C_old(km1,l))/2/dx...
%                        -deltat*Vg_k_l*(C_old(k,lp1)-C_old(k,lm1))/2/dy...
%                        +deltat^2/2/dx^2*Ug_k_l^2*(C_old(kp1,l)-2*C_old(k,l)+C_old(km1,l))...
%                        +deltat^2/2/dy^2*Vg_k_l^2*(C_old(k,lp1)-2*C_old(k,l)+C_old(k,lm1));
%                    
                   C(k,l)=C_old(k,l)...
                       -deltat*(Ug_kp1_l*C_old(kp1,l)-Ug_km1_l*C_old(km1,l))/2/dx...
                       -deltat*(Vg_k_lp1*C_old(k,lp1)-Vg_k_lm1*C_old(k,lm1))/2/dy;
%                    
%                    C(k,l) = C_old(k,l)...
%                        -deltat*Ug_k_l*(C_old(kp1,l)-C_old(km1,l))/2/dx...
%                        -deltat*Vg_k_l*(C_old(k,lp1)-C_old(k,lm1))/2/dy;
                  
%                     C(k,l)=C_old(k,l)...
%                        -deltat*(Ug_k_l*C_old(k,l)-Ug_km1_l*C_old(km1,l))/dx...
%                        -deltat*(Vg_k_l*C_old(k,l)-Vg_k_lm1*C_old(k,lm1))/dy;
               %diffusion
                  C(k,l)=C(k,l)+deltat*D*(C_old(kp1,l)+C_old(km1,l)+C_old(k,lp1)+C_old(k,lm1)-4*C_old(k,l))/dx^2;
                %source
                  C(k,l)=C(k,l)+alpha*density(l,k).^beta;  
            end
        end
        
        imagesc(C)
        %Figure intermédiaire du champ scalaire
        if i == 600
            figure
           imagesc(C)
            colorbar
            title('Scalar field at i=100')
        end
end

t=linspace(1,Npas,Npas);

% % Particules initiales et finales
% figure
% scatter(X(:,1),Y(:,1),5,'filled')
% title('Initial distribution')
figure
scatter(X(:,1),Y(:,1),5,'filled')
title('Final distribution')


% Densité de particules
figure
imagesc(density)
colorbar
title('Particles density')

%Champ scalaire final
figure
imagesc(C)
colorbar
title('Final Scalar Field')