function C = Scalar_Field(Stc,N,Tend)

%Paramètres Taylor Green
tau=Stc/8/3.14;
% N=10000;

%Parametres Champ scalaire
D=1d-2;
alpha=0.1;
beta=2;

%Taille du maillage pour la projection
resolution=200; 
x=linspace(0,1,resolution+1); Nx = resolution; dx = 1/Nx;
y=linspace(0,1,resolution+1); Ny = resolution; dy = 1/Ny;

%Pas de temps
deltat=0.01;
Fo=0.25;
CFL=1.0;
deltat=min([deltat,Fo*dx^2/D,CFL*dx/1 ]);

% Tend=4;
Npas=Tend/deltat;


%Initialisation Particules Taylor Green    
X=zeros(N,Npas);
U=zeros(N,Npas);
Y=zeros(N,Npas);
V=zeros(N,Npas);
Ug=zeros(N,Npas);
Vg=zeros(N,Npas);


%initialisation particules
X(:,1)=random('uniform',0,1,N,1);
U(:,1)=random('uniform',0,0,N,1);
Y(:,1)=random('uniform',0,1,N,1);
V(:,1)=random('uniform',0,0,N,1);

%initialisation champ scalaire
C = zeros(Nx,Ny);
C_old = zeros(Nx,Ny); 


%Boucle
for i=1:Npas-1
    for k=1:N
        Ug(k,i)=-sin(2*pi*X(k,i))*cos(2*pi*Y(k,i));
        Vg(k,i)=cos(2*pi*X(k,i))*sin(2*pi*Y(k,i));
    end
    
    X(:,i+1)=mod(X(:,i)+deltat*U(:,i),1);
    Y(:,i+1)=mod(Y(:,i)+deltat*V(:,i),1);
    
    
    U(:,i+1)=U(:,i)+deltat/tau*(Ug(:,i)-U(:,i));
    V(:,i+1)=V(:,i)+deltat/tau*(Vg(:,i)-V(:,i));
   
    %Densité de particules
    mX=zeros(N,2);
    mX(:,1)=X(:,i);
    mX(:,2)=Y(:,i);
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
                
                Ug_kp1_l = -sin(2*pi*kp1/Nx)*cos(2*pi*l/Ny);
                Ug_km1_l = -sin(2*pi*km1/Nx)*cos(2*pi*l/Ny);
                Vg_k_lp1 = cos(2*pi*k/Nx)*sin(2*pi*lp1/Ny);
                Vg_k_lm1 = cos(2*pi*k/Nx)*sin(2*pi*lm1/Ny);
                
                Ug_k_l = -sin(2*pi*k/Nx)*cos(2*pi*l/Ny);
                Vg_k_l = cos(2*pi*k/Nx)*sin(2*pi*l/Ny);
                
                %C(k,l) = C(yi,xj) ??
%                 %advection
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
        
end

return
end