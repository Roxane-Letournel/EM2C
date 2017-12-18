function C = TG1D_Scalar_Field(N)
%PARAMETRES
taup=10;

%Parametres Champ scalaire
D=1d-2;
alpha=0.1;
beta=1;

%Taille du maillage pour la projection
resolution=100; 
x=linspace(0,1,resolution+1); Nx = resolution; dx = 1/Nx;

%Pas de temps
deltat=0.01;
Fo=0.25;
CFL=1.0;
deltat=min([deltat,Fo*dx^2/D,CFL*dx/1 ]);

Tend=0.1;
Npas=Tend/deltat;
% N = 5;


%Initialisation Particules Taylor Green    
X=zeros(N,Npas);
U=zeros(N,Npas);
Ug=zeros(N,Npas);


%initialisation particules
X(:,1)=random('uniform',0,1,N,1);
U(:,1)=random('uniform',0,0,N,1);

%initialisation champ scalaire
C = zeros(Nx,Npas);

for i=1:Npas-1

    for k=1:N
        Ug(k,i)=sin(2*pi*X(k,i));
    end
    
    X(:,i+1) = X(:,i)+deltat*U(:,i);
    U(:,i+1) = U(:,i)+deltat/taup*(Ug(:,i)-U(:,i));
    
    X(:,i+1) = mod( X(:,i+1) , 1);
    
    
    %Densité de particules
    A = histogram(X(:,i),x);
    density= A.Values / N;
    
    % champ scalaire
            
        for k=1:Nx           
            
            kp1 = (k+1<Nx+1)*(k+1)+(k+1==Nx+1)*1;
            km1 = (k-1>0)*(k-1)+(k-1==0)*Nx;

           %advection
              C(k,i+1) = C(k,i) - deltat*(sin(2*pi*kp1/Nx)*C(kp1,i) - sin(2*pi*k/Nx)*C(k,i)) /dx;

           %diffusion
              C(k,i+1) = C(k,i+1) + deltat*D*(C(kp1,i) + C(km1,i) - 2*C(k,i)) /dx^2;
            %source
              C(k,i+1) = C(k,i+1) + deltat*alpha*density(k).^beta;  
            
        end
      
    
end
Ug(:,Npas)=sin(2*pi*X(:,Npas));

return
end

