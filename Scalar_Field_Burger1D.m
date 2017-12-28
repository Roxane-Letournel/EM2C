function C = Scalar_Field_Burger1D(omega1,couplage,N,Tend)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES PHYSIQUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parametres Burger gas
% omega1=2*pi;
omega2=0*pi;
omega3=0*pi;  
nu = 0.01;
mu=1.9*10^(-8); %-5 en vrai mais stokes trop faible sinon
% couplage = 1;

% Paramètres Particules 
densite=10^10;
rho=8900;
dp0=1.2*10^(-5);
taup0=rho*dp0^2/18/mu;
mp0=rho*pi/6*dp0^3*densite; %on multiplie par 2*10^10 pour avoir le meme particle number density

Stokes=taup0*omega1 ;

taup = taup0;
mp = mp0;

%Parametres Champ scalaire
D=1d-2;
alpha=0.1;
beta=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES NUMERIQUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Taille du maillage pour la projection
resolution=100; 
x=linspace(0,1,resolution+1); Nx = resolution; dx = 1/Nx;

%Pas de temps
deltat=0.01;
Fo=0.25;
CFL=1.0;
deltat=min([deltat,Fo*dx^2/D,CFL*dx/1 ]);

% Tend=1;
Npas=Tend/deltat;
% N = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialisation Particules    
X=zeros(N,Npas);
U=zeros(N,Npas);
X(:,1)=random('uniform',0,1,N,1);
U(:,1)=random('uniform',0,1,N,1);
F = zeros(Nx,Npas); %coupling term

%Initialisation Burger
Ug=zeros(Nx,Npas);
u = zeros(Nx,1);
unew = zeros(Nx,1);


%initialisation champ scalaire
C = zeros(Nx,Npas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUCLE PRINCIPALE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:Npas-1

    % Calcul du champ de vitesse à l'instant i Equation de Burger
    
    %  First node.
          unew(1) = u(1);
          
    %  Interior nodes.
        for j=2:Nx-1
              unew(j) =                u(j) + deltat * ( ...
              nu    * (   u(j+1) - 2.0 * u(j) + u(j-1)  ) / dx^2 ...
              - 0.5 * ( f(u(j+1))                - f(u(j-1)) ) / dx ...
              +1*(sin(omega1*x(j))+sin(omega2*x(j))+sin(omega3*x(j))) + F(j,i));
        end

    %  Last node.
           ur = 2.0 * u(Nx) - u(Nx-1);
           unew(Nx) =               u(Nx) + deltat * ( ...
            nu    * (   ur - 2.0 * u(Nx) + u(Nx-1)  ) / dx^2 ...
            -       ( f(ur)            - f(u(Nx-1)) ) / dx ...
            +1*(sin(omega1*x(Nx))+sin(omega2*x(Nx))+sin(omega3*x(Nx))) + F(Nx,i));
        
        u(:) = unew(:);
        Ug(:,i) = u(:);

    
    X(:,i+1) = X(:,i)+deltat*U(:,i);
    X(:,i+1) = mod( X(:,i+1) , 1);
    
    ind = floor(100*X(:,i))+1; % l'indice donne la position de la particule dans le maillage
    for k=1:N
        U(k,i+1) = U(k,i)+deltat/taup*(Ug(ind(k),i)-U(k,i));
    end
    
    if couplage == 1
        for k = 1:N
            j = ind(k);
            F(j,i+1) = F(j,i+1) + mp*(Ug(j,i)-U(k,i))/taup ;
        end
    end
    
    %Densité de particules
    A = histogram(X(:,i),x);
    density= A.Values / N;
    
    % champ scalaire
            
        for k=1:Nx           
            
            kp1 = (k+1<Nx+1)*(k+1)+(k+1==Nx+1)*1;
            km1 = (k-1>0)*(k-1)+(k-1==0)*Nx;

           %advection
              C(k,i+1) = C(k,i) - deltat*(Ug(kp1,i)*C(kp1,i) - Ug(k,i)*C(k,i)) /dx;
           %diffusion
              C(k,i+1) = C(k,i+1) + deltat*D*(C(kp1,i) + C(km1,i) - 2*C(k,i)) /dx^2;
            %source
              C(k,i+1) = C(k,i+1) + deltat*alpha*density(k).^beta;  
            
        end
    
end

return 
end