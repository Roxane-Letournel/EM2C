% clear all
% close all

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
N = 5;


%Initialisation Particules Taylor Green    
X=zeros(N,Npas);
U=zeros(N,Npas);
Ug=zeros(N,Npas);


%initialisation particules
X(:,1)=random('uniform',0,1,N,1);
U(:,1)=random('uniform',0,0,N,1);

%initialisation champ scalaire
C = zeros(Nx,Npas);

figure
plot(linspace(0,1,resolution),C(:,1))
axis([0 1 0 0.03])
ax = gca;
ax.NextPlot = 'replaceChildren';

for i=1:Npas-1

    for k=1:N
        Ug(k,i)=sin(2*pi*X(k,i));
    end
    
    X(:,i+1) = X(:,i)+deltat*U(:,i);
    U(:,i+1) = U(:,i)+deltat/taup*(Ug(:,i)-U(:,i));
    
    X(:,i+1) = mod( X(:,i+1) , 1);
    
    
    %Densit� de particules
    A = histogram(X(:,i),x);
    density= A.Values / N;
    
    % champ scalaire
            
        for k=1:Nx           
            
            kp1 = (k+1<Nx+1)*(k+1)+(k+1==Nx+1)*1;
            km1 = (k-1>0)*(k-1)+(k-1==0)*Nx;

           %advection
              C(k,i+1) = C(k,i) - deltat*(f(C(kp1,i)) - f(C(k,i))) /dx;

           %diffusion
              C(k,i+1) = C(k,i+1) + deltat*D*(C(kp1,i) + C(km1,i) - 2*C(k,i)) /dx^2;
            %source
              C(k,i+1) = C(k,i+1) + deltat*alpha*density(k).^beta;  
            
        end
        plot(linspace(0,1,resolution),C(:,i+1),'LineWidth',2)
        title(['Scalar Field for N = ',num2str(N),' particles and t =',num2str(i*deltat),' s'])
        drawnow
        F(i) = getframe(gcf);
    
end
Ug(:,Npas)=sin(2*pi*X(:,Npas));

% fig = figure;
% movie(fig,F,1)

figure

subplot(1,2,1)
plot(linspace(0,1,resolution),C(:,Npas))
xlabel('Position')
title(['Scalar Field for N = ',num2str(N),' particles'])

subplot(1,2,2)
for k=1:10:N
    scatter(X(k,:),U(k,:),1.5,'b','filled')
    hold on
end
scatter(X(1,:),Ug(1,:),1.5,'r','filled')
% legend('One particle','Gas')
xlabel('Position')
ylabel('Velocity')

% G�n�rer plusieurs simulations de Csim(k,x,t)
% KC(s,t) = mean(Csim(:,x,s).*Csim(:,x,t)

% figure
% M = cov(C);
% pcolor(M)
% 
% figure
% bar(eig(M))

K = zeros(40); %Sur temporel
for s = 1:40
    for t = 1:40
        K(s,t) = mean(Csim(:,5,s).*Csim(:,5,t));
    end
end

K = zeros(100); %sur spatial
for s = 1:100
    for t = 1:100
        K(s,t) = mean(Csim(:,s,Npas).*Csim(:,t,Npas));
    end
end
figure
bar(eig(K))