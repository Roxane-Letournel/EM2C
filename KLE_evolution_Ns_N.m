%% Eigenvalues for Ns and N
Ns = 10;
N = 100;
% T_end = 15;
figure
normalized_eigenvalues = KLE1D(Ns,N); % Pour Taylor Green 1D
% normalized_eigenvalues = KLE_Burger1D(2*pi,1,N,0.2,Ns); % Pour Burger1D
% [normalized_eigenvalues] = KLE_Turbulent2D(N,Ns,T_end); % Pour champ 2D turbulent

figure
bar(normalized_eigenvalues(1:30))
title(['Eigenvalues for Ns =',num2str(Ns),' and Np =',num2str(N)])
xlabel(['Sum of 3 first eigenvalues =',num2str(sum(normalized_eigenvalues(1:3)))])
% xlabel(['1st eigenvalue =',num2str(normalized_eigenvalues(1)),', 2nd eigenvalue =',num2str(normalized_eigenvalues(2)),', 3rd eigenvalue =',num2str(normalized_eigenvalues(3))])

x = linspace(1,100,100)
y = normalized_eigenvalues

% 
figure(7)
hold on
plot(x,y)
hold off
%% Evolution of eigenvalues for N fixed --> Peu important
N = 100;
Ns = [4,10,20,30,50,100,200,400,600,1000];
Np = 10000/Ns
Eigenvalues = zeros(size(Ns,2),100,1);

for k=1:size(Ns,2)
    k
%     Eigenvalues(k,:,:) = KLE1D(Ns(k),N);
%     Eigenvalues(k,:,:) = KLE_Burger1D(2*pi,1,N,0.2,Ns(k));
    Eigenvalues(k,:,:) = KLE_Turbulent2D(N,Ns(k),T_end);
end

figure
for k=1:50
    semilogx(Ns,Eigenvalues(:,k),'LineWidth',2)
    hold on
end
xlabel('Number of simulation')
ylabel('Eigenvalues')
title(['Evolution of eigenvalues for N =',num2str(N)])


%% Evolution of eigenvalues for Ns fixed
Ns = 10;
N = [1,10,100,1000,10000];
Eigenvalues = zeros(size(N,2),size(normalized_eigenvalues,1),1);

for k=1:size(N,2)
%     Eigenvalues(k,:) = KLE1D(Ns,N(k));
%     Eigenvalues(k,:,:) = KLE_Burger1D(2*pi,1,N(k),0.2,Ns);
    normalized_eigenvalues = KLE_Turbulent2D(N(k),Ns,T_end);
    Eigenvalues(k,:,:) = normalized_eigenvalues(:,:);
end

figure
for k=1:size(N,2)
    semilogx(N,Eigenvalues(:,k),'LineWidth',2)
    hold on
end
xlabel('Number of particles')
ylabel('Relative weight of eigenvalues')
title(['Evolution of eigenvalues for Ns = ',num2str(Ns)])
legend('1st eigenvalue','2nd eigenvalue','3rd eigenvalue','4th eigenvalue','5th eigenvalue')

%% Evolution of eigenvalues with Ns and N
% N = 100;
Ns = [4,10,20,50,80,100,200,400,600,800];
Np = Ns;
for i=1:size(Ns,2)
    Np(i) = 10000/Ns(i);
end

Eigenvalues = zeros(size(Ns,2),size(normalized_eigenvalues,1),1);

for k=1:size(Ns,2)
    k
    Eigenvalues(k,:,:) = KLE1D(Ns(k),Np(k));
%     Eigenvalues(k,:,:) = KLE_Burger1D(2*pi,1,N,0.2,Ns(k));
%     Eigenvalues(k,:,:) = KLE_Turbulent2D(Np(k),Ns(k),T_end);
end

figure
for k=1:50
    semilogx(Ns,Eigenvalues(:,k),'LineWidth',2)
    hold on
end
xlabel('Number of simulation')
ylabel('Eigenvalues')
title(['Evolution of eigenvalues for N*Ns =',num2str(N(1)*Ns(1))])

figure
for k=1:50
    plot(Ns,Eigenvalues(:,k),'LineWidth',2)
    hold on
end
xlabel('Number of simulation')
ylabel('Eigenvalues')
title(['Evolution of eigenvalues for N =',num2str(N(1)*Ns(1))])