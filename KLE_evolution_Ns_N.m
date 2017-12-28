%% Eigenvalues for Ns and N
Ns = 30;
N = 50;
% normalized_eigenvalues = KLE1D(Ns,N); % POur Taylor Green
normalized_eigenvalues = KLE_Burger1D(2*pi,1,N,0.2,Ns); % Pour Burger1D

figure
bar(normalized_eigenvalues(1:Ns+10))
title(['Eigenvalues for Ns =',num2str(Ns),' and Np =',num2str(N)])
xlabel(['Sum of 3 first eigenvalues =',num2str(sum(normalized_eigenvalues(1:3)))])


%% Evolution of eigenvalues for N fixed --> Peu important
N = 1;
Ns = [3,5,10,50,500];
Eigenvalues = zeros(5,100,1);

for k=1:5
    Eigenvalues(k,:,:) = KLE1D(Ns(k),N);
end

figure
for k=1:50
    plot(Ns,Eigenvalues(:,k),'LineWidth',2)
    hold on
end
xlabel('Number of simulation')
ylabel('Eigenvalues')
title(['Evolution of eigenvalues for N =',num2str(N)])


%% Evolution of eigenvalues for Ns fixed
Ns = 20;
N = [1,10,100,1000,10000,100000];
Eigenvalues = zeros(6,100,1);

for k=1:6
    Eigenvalues(k,:,:) = KLE1D(Ns,N(k));
end

figure
for k=1:Ns
    plot(log10(N),Eigenvalues(:,k),'LineWidth',2)
    hold on
end
xlabel('Number of particle in log scale')
ylabel('Eigenvalues')
title(['Evolution of eigenvalues for Ns =',num2str(Ns)])


