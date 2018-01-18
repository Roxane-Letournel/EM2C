function [normalized_eigenvalues] = KLE_Burger1D(omega1,couplage,N,Tend,Ns)

% Ns = 3;
% N = 10000;

C = Scalar_Field_Burger1D(omega1,couplage,N,Tend);

Csim = zeros(Ns,size(C,1),size(C,2));
Csim(1,:,:) = C(:,:);

parfor k = 2:Ns
    k
    Csim(k,:,:) = Scalar_Field_Burger1D(omega1,couplage,N,Tend);
end

champ_moyen = zeros(size(C,1),size(C,2));
for i=1:size(C,1)
    for j = 1:size(C,2)
        champ_moyen(i,j) = mean(Csim(:,i,j));
    end
end

for k = 1:Ns
    for i=1:size(C,1)
        for j = 1:size(C,2)
            Csim(k,i,j) = Csim(k,i,j) - champ_moyen(i,j);
        end
    end
end

K = zeros(100); %sur spatial
for s = 1:100
    for t = 1:100
        K(s,t) = mean(Csim(:,s,end).*Csim(:,t,end));
    end
end
normalized_eigenvalues = eig(K)/sum(eig(K));
normalized_eigenvalues = flipud(normalized_eigenvalues);

return 
end

