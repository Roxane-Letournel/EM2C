function [normalized_eigenvalues] = KLE_Burger1D(omega1,couplage,N,Tend,Ns)

% Ns = 3;
% N = 10000;

C = Scalar_Field_Burger1D(omega1,couplage,N,Tend);

Csim = zeros(Ns,size(C,1),size(C,2));

parfor k = 1:Ns
    Csim(k,:,:) = Scalar_Field_Burger1D(omega1,couplage,N,Tend);
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

