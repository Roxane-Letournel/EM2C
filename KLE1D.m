function [normalized_eigenvalues] = KLE1D(Ns,N)

% Ns = 3;
% N = 10000;

C = TG1D_Scalar_Field(N);

Csim = zeros(Ns,size(C,1),size(C,2));

parfor k = 1:Ns
    Csim(k,:,:) = TG1D_Scalar_Field(N);
end
    
K = zeros(100); %sur spatial
for s = 1:100
    for t = 1:100
        K(s,t) = mean(Csim(:,s,end).*Csim(:,t,end));
    end
end
normalized_eigenvalues = eig(K)/sum(eig(K));
normalized_eigenvalues = flipud(normalized_eigenvalues)

return 
end

