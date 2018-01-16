function [normalized_eigenvalues] = KLE1D(Ns,N)

% Ns = 3;
% N = 10000;

C = TG1D_Scalar_Field(N);

Csim = zeros(Ns,size(C,1),size(C,2));

parfor k = 1:Ns
    Csim(k,:,:) = TG1D_Scalar_Field(N);
end

champ_moyen = zeros(size(C,1),size(C,2));
for i=1:size(C,1)
    for j = 1:size(C,2)
        champ_moyen(i,j) = mean(Csim(:,i,j));
    end
end
% x = linspace(1,100,100);
% figure
% plot(x,champ_moyen(:,end))
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
somme = sum(eig(K))
normalized_eigenvalues = flipud(normalized_eigenvalues);

return 
end

