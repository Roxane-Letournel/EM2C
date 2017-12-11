close all
clear all

Tend = 4;
Stc = 5;
N = 10;

Ns = 5;


for k=1:Ns
    C=Scalar_Field(Stc,N,Tend);
    Ck(k,:,:) = C(:,:);
end
%% 

[coeff,score,latent] = pca(Ck(:,:,2));

figure
bar(latent/sum(latent))
xlim([0.5 5])
xlabel(['sum of the 3 first eigenvalues = ',num2str((latent(1)+latent(2)+latent(3))/sum(latent))])
title(['Principal Analysis Components for Ns =',num2str(Ns)])
