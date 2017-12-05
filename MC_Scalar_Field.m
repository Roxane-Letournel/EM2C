clear all
close all
%% Convergence du champ scalaire avec N

Tend = 1;
Stc = 5;
N=[10,100,1000,1000];
CN = zeros(200,200,4);

for k=1:4
    C = Scalar_Field(Stc,N(k),Tend);
    CN(:,:,k) = C(:,:);
end

%Champ scalaire final
for k=1:4
    figure
    imagesc(CN(:,:,k))
    colorbar
    title('Final Scalar Field')
end

%% Etude du RMS en fonctions de Ns
Tend = 1;
Stc = 5;
N=333;
Ns = 30;

for k=1:Ns
    k
    C = Scalar_Field(Stc,N,Tend);
    CN(:,:,k) = C(:,:);
end

for i = 1:200
    for j = 1:200
        m = zeros(1,Ns);
        for k = 1:Ns
            m(k) = CN(i,j,k)^2;
        end
        C_RMS (i,j) = sqrt( mean(m) - mean(CN(i,j,:))^2);
    end
end

figure
imagesc(C_RMS)
colorbar
