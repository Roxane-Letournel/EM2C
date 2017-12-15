clear all
close all

Ns = 50;
N = 10;

C = TG1D_Scalar_Field(N);

Csim = zeros(Ns,size(C,1),size(C,2));

for k = 1:Ns
    Csim(k,:,:) = TG1D_Scalar_Field(N);
end
    
K = zeros(100); %sur spatial
for s = 1:100
    for t = 1:100
        K(s,t) = mean(Csim(:,s,end).*Csim(:,t,end));
    end
end

figure
bar(eig(K))