Ns = 100;
N = 5;

for k = 1:Ns
    Csim(k,:,:) = TG1D_Scalar_Field(N);
end
    
K = zeros(100); %sur spatial
for s = 1:100
    for t = 1:100
        K(s,t) = mean(Csim(:,s,Npas).*Csim(:,t,Npas));
    end
end

figure
bar(eig(K))