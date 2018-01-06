function [normalized_eigenvalues] = KLE_Turbulent2D(N,Ns,T_end)

% Ns = 3;
% N = 10000;

[x,y,C] = Scalar_Field_turbulent2D(N,T_end);

Csim = zeros(Ns,size(C,1),size(C,2));
Csim(1,:,:) = C(:,:);

parfor k = 2:Ns
    [x,y,C] = Scalar_Field_turbulent2D(N,T_end);
    Csim(k,:,:) = C(:,:);
end
    
K = zeros(size(C,1),size(C,2)); %sur spatial

for s = 1:size(C,1)
    for t = 1:size(C,2)
        i = (s-1)*size(C,2) + t;
        
        for u = 1:size(C,1)
            for v = 1:size(C,2)
                j = (u-1)*size(C,2) + v ;
                
                if j>=i
                    K(i,j) = mean(Csim(:,s,t).*Csim(:,u,v));
                    K(j,i) = K(i,j);
                end
                
            end
        end
    end
end
normalized_eigenvalues = eig(K)/sum(eig(K));
normalized_eigenvalues = flipud(normalized_eigenvalues);

return 
end

