close all
clear all

Ns = 100;
figure
for k=1:Ns
    Xp0 = rand*7;
    Up0 = rand;
    [Xp,Up,Ug]=TG1D(Xp0,Up0,10);
    Xpk(k,:) = Xp(:);
    Upk(k,:) = Up(:);
    plot(Xp,Up)
    hold on
end
plot(Xp,Ug)
xlabel('Position')
ylabel('Velocity')

x=linspace(1,Ns-1,Ns-1);
[coeff,score,latent] = pca(Xpk);

figure
bar(latent/sum(latent))
xlim([0.5 5])
xlabel(['sum of the 3 first eigenvalues = ',num2str((latent(1)+latent(2)+latent(3))/sum(latent))])
title(['Principal Analysis Components for Ns =',num2str(Ns)])

% figure
% plot(Xp,Up)
% hold on
% plot(Xp,Ug)
% xlabel('Position')
% ylabel('Velocity')
% ylim([-1 1])