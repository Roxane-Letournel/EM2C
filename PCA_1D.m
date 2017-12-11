close all
clear all

Ns = 100;
figure
for k=1:Ns
    Xp0 = rand;
    Up0 = 0;
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
scatter(x,latent,40,'filled')
xlim([0 5])
% figure
% plot(Xp,Up)
% hold on
% plot(Xp,Ug)
% xlabel('Position')
% ylabel('Velocity')
% ylim([-1 1])