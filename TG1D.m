function [Xp,Up,Ug]=TG1D(Xp0,Up0,taup)

%PARAMETRES
% taup=20;

deltat=0.01;
Tend=1000;
Npas=Tend/deltat;

%Initialisation
% Xp0=0;
% Up0=1;


%copie
Xp(1)=Xp0;
Up(1)=Up0;

for i=1:Npas-1

    Xp(i+1)=Xp(i)+deltat*Up(i);
    Ug(i)=sin(Xp(i));
    
    Up(i+1)=Up(i)+deltat*(Ug(i)-Up(i))/taup;
    
end
Ug(Npas)=sin(Xp(Npas));
return
end
