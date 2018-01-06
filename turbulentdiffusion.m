%clear


for lll=6%6:8
    
    
    for kkk=2
    
    
        scheme=kkk; %0=upwind, 1=LW, 2=RK4-DF4, 3=Upwind second order

        %gas field
        L=3;
         %load Gas_filt50.mat
        load Ugfull.dat
        load Vgfull.dat
         Nxin=length(Ugfull);
        %Choose mesh
        % Nx=2^k+1
        Nx=2^lll+1;
        ratio=(Nxin-1)/(Nx-1);
        Ug=zeros(Nx,Nx);
        Vg=zeros(Nx,Nx);




        %load x.dat
        %load y.dat

        Ug(1:Nx,1:Nx)=Ugfull(1:ratio:Nxin,1:ratio:Nxin);
        Vg(1:Nx,1:Nx)=Vgfull(1:ratio:Nxin,1:ratio:Nxin);
        % 
        % Ug=Ug';
        % Vg=Vg';
        %Ug(:,:)=1;
        %Vg(:,:)=0;
        % 
         dx=L/(Nx-1);
        % 
         for i=1:Nx
             for j=1:Nx
                 x(i,j)=(i-1)*dx;
                 y(i,j)=(j-1)*dx;
             end
         end


        %Ug(:,:)=cos(2*pi/3*x).*sin(2*pi/3*y);
        %Vg(:,:)=-sin(2*pi/3*x).*cos(2*pi/3*y);

        T=150;



        Umax=max(max([abs(Ug); abs(Vg)]));


        CFL=0.4;
        Fo=0.2;

        %molecular diffusivity
        D=1.d-3;
        dt=min([Fo*dx^2/D,CFL*dx/Umax]);
        ReCell=Umax*dx/D;


        Nt=floor(T/dt);


        Ntest=1;
        for l=1:Ntest
           l
        Ynew=zeros(Nx,Nx);
        Yold=zeros(Nx,Nx);
        if scheme==2
        Yint1=zeros(Nx,Nx);
        Yint2=zeros(Nx,Nx);
        end
        posinitx=1.5;%rand(1,1)*3;
        posinity=1.5;%rand(1,1)*3;
        r=sqrt((x-posinitx).^2+(y-posinity).^2);
        rlim=0.5;
        rsave(1,l)=rlim;
        for i=1:Nx
            for j=1:Nx
                if (r(i,j)<rlim);
                    Yold(i,j)=1.0d0;
                end
            end
        end
        mass(1)=sum(sum(Yold))/Nx^2;
        RMS(1)=std(std(Yold))/mean(mean(Yold));
        Ynew=Yold;
        %surfc(Yold),shading interp,view(2)

        t(1)=0;

        veci=[1:Nx];
        vecip1=[2:Nx, 2];
        vecim1=[Nx-1, 1:Nx-1];



        tic
        for k=1:Nt
        %     for i=1:Nx
        %         for j=1:Nx
        %             ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
        %             jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
        %             im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
        %             jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
        %             Ynew(i,j)=Yold(i,j) - max(Ug(i,j),0.d0)*dt/dx*(Yold(i,j)-Yold(im1,j)) ...
        %                                 - min(Ug(i,j),0.d0)*dt/dx*(Yold(ip1,j)-Yold(i,j)) ...
        %                                 - max(Vg(i,j),0.d0)*dt/dx*(Yold(i,j)-Yold(i,jm1)) ...
        %                                 - min(Vg(i,j),0.d0)*dt/dx*(Yold(i,jp1)-Yold(i,j)) ;
        %         end
        %     end


            if scheme==0 
                for i=1:Nx
                    for j=1:Nx
                        ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
                        jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
                        im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
                        jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
                        %UPWIND + CENTERED
                       % Ynew(i,j)=Yold(i,j) - (Ug(i,j)>0.d0)*dt/dx*(Ug(i,j)*Yold(i,j)-Ug(im1,j)*Yold(im1,j)) ...
                       %                     - (Ug(i,j)<0.d0)*dt/dx*(Ug(ip1,j)*Yold(ip1,j)-Ug(i,j)*Yold(i,j)) ...
                       %                     - (Vg(i,j)>0.d0)*dt/dx*(Vg(i,j)*Yold(i,j)-Vg(i,jm1)*Yold(i,jm1)) ...
                       %                     - (Vg(i,j)<0.d0)*dt/dx*(Vg(i,jp1)*Yold(i,jp1)-Vg(i,j)*Yold(i,j)) ...
                       %                     + D*dt/dx^2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ...
                       %                                   + Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1)   );
                        %UPWIND + CENTERED
                        Ynew(i,j)=Yold(i,j) - max(Ug(i,j),0.d0)*dt/dx*(Yold(i,j)-Yold(im1,j)) ...
                                            - min(Ug(i,j),0.d0)*dt/dx*(Yold(ip1,j)-Yold(i,j)) ...
                                            - max(Vg(i,j),0.d0)*dt/dx*(Yold(i,j)-Yold(i,jm1)) ...
                                            - min(Vg(i,j),0.d0)*dt/dx*(Yold(i,jp1)-Yold(i,j)) ...
                                            + D*dt/dx^2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ...
                                                          + Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1)   );
                    end
                end

                %vectorized...

        %         Ynew(veci,veci)=Yold(veci,1:Nx)-max(Ug(veci,1:Nx),0.d0)*dt/dx*(Yold(veci,1:Nx)-Yold(vecim1,1:Nx))...
        %                                        -min(Ug(veci,1:Nx),0.d0)*dt/dx*(Yold(vecip1,1:Nx)-Yold(veci,1:Nx))...
        %                                        -max(Vg(1:Nx,veci),0.d0)*dt/dx*(Yold(1:Nx,veci)-Yold(1:Nx,vecim1))...
        %                                        -min(Vg(1:Nx,veci),0.d0)*dt/dx*(Yold(1:Nx,vecip1)-Yold(1:Nx,veci))...                  
        %                                        + D*dt/dx^2*( Yold(vecip1,veci)-2*Yold(veci,veci)+Yold(vecim1,veci) ...
        %                                                     + Yold(veci,vecip1)-2*Yold(veci,vecj)+Yold(veci,vecim1)   );
        %         
            elseif scheme==1        
                for i=1:Nx
                    for j=1:Nx
                        ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
                        jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
                        im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
                        jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
                        %LAX WENDROFF + CENTERED
        %                 Ynew(i,j) = Yold(i,j) - dt/2/dx*(Ug(ip1,j)*Yold(ip1,j)-Ug(im1,j)*Yold(im1,j)) ...
        %                                             - dt/2/dx*(Vg(i,jp1)*Yold(i,jp1)-Vg(i,jm1)*Yold(i,jm1)) ...
        %                                             + dt^2/2/dx^2*( 0.5*(Ug(ip1,j)+Ug(i,j))*(Ug(ip1,j)*Yold(ip1,j)-Ug(i,j)*Yold(i,j)) ...
        %                                                            -0.5*(Ug(i,j)+Ug(im1,j))*(Ug(i,j)*Yold(i,j)-Ug(im1,j)*Yold(im1,j)) ) ...                 );
        %                                             + dt^2/2/dx^2*( 0.5*(Vg(i,jp1)+Vg(i,j))*(Vg(i,jp1)*Yold(i,jp1)-Vg(i,j)*Yold(i,j)) ...
        %                                                            -0.5*(Vg(i,j)+Vg(i,jm1))*(Vg(i,j)*Yold(i,j)-Vg(i,jm1)*Yold(i,jm1)) )...
        %                                             + dt^2/8/dx^2* ( Ug(ip1,j)*(Vg(ip1,jp1)*Yold(ip1,jp1)-Vg(im1,jp1)*Yold(im1,jp1))  ...
        %                                                             -Ug(im1,j)*(Vg(ip1,jm1)*Yold(ip1,jm1)-Vg(im1,jm1)*Yold(im1,jm1)) ...
        %                                                             +Vg(i,jp1)*(Ug(ip1,jp1)*Yold(ip1,jp1)-Ug(ip1,jm1)*Yold(ip1,jm1))  ...
        %                                                             -Vg(i,jm1)*(Ug(im1,jp1)*Yold(im1,jp1)-Ug(im1,jm1)*Yold(im1,jm1)) )... 
        %                                             +D*dt/dx^2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ...
        %                                                         + Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1)   );  
        %                LAX WENDROFF + CENTERED
                        Ynew(i,j) = Yold(i,j)       - Ug(i,j)*dt/2/dx*(Yold(ip1,j)-Yold(im1,j)) ...
                                                    - Vg(i,j)*dt/2/dx*(Yold(i,jp1)-Yold(i,jm1)) ...
                                                    + Ug(i,j)^2*dt^2/dx^2/2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ) ...
                                                    + Vg(i,j)^2*dt^2/dx^2/2*( Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1) )...
                                                    + D*dt/dx^2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ...
                                                                + Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1)   ); 
        %                 Ynew(i,j) = Yold(i,j)       - Ug(i,j)*dt/2/dx*(Yold(ip1,j)-Yold(im1,j)) ...
        %                                             - Vg(i,j)*dt/2/dx*(Yold(i,jp1)-Yold(i,jm1)) ...
        %                                             + Umax^2*dt^2/dx^2/2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ) ...
        %                                             + Umax^2*dt^2/dx^2/2*( Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1) )...
        %                                             + D*dt/dx^2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ...
        %                                                         + Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1)   ); 
                    end
                end     
            elseif scheme==2
                %%RK4-DF4

                %%step 1 
                for i=1:Nx
                    for j=1:Nx
                        ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
                        jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
                        im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
                        jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
                        ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
                        jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
                        im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
                        jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


                        dYdx=primederivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(im1,j),Yold(im2,j),dx);
                        dYdy=primederivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,jm1),Yold(i,jm2),dx);
                        d2Ydx2=secondderivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(i,j),Yold(im1,j),Yold(im2,j),dx);
                        d2Ydy2=secondderivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,j),Yold(i,jm1),Yold(i,jm2),dx);

                        %LAX WENDROFF + CENTERED
                        Yint1(i,j) = Yold(i,j)   +0.5d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                                    + D*dt*(d2Ydx2+d2Ydy2 )  );      
                    end
                end  
                %%step 2
                for i=1:Nx
                    for j=1:Nx
                        ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
                        jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
                        im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
                        jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
                        ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
                        jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
                        im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
                        jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


                        dYdx=primederivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(im1,j),Yint1(im2,j),dx);
                        dYdy=primederivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,jm1),Yint1(i,jm2),dx);
                        d2Ydx2=secondderivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(i,j),Yint1(im1,j),Yint1(im2,j),dx);
                        d2Ydy2=secondderivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,j),Yint1(i,jm1),Yint1(i,jm2),dx);

                        %LAX WENDROFF + CENTERED
                        Yint2(i,j) = Yold(i,j)   +0.5d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                                    + D*dt*(d2Ydx2+d2Ydy2 )  );     
                    end
                end 
                %%step 3
                for i=1:Nx
                    for j=1:Nx
                        ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
                        jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
                        im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
                        jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
                        ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
                        jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
                        im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
                        jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


                        dYdx=primederivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(im1,j),Yint2(im2,j),dx);
                        dYdy=primederivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,jm1),Yint2(i,jm2),dx);
                        d2Ydx2=secondderivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(i,j),Yint2(im1,j),Yint2(im2,j),dx);
                        d2Ydy2=secondderivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,j),Yint2(i,jm1),Yint2(i,jm2),dx);

                        %LAX WENDROFF + CENTERED
                        Yint3(i,j) = Yold(i,j)   +(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                                    + D*dt*(d2Ydx2+d2Ydy2 )  );     
                    end
                end 
                %%step 4 (final)
                for i=1:Nx
                    for j=1:Nx
                        ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
                        jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
                        im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
                        jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
                        ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
                        jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
                        im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
                        jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);


                        dYdx=primederivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(im1,j),Yold(im2,j),dx);
                        dYdy=primederivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,jm1),Yold(i,jm2),dx);
                        d2Ydx2=secondderivative_4th(Yold(ip2,j),Yold(ip1,j),Yold(i,j),Yold(im1,j),Yold(im2,j),dx);
                        d2Ydy2=secondderivative_4th(Yold(i,jp2),Yold(i,jp1),Yold(i,j),Yold(i,jm1),Yold(i,jm2),dx);

                        %LAX WENDROFF + CENTERED
                        Ynew(i,j) = Yold(i,j)   +1.d0/6.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                                    + D*dt*(d2Ydx2+d2Ydy2 )  ); 


                        dYdx=primederivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(im1,j),Yint1(im2,j),dx);
                        dYdy=primederivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,jm1),Yint1(i,jm2),dx);
                        d2Ydx2=secondderivative_4th(Yint1(ip2,j),Yint1(ip1,j),Yint1(i,j),Yint1(im1,j),Yint1(im2,j),dx);
                        d2Ydy2=secondderivative_4th(Yint1(i,jp2),Yint1(i,jp1),Yint1(i,j),Yint1(i,jm1),Yint1(i,jm2),dx);


                        Ynew(i,j) = Ynew(i,j)   +1.d0/3.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                                    + D*dt*(d2Ydx2+d2Ydy2 )  ); 

                        dYdx=primederivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(im1,j),Yint2(im2,j),dx);
                        dYdy=primederivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,jm1),Yint2(i,jm2),dx);
                        d2Ydx2=secondderivative_4th(Yint2(ip2,j),Yint2(ip1,j),Yint2(i,j),Yint2(im1,j),Yint2(im2,j),dx);
                        d2Ydy2=secondderivative_4th(Yint2(i,jp2),Yint2(i,jp1),Yint2(i,j),Yint2(i,jm1),Yint2(i,jm2),dx);

                        Ynew(i,j) = Ynew(i,j)   +1.d0/3.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                                    + D*dt*(d2Ydx2+d2Ydy2 )  );

                        dYdx=primederivative_4th(Yint3(ip2,j),Yint3(ip1,j),Yint3(im1,j),Yint3(im2,j),dx);
                        dYdy=primederivative_4th(Yint3(i,jp2),Yint3(i,jp1),Yint3(i,jm1),Yint3(i,jm2),dx);
                        d2Ydx2=secondderivative_4th(Yint3(ip2,j),Yint3(ip1,j),Yint3(i,j),Yint3(im1,j),Yint3(im2,j),dx);
                        d2Ydy2=secondderivative_4th(Yint3(i,jp2),Yint3(i,jp1),Yint3(i,j),Yint3(i,jm1),Yint3(i,jm2),dx);

                        Ynew(i,j) = Ynew(i,j)   +1.d0/6.d0*(  - Ug(i,j)*dt*dYdx- Vg(i,j)*dt*dYdy ...
                                                    + D*dt*(d2Ydx2+d2Ydy2 )  );                               
                    end
                end 
            elseif scheme==3
                %%Upwind second order

                %%step 1 
                for i=1:Nx
                    for j=1:Nx
                        ip1=  (i==Nx)*(2) +(i~=Nx)*(i+1);
                        jp1=  (j==Nx)*(2) +(j~=Nx)*(j+1);
                        im1=  (i==1)*(Nx-1) +(i~=1)*(i-1);
                        jm1=  (j==1)*(Nx-1) +(j~=1)*(j-1);
                        ip2=  (i==Nx-1)*(2) +(i==Nx)*(3) +((i~=Nx)&(i~=Nx-1))*(i+2);
                        jp2=  (j==Nx-1)*(2) +(j==Nx)*(3) +((j~=Nx)&(j~=Nx-1))*(j+2);
                        im2=  (i==1)*(Nx-2) + (i==2)*(Nx-1) +((i~=1)&(i~=2))*(i-2);
                        jm2=  (j==1)*(Nx-2) + (j==2)*(Nx-1) +((j~=1)&(j~=2))*(j-2);



                        Ynew(i,j) = Yold(i,j)       - max(Ug(i,j),0.d0)*dt/2/dx*(3*Yold(i,j)-4*Yold(im1,j)+Yold(im2,j)) ...
                                                    - min(Ug(i,j),0.d0)*dt/2/dx*(-3*Yold(i,j)+4*Yold(ip1,j)-Yold(ip2,j)) ...
                                                    - max(Vg(i,j),0.d0)*dt/2/dx*(3*Yold(i,j)-4*Yold(i,jm1)+Yold(i,jm2)) ...
                                                    - min(Vg(i,j),0.d0)*dt/2/dx*(-3*Yold(i,j)+4*Yold(i,jp1)-Yold(i,jp2)) ...
                                                    + D*dt/dx^2*( Yold(ip1,j)-2*Yold(i,j)+Yold(im1,j) ...
                                                                + Yold(i,jp1)-2*Yold(i,j)+Yold(i,jm1)   );  
                    end
                end 
            end

            t(k+1)=t(k)+dt;
            t(k+1)
            Yold=Ynew;
            mass(k+1)=sum(sum(Ynew))/Nx^2;
            %RMS(k+1)=std(std(Ynew))/mean(mean(Ynew));
            RMS(k+1)=std(Ynew(:))/mean(mean(Ynew));
            if RMS(k+1)<0.05
                fprintf('mixing reached \n')
                fprintf('time=%f s\n',t(k+1))
                break
            end
            %gaussian fit to get variance
            %optimisation
            rtest=0;
            test=1;

        %     while test
        %         Ysum=0;
        %         for i=1:Nx
        %             for j=1:Nx
        %                 if (r(i,j)<rtest);
        %                     Ysum=Ysum+Ynew(i,j)/Nx^2;
        %                 end
        %             end
        %         end
        %         if (Ysum>=0.95*mass(k+1) ) 
        %             test=0;
        %             rsave(k+1,l)=rtest;
        %         else
        %             rtest=rtest+0.01*rlim;
        %             if rtest>L
        %                 test=0;
        %             end
        %         end
        %     end

        end
        toc

        figure(1);subplot(1,2,1),plot(t,mass),hold on
        figure(1);subplot(1,2,2),plot(t,RMS),hold on
        %figure;surfc(x,y,Ynew),shading interp,view(2),colormap('jet')

        end
    end
end


