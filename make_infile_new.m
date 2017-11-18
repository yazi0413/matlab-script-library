%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make the .in files for ram

% simulate the internal wave
%多频的情况
%改进后的程序，效率有一定提高
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;clear;clc;
cd('过渡数据')
for ff=100:10:100;
     %%%%%%%%%%%%%%%%    ssp
        %add pertubation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H=3000;R=2;
        B=1000;M=(pi*3-1)/(2*9);N0=5.2e-3;E=6.3e-5;fi=7.3e-5;
        q=pi*fi/(B*N0);wi=2*pi*fi;deltaz=100;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %The perturbation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  还没处理好边界条件。
            kk=1;deltak=6.28e-6;
            for k=6.28e-6:deltak:6.28e-5;
                C(1,:)=[(2+1e4*k^2)/(1e4*(N0^2*exp(2*deltaz/B)-wi^2)) 1/(1e4*sqrt((N0^2*...
                    exp(2*2*deltaz/B)-wi^2)*(N0^2*exp(2*deltaz/B)-wi^2))) zeros(1,29)];
                
                b(1)=sqrt(N0^2*exp(2*deltaz/B)-wi^2)*deltaz;   
                for i=2:30;
                    C(i,:)=[zeros(1,i-2) 1/(1e4*sqrt((N0^2*exp(2*i*deltaz/B)-wi^2)*...
                        (N0^2*exp(2*(i-1)*deltaz/B)-wi^2))) (2+1e4*k^2)/(1e4*(N0^2*exp(2*i*deltaz/B)-wi^2)) (1e4*sqrt((N0^2*exp(2*i*deltaz/B)-wi^2)*...
                        (N0^2*exp(2*(i+1)*deltaz/B)-wi^2))) zeros(1,30-i)];
                    
                    b(i)=sqrt(N0^2*exp(2*i*deltaz/B)-wi^2)*deltaz;
                end
                
                C(31,:)=[zeros(1,29) 1/(1e4*sqrt((N0^2*exp(2*30*deltaz/B)-wi^2)*(N0^2*exp(2*29*deltaz/B)-wi^2))) ...
                    (2+1e4*k^2)/(1e4*(N0^2*exp(2*30*deltaz/B)-wi^2))];
                
                b(31)=sqrt(N0^2*exp(2*30*deltaz/B)-wi^2)*deltaz;  
                [V,lamda]=eig(C);
                lamda=eig(C);
                W=V/b;g=0;
                for j=1:50;
                    g=g+sqrt(2*B^3*N0^2*E*j*q*k/(pi^2*M*(j^2+9)*(k^2+j^2*q^2)^2));  %%%存在问题，应该具体随机性，每次运行的结果是不一样的。
                end
                syms rr tt;
                for i=1:31;
                    rt(i)=cos(k*rr-sqrt(k^2/abs(lamda(i))+wi^2)*tt);
                    kesi(i,kk)=g*W(i)*rt(i);
                end
                kk=kk+1;
            end
            Kesi=sum(kesi,2);
            
            
    deltaT=10;Ft=600;  %%%%%是以秒为单位吗？这是一个问题
    for t=0:deltaT:Ft;  
        infilname=['temcoh' num2str(t) 'f' num2str(ff) '.in'];
        fid=fopen(infilname,'w');
        fs = -1; nfreq = 1; freq0 = ff; dfreq = ff; 
        nsd = 1;
        sd = 807.0;
        zr = 890;
        rmax=1200000.0;dr=1000.0;ndr=100; %试一下：dr取１公里的步长；grid二维网格取100公里的步长，取到1200公里的范围
        zmax=3500.0;dz=10;ndz=1;zmplt=3000.0;
        c0=1500;np=10;ns=1;rs=0.0;   
        rsrf = 0.0;  zsrf = 0.0;
        fprintf(fid,'%s\n',['temcoh' num2str(t) 'f' num2str(ff) '.in']);
        fprintf(fid,'%s\n',['temcoh' num2str(t) 'f' num2str(ff) '.tl']);
        fprintf(fid,'%s\n',['temcoh' num2str(t) 'f' num2str(ff) '.grid']);
        fprintf(fid,'%2d %8.2f %8.2f %8.2f\n',[fs,nfreq,freq0,dfreq]);
        fprintf(fid,'%2d\n',[nsd]);
        fprintf(fid,'%8.2f\n',[sd]);
        fprintf(fid,'%8.2f\n',[zr]);
        fprintf(fid,'%8.2f %6.1f %2d\n',[rmax, dr, ndr]);
        fprintf(fid,'%8.2f %8.2f %2d %8.2f\n',[zmax,dz, ndz, zmplt]);
        fprintf(fid,'%8.2f %2d %2d %8.2f\n',[ c0, np, ns, rs]);
        fprintf(fid,'%8.2f %8.2f\n',[rsrf zsrf]);
        fprintf(fid,'%2d  %2d\n',[-1 -1]);
        fprintf(fid,'%8.2f  %8.2f\n',[0,zmplt]);
        fprintf(fid,'%8.2f  %8.2f\n',[rmax,zmplt]);
        fprintf(fid,'%2d  %2d\n',[-1 -1]);
        
       
        for r=0:rmax/R:rmax;
            if(r>0)
                fprintf(fid,'%8f\n',r);
            end
            
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Above the perturbation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            i=1;
            for z=0:100:H;
                %             zz=2*(z-1300)/1300;
                %             c0(i)=1500*(1+0.00737*(zz-1+exp(-zz)));
                %             Nz=N0*exp(-z/H);
                %             deltac=c0(i)*(24.5/9.8*Nz^2);
                %             kesi=(1+0.5*t*abs(a(r+1)));
                %             c(i)=c0(i)+deltac*kesi;
                Nz=N0*exp(z/B);
                
                c0(i)=1500*(1+0.0057*(exp(-2*(z-890)/B)+2*(z-890)/B-1));
                %kesi=cos(k*r-w*t)*sqrt(G)*cos(k*sqrt(N0^2/w^2-1)*z);
                
                deltac(i)=c0(i)*24.5*Nz^2*subs(Kesi(i),[rr,tt],[r,t])/9.8;
                c(i)=c0(i)+deltac(i);
                
                fprintf(fid,'%6.1f   %6.2f\n',z,c(i));
                i=i+1;
            end
            fprintf(fid,'%2d  %2d\n',[-1 -1]);
            fprintf(fid,'%8.2f  %8.2f\n',[0.00   1600.00]);
            fprintf(fid,'%2d  %2d\n',[-1 -1]);  
            fprintf(fid,'%8.2f  %8.2f\n',[0.00      1.80 ]);           
            fprintf(fid,'%2d  %2d\n',[-1 -1]);
            fprintf(fid,'%8.2f  %8.2f\n',[0.00      0.8]);
            fprintf(fid,'%8.2f  %8.2f\n',[3100.00      0.8]);
            fprintf(fid,'%8.2f  %8.2f\n',[3101.00     10.00]);
            fprintf(fid,'%8.2f  %8.2f\n',[3500.00     10.00]);
            fprintf(fid,'%2d  %2d\n',[-1 -1]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fclose(fid);
    end
    
end
cd ..