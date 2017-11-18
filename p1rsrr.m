%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 声源不移动 in 文件 p1rsrr p0rsrr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
cd workdir
deltac=1;
for ff=60:10:200;
    %%%%%%%%%%%%%%%parameters
    H=102;
    deltaz=H/30;
    fs = -1; nfreq = 1; freq0 = ff; dfreq = ff;
    nsd = 1;
    sd = 20.0;
    zr = 20;
    rmax=10000.0;dr=100.0;ndr=1; 
    zmax=H+300;dz=1;ndz=1;zmplt=H;
    c0=1500;np=10;ns=1;rs=0.0;
    rsrf = 0.0;  zsrf = 0.0;
    R=20;
    dt=0;
    
    for time=dt;
        
        infilname=['temcoh' num2str(time) 'f' num2str(ff) '.in'];
        fid=fopen(infilname,'w');
        
        fprintf(fid,'%s\n',['temcoh' num2str(time) 'f' num2str(ff) '.in']);
        fprintf(fid,'%s\n',['temcoh' num2str(time) 'f' num2str(ff) '.tl']);
        fprintf(fid,'%s\n',['temcoh' num2str(time) 'f' num2str(ff) '.grid']);
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
            
            if(r==0)
                c(1:31)=1500*ones(31,1);
            elseif(r<5000)
                c(1:31)=1500*ones(31,1);
                fprintf(fid,'%8f\n',r);
            elseif(r<=5500)
                
                
                c=(1500+deltac)*ones(31,1);
%                 c=1500*ones(31,1);
                fprintf(fid,'%8f\n',r);
            else
                c(1:31)=1500*ones(31,1);
                fprintf(fid,'%8f\n',r);
            end
            
            
            for i=1:31;
                
                fprintf(fid,'%6.1f   %6.2f\n',deltaz*(i-1),c(i));
                
            end
            fprintf(fid,'%2d  %2d\n',[-1 -1]);
            fprintf(fid,'%8.2f  %8.2f\n',[0.00   1600.00]);
            fprintf(fid,'%2d  %2d\n',[-1 -1]);
            fprintf(fid,'%8.2f  %8.2f\n',[0.00      1.80 ]);
            fprintf(fid,'%2d  %2d\n',[-1 -1]);
            fprintf(fid,'%8.2f  %8.2f\n',[0.00      0.1]);
            fprintf(fid,'%8.2f  %8.2f\n',[zmax-10      0.1]);
            fprintf(fid,'%8.2f  %8.2f\n',[zmax-9     10.00]);
            fprintf(fid,'%8.2f  %8.2f\n',[zmax     10.00]);
            fprintf(fid,'%2d  %2d\n',[-1 -1]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fclose(fid);
    end
    
end
cd ..