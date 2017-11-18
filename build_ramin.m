close all;clear;clc;

load svptest.mat;

fid=fopen('ramsurf_test.in','w');
fs = -1; nfreq = 101; freq0 = 150.0; dfreq = 100.0; 
nsd = 1;
sd = 50.0;
zr = 15.4;
rmax=10800.0;dr=2.0;ndr=10;
zmax=1000.0;dz=0.2;ndz=7;zmplt=95.0;
c0=1530.0;np=8;ns=1;rs=0.0;   
rsrf = 0.0;  zsrf = 0.0;
fprintf(fid,'%s\n',['ramsurf.in']);
fprintf(fid,'%s\n',['svptest.tl']);
fprintf(fid,'%s\n',['svptest.grid']);
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
fprintf(fid,'%2d  %2d\n',[-1 -1]);

nr = length(range); ndepth = length(z_ctd);
for i=1:nr
    for j=1:ndepth
        fprintf(fid,'%8.2f  %8.2f\n',[z_ctd(j) svp1(j,i)]);
    end
    fprintf(fid,'%2d  %2d\n',[-1 -1]);
    fprintf(fid,'%8.2f  %8.2f\n',[0.00   1600.00]);
  %  fprintf(fid,'%8.2f  %8.2f\n',[104.00   1650.00]);
    fprintf(fid,'%2d  %2d\n',[-1 -1]);
    fprintf(fid,'%8.2f  %8.2f\n',[0.00      1.80 ]);
   % fprintf(fid,'%8.2f  %8.2f\n',[104.00      1.80 ]);
    fprintf(fid,'%2d  %2d\n',[-1 -1]);
    fprintf(fid,'%8.2f  %8.2f\n',[0.00      0.15]);
    fprintf(fid,'%8.2f  %8.2f\n',[500.00      0.15]);
    fprintf(fid,'%8.2f  %8.2f\n',[500.00     10.00]);
    fprintf(fid,'%8.2f  %8.2f\n',[600.00     10.00]);
    fprintf(fid,'%2d  %2d\n',[-1 -1]);
     if(i~=nr)
        fprintf(fid,'%8.2f\n',range(i+1));
    end
end
fclose(fid);