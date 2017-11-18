%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the ram's output files
%-- 有扰动接收处的数据读取出来
%---------p1(rs,rr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;close all;fclose all;
tic;   %开始计时
cd workdir
ff1='deltac.grid';

R1=5000;R2=6000;  %---代入反演计算的扰动区水平范围
Ndz=10;  %-----深度方向扰动区分的网格；水平方向按dr的来分


fname1=fopen(ff1,'rb');
nf=fread(fname1,1,'int');
f0=fread(fname1,1,'double');
df=fread(fname1,1,'double');
nsd=fread(fname1,1,'int');
sd=fread(fname1,nsd,'double');
rmax=fread(fname1,1,'double');
nz1=fread(fname1,1,'int');       %nz1=zmplt/dz 相当于输出深度方向上的网格数
ndz=fread(fname1,1,'int');
ndr=fread(fname1,1,'int');
dz=fread(fname1,1,'double');
dr=fread(fname1,1,'double');
headersize = ftell(fname1);  %-------------保存此时的位置，代表GRID文件头（存一些参数）的大小
drr = dr * ndr;
nr = floor(rmax / drr);
r = [1:nr]*drr;
z = [0:nz1-1]*ndz*dz;
f = f0+ (0:nf-1)*df;
blocksize = 2*nsd*nz1*nr*8;
thesd=sd;                  %---------------声源的位置，初期只用一个声源
for ifreq = 1:nf %nf
    %             loc = (isd-1)*blocksize*nf + (ifreq-1)*blocksize + headersize;
    loc = (ifreq-1)*blocksize + headersize;
    %loc = headersize;
    %             fseek(fname1,loc,-1);
    %             press=fread(fname1,[2 nz1*nr],'double');
    fseek(fname1,loc,-1);
    
    press=fread(fname1,[2 nsd*nz1*nr],'double');
    press=press(1,:)+1i*press(2,:);
    p=reshape(squeeze(press),nz1,nsd,nr);
    p=squeeze(p);
    ff=(ifreq-1)*df+f0;
    eval(['p1_' num2str(ff) '=p(6:5:end,end)/sqrt(2*pi/(rmax));']);
    
    %----------p0rsr(zr,r,frequency)
    for ir=R2:-dr:R1+dr;
        p0rsr(:,(ir-R1)/dr,ifreq)=p(Ndz+1:Ndz:end,ir/dr)/sqrt(2*pi/(ir));    %!!!!!!!!!!!!!!!!!!!!!!!
    end
end
fclose(fname1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------GP&p0(rs,rr)---------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff1='nodeltac.grid';
fname1=fopen(ff1,'rb');
nf=fread(fname1,1,'int');
f0=fread(fname1,1,'double') ;
df=fread(fname1,1,'double');
nsd=fread(fname1,1,'int');
sd=fread(fname1,nsd,'double');
rmax=fread(fname1,1,'double');
nz1=fread(fname1,1,'int');       %nz1 输出深度方向上的网格数
ndz=fread(fname1,1,'int');
ndr=fread(fname1,1,'int');
dz=fread(fname1,1,'double');
dr=fread(fname1,1,'double');
headersize = ftell(fname1);
drr = dr * ndr;
nr = floor(rmax / drr);
r = [1:nr]*drr;
z = [0:nz1-1]*ndz*dz;
f = f0+ (0:nf-1)*df;
blocksize = 2*nsd*nz1*nr*8;

for ifreq =1:nf; %nf
    %             loc = (isd-1)*blocksize*nf + (ifreq-1)*blocksize + headersize;
    loc = (ifreq-1)*blocksize + headersize;
    %loc = headersize;
    %             fseek(fname1,loc,-1);
    %             press=fread(fname1,[2 nz1*nr],'double');
    fseek(fname1,loc,-1);
    
    press=fread(fname1,[2 nsd*nz1*nr],'double');
    press=press(1,:)+1i*press(2,:);
    p=reshape(squeeze(press),nz1,nsd,nr);
    ff=(ifreq-1)*df+f0;
    %          %p(z,sd,r)
    
    
    
    for zr=6:5:101;
        % for zr=5+1;
        for ir=rmax-R1-dr:-dr:rmax-R2;
            p0rrr(:,(rmax-R1-ir)/dr,ifreq,(zr-1)/5)=p(zr,1+Ndz:Ndz:end,ir/dr)/sqrt(2*pi/(ir)); %!!!!!!!!!!!!!!!!!!!!!!
        end
        %------------p0rrr(zs,r,frequency,zr)
    end
    
    eval(['p0_' num2str(ff) '=p(6:5:end,thesd+1,end)/sqrt(2*pi/(rmax));']);   %-----5m一个接收器，一个频率20个方程
    
    
end
fclose(fname1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  save as A B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(Ndz*(R2-R1)/dr,Ndz*(R2-R1)/dr);B=zeros(Ndz*(R2-R1)/dr,1);

c0=1500;
for ff=f0:df:(nf-1)*df+f0;
    w=2*pi*ff;
    
    for zr=6:5:101;

        G=p0rsr(:,:,(ff-f0)/df+1).*p0rrr(:,:,(ff-f0)/df+1,(zr-1)/5)*dr*10/(2*sqrt(w/c0))*w^2;
        A((zr-1)/5+((ff-f0)/df)*20,:)=reshape(G,1,Ndz*(R2-R1)/dr);
        
    end
    
    eval(['dp=p1_' num2str(ff) '-p0_' num2str(ff) ';']);
    
    B(1+((ff-f0)/df)*20:((ff-f0)/df+1)*20,1)=dp;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load dc
c0=1500;
x=1/c0^2-1./(c0+dc).^2;
stem(abs(A*x./B))
%
format long;  %修改默认数据类型
n=Ndz*(R2-R1)/dr;   %矩阵规模
e=1e-1;  %精度
gama=1e-19;
D=A'*A+gama*eye(n);
B=A'*B;
A=D;

X=zeros(Ndz*(R2-R1)/dr,1);  %初始
r=B-A*X;
p=r;
k=1;%迭代次数
fan=norm(r);
% errcot=abs(mean(real(A*X./B)-1));
while fan>e
    zr=r'*r;
    afa=zr/(p'*(A*p));
    X=X+afa*p;
    r=r-afa*(A*p);
    fan=norm(r);
    % errcot=abs(mean(real(A*X./B)-1));
    beta=r'*r/zr;
    p=r+beta*p;
    k=k+1;
end

dc=1./sqrt(1/c0^2.-abs(X))-c0;
% dc=abs(X)*c0^3/2;
Dc=1./sqrt(1/c0^2.-(x))-c0;


figure;pcolor(dr:dr:R2-R1,11:10:101,reshape(dc,10,(R2-R1)/dr))
shading interp;xlabel(' r(m)','fontsize',15);ylabel('Depth','fontsize',15);set(gca,'ydir','reverse')
title('反演结果');colorbar;
% figure;stem(Dc);hold on;plot(dc,'r*-')
% dc=dc(101:end);
save dc dc;
cd ..

toc