%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the ram's output files
%-- 有扰动接收处的数据读取出来
%---------p1(rs,rr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;close all;fclose all;
% tic;   %开始计时
cd workdir
ff1='deltac.grid'; %--------此时的deltac.in作为背景声速剖面

R1=4000;R2=6000;  %---代入反演计算的扰动区水平范围
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
    
    %          %p(z,sd,r)
    for ir=R2:-dr:R1+dr;
        p0rsr(:,(ir-R1)/dr,ifreq)=p(Ndz+1:Ndz:end,ir/dr)/sqrt(2*pi/(ir));    %!!!!!!!!!!!!!!!!!!!!!!!
    end
    %----------p0rsr(zr,r,frequency)
    
    eval(['p11_' num2str(ff) '=p(6:5:end,end)/sqrt(2*pi/(rmax));']);  %-----5m一个接收器，一个频率20个方程
    
end
fclose(fname1);



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
    
    for zr=6:5:101;
        % for zr=5+1;
        for ir=rmax-R1-dr:-dr:rmax-R2;
            p0rrr(:,(rmax-R1-ir)/dr,ifreq,(zr-1)/5)=p(zr,1+Ndz:Ndz:end,ir/dr)/sqrt(2*pi/(ir)); %!!!!!!!!!!!!!!!!!!!!!!
        end
        %------------p0rrr(zs,r,frequency,zr)
    end
end
fclose(fname1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  save as A B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(Ndz*(R2-R1)/dr,Ndz*(R2-R1)/dr);B=zeros(Ndz*(R2-R1)/dr,1);

c0=1500;dp_fan=0;p1_fan=0;
for ff=f0:df:(nf-1)*df+f0;
    w=2*pi*ff;
    
    for zr=6:5:101;

        G=p0rsr(:,:,(ff-f0)/df+1).*p0rrr(:,:,(ff-f0)/df+1,(zr-1)/5)*dr*10/(2*sqrt(w/c0))*w^2;
        A((zr-1)/5+((ff-f0)/df)*20,:)=reshape(G,1,Ndz*(R2-R1)/dr);
        
    end
    
    load(['p1_' num2str(ff)])
    eval(['dp=p11_' num2str(ff) '-p1_' num2str(ff) ';']);
    B(1+((ff-f0)/df)*20:((ff-f0)/df+1)*20,1)=dp;
    %-------||dp||/||p1||
    dp_fan = dp_fan+sum(abs(dp).^2);
    eval(['p1_fan = p1_fan+sum(abs(p1_' num2str(ff) ').^2);'])
end

%-------------------由上一次反演出来的剖面作为背景剖面计算出来的声压与实测量的相对误差
% dp_fan=sum(abs(B).^2);
dp_p1=dp_fan/p1_fan


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0=1500;

format long;  %修改默认数据类型
n=Ndz*(R2-R1)/dr;   %矩阵规模
e=1e4;     %================================================================精度
gama=1e-10;
D=A'*A+gama*eye(n);
B=A'*B;
A=D;

X=zeros(Ndz*(R2-R1)/dr,1);  %初始
r=B-A*X;
pp=r;
k=1;%迭代次数
fan=norm(r);

while fan>e
    zr=r'*r;
    afa=zr/(pp'*(A*pp));
    X=X+afa*pp;
    r=r-afa*(A*pp);
    fan=norm(r);
    beta=r'*r/zr;
    pp=r+beta*pp;
    k=k+1;
end


load dc;
Dc=(1./sqrt(1/c0^2.-real(X))-c0);


figure;pcolor(dr:dr:R2-R1,11:10:101,reshape(dc,10,(R2-R1)/dr))
shading interp;xlabel(' r(m)','fontsize',15);ylabel('Depth','fontsize',15);set(gca,'ydir','reverse')
title('上一次反演结果');colorbar;

dc=dc+Dc;

% %------------修正
% for k=1:length(dc);
%     if(dc(k)>=5)
%         dc(k)=5;
%     elseif(dc(k)<=0)
%         dc(k)=0;
%     end
% end
figure;pcolor(dr:dr:R2-R1,11:10:101,reshape(dc,10,(R2-R1)/dr))
shading interp;xlabel(' r(m)','fontsize',15);ylabel('Depth','fontsize',15);set(gca,'ydir','reverse')
title('累加反演结果');colorbar;


figure;pcolor(dr:dr:R2-R1,11:10:101,reshape(Dc,10,(R2-R1)/dr))
shading interp;xlabel(' r(m)','fontsize',15);ylabel('Depth','fontsize',15);set(gca,'ydir','reverse')
title('本次迭代的反演结果');colorbar;



%------------反演结果与真实扰动的相对误差
x=[zeros(100,1);5*ones(100,1)];
errdk=sum((dc-x).^2)/sum(x.^2)

% save dc dc;
cd ..
% toc