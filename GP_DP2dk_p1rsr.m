%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the ram's output files
%-- ���Ŷ����մ������ݶ�ȡ����
%---------p1(rs,rr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;close all;fclose all;

tic;   %��ʼ��ʱ
cd workdir
ff1='deltac.grid';

deltac=5;
R1=4000;R2=6000;  %---���뷴�ݼ�����Ŷ���ˮƽ��Χ
Ndz=10;  %-----��ȷ����Ŷ����ֵ�����ˮƽ����dr������

fname1=fopen(ff1,'rb');
nf=fread(fname1,1,'int');
f0=fread(fname1,1,'double');
df=fread(fname1,1,'double');
nsd=fread(fname1,1,'int');
sd=fread(fname1,nsd,'double');
rmax=fread(fname1,1,'double');
nz1=fread(fname1,1,'int');       %nz1=zmplt/dz �൱�������ȷ����ϵ�������
ndz=fread(fname1,1,'int');
ndr=fread(fname1,1,'int');
dz=fread(fname1,1,'double');
dr=fread(fname1,1,'double');
headersize = ftell(fname1);  %-------------�����ʱ��λ�ã�����GRID�ļ�ͷ����һЩ�������Ĵ�С
drr = dr * ndr;
nr = floor(rmax / drr);
r = [1:nr]*drr;
z = [0:nz1-1]*ndz*dz;
f = f0+ (0:nf-1)*df;
blocksize = 2*nsd*nz1*nr*8;
thesd=sd;                  %---------------��Դ��λ�ã�����ֻ��һ����Դ
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
    eval(['p1_' num2str(ff) '=p(2:5:end,end)/sqrt(2*pi/(rmax));']);
    
    
    for ir=R2:-dr:R1+dr;
        p1rsr(:,(ir-R1)/dr,ifreq)=p(2:Ndz:end,ir/dr)/sqrt(2*pi/(ir));    %!!!!!!!!!!!!!!!!!!!!!!!
    end
    %----------p1rsr(zr,r,frequency)
    
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
nz1=fread(fname1,1,'int');       %nz1 �����ȷ����ϵ�������
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
    
    for zr=2:5:99;
        % for zr=5+1;
        for ir=rmax-R1-dr:-dr:rmax-R2;
            p0rrr(:,(rmax-R1-ir)/dr,ifreq,(zr+3)/5)=p(zr,2:Ndz:end,ir/dr)/sqrt(2*pi/(ir)); %!!!!!!!!!!!!!!!!!!!!!!
        end
        %------------p0rrr(zs,r,frequency,zr)
    end
    
    eval(['p0_' num2str(ff) '=p(2:5:end,thesd+1,end)/sqrt(2*pi/(rmax));']);   %-----5mһ����������һ��Ƶ��20������
    
    
end
fclose(fname1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  save as A B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(Ndz*(R2-R1)/dr,Ndz*(R2-R1)/dr);B=zeros(Ndz*(R2-R1)/dr,1);

c0=1500;
for ff=f0:df:(nf-1)*df+f0;
    w=2*pi*ff;
    
    for zr=2:5:99;
        
        
        G=p1rsr(:,:,(ff-f0)/df+1).*p0rrr(:,:,(ff-f0)/df+1,(zr+3)/5)*dr*10/(2*sqrt(w/c0))*w^2;
        A((zr+3)/5+((ff-f0)/df)*20,:)=reshape(G,1,Ndz*(R2-R1)/dr);
        
    end
    
    eval(['dp=p1_' num2str(ff) '-p0_' num2str(ff) ';']);
    
    B(1+((ff-f0)/df)*20:((ff-f0)/df+1)*20,1)=dp;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c0=1500;c=1500+deltac;
% if R2==7000;
%     x=[zeros(100,1);(1/c0^2-1/c^2)*ones(100,1);zeros(100,1);];%the true answer
% elseif R2==6000;
%     x=[zeros(100,1);(1/c0^2-1/c^2)*ones(100,1)];
% else
%     x=[(1/c0^2-1/c^2)*ones(100,1)];
% end

% stem(abs(A*x./B))
format long;  %�޸�Ĭ����������
n=Ndz*(R2-R1)/dr;   %�����ģ
e=1e-4;  %����
gama=1e-13;
D=A'*A+gama*eye(n);
B=A'*B;
A=D;
x=[zeros(100,1);(1/c0^2-1/c^2)*ones(100,1)];
stem(abs(A*x./B))
X=zeros(Ndz*(R2-R1)/dr,1);  %��ʼ
r=B-A*X;
pp=r;
k=1;%��������
fan=norm(r);
% errcot=abs(mean(real(A*X./B)-1));
while fan>e
    zr=r'*r;
    afa=zr/(pp'*(A*pp));
    X=X+afa*pp;
    r=r-afa*(A*pp);
    fan=norm(r);
    % errcot=abs(mean(real(A*X./B)-1));
    beta=r'*r/zr;
    pp=r+beta*pp;
    k=k+1;
end

toc;
dc=1./sqrt(1/c0^2.-abs(X))-c0;
% Dc=1./sqrt(1/c0^2.-(x))-c0;
% figure;pcolor(dr:dr:2000,11:10:101,reshape(Dc,10,20))
% shading interp;xlabel(' r(m)','fontsize',15);ylabel('Depth','fontsize',15);set(gca,'ydir','reverse')
% title('�������');colorbar;
figure;pcolor(dr:dr:(R2-R1),11:10:101,reshape(dc,10,(R2-R1)/dr))
shading interp;xlabel(' r(m)','fontsize',15);ylabel('Depth','fontsize',15);set(gca,'ydir','reverse')
title('���ݽ��');colorbar;
% figure;
% stem(Dc);hold on;
% plot(dc,'r*-')
save dc dc
cd ..
