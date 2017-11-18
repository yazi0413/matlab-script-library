%%%%%% EOF

clear;
clc;fclose all;
%------------------ 自己生成的
nscatterz = 2;
%depth=input('input the depth (m) ');
depth=100;
%
% deltac=zeros(depth/nscatterz,3600);
% for t=1:3600;
%     deltac(:,t)=5*randn(depth/nscatterz,1);
% end
%
% R=deltac*deltac';
% [u,s,v]=svd(R);
% fid=fopen('eof.txt','w');
% fprintf(fid,'%6.2f',v);
% fclose(fid);
%
% cw=[1520*ones(1,21) 1520:-1:1500 1500*ones(1,80-42)];
% cw=cw'+deltac;
% plot(cw(:,1),1:depth);
% set(gca,'ydir','reverse');

%------------------ 借用历史资料
cd data;
load sp1_30.mat;
cm=mean(svp1_30,1);c=zeros(10800,201);
for ii=1:10800;
    c(ii,:)=svp1_30(ii,:)-cm(1,:);
end
R=c'*c;
[u,s,v]=svd(R);
cd ..
fid=fopen('eof.txt','w');
fprintf(fid,'%6.2f',v(1:4:depth*2,1:4:depth*2));
for ii=1:2:20;
    figure;
    plot(v(1:4:depth*2,ii),1:2:depth);set(gca,'ydir','reverse');
end
pause;
close all;
fclose(fid);



