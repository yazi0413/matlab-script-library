clear all; close all;
format long;  %修改默认数据类型 
tic;   %开始计时 
n=300;   %矩阵规模 
l=10; %预处理子矩阵的个数 
e=1e-10;  %精度 
load A;
load B;
gama=1e7;
D=A'*A+gama*eye(n);
B=A'*B;
A=D;

c0=1500;c=1500+1;

K=cond(A);  %条件数 
M=zeros(n); 
for i=1:1:l 
   for j=1:1:(n/l) 
       for k=1:1:(n/l) 
           line=j+(i-1)*(n/l); 
           low=k+(i-1)*(n/l); 
           M(line,low)=A(line,low); 
       end 
   end 
end 
x=[zeros(100,1);(1/c0^2-1/c^2)*ones(100,1);zeros(100,1)];%the true answer

b=B;   %生成精确解对应的b 
X=zeros(n,1);  %初始 
ek1=sqrt((X-x)'*A*(X-x)); %开始是不会知道精确解的！ 
sulv=0;  % 评价指标 
r=b-A*X; 
ni=inv(M); 
z=ni*r; 
p=z; 
k=1;%迭代次数 
fan=norm(r); 
while (fan>e) 
zr=r'*z; 
afa=zr/(p'*(A*p)); 
X=X+afa*p; 
ek2=sqrt((X-x)'*A*(X-x)); 
sulv=sulv+ek2/ek1; 
ek1=ek2; 
r=r-afa*(A*p); 
z=ni*r; 
fan=norm(r); 
beta=r'*z/zr; 
p=z+beta*p; 
k=k+1; 
end 
en=norm(X-x)  %求解的误差 
sulv=sulv/(k-1)  %下降速率的平均值 

toc; 
dc=1./sqrt(1/c0^2.-abs(X))-c0;

Dc=1./sqrt(1/c0^2.-abs(x))-c0;

stem(Dc);hold on;plot((dc),'r*-')