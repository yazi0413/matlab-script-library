clear all; close all;
format long;  %�޸�Ĭ���������� 
tic;   %��ʼ��ʱ 
n=300;   %�����ģ 
l=10; %Ԥ�����Ӿ���ĸ��� 
e=1e-10;  %���� 
load A;
load B;
gama=1e7;
D=A'*A+gama*eye(n);
B=A'*B;
A=D;

c0=1500;c=1500+1;

K=cond(A);  %������ 
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

b=B;   %���ɾ�ȷ���Ӧ��b 
X=zeros(n,1);  %��ʼ 
ek1=sqrt((X-x)'*A*(X-x)); %��ʼ�ǲ���֪����ȷ��ģ� 
sulv=0;  % ����ָ�� 
r=b-A*X; 
ni=inv(M); 
z=ni*r; 
p=z; 
k=1;%�������� 
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
en=norm(X-x)  %������� 
sulv=sulv/(k-1)  %�½����ʵ�ƽ��ֵ 

toc; 
dc=1./sqrt(1/c0^2.-abs(X))-c0;

Dc=1./sqrt(1/c0^2.-abs(x))-c0;

stem(Dc);hold on;plot((dc),'r*-')