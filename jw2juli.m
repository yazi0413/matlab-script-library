function [S]=jw2juli(L1,B1,L2,B2)
%����֪�ľ�γ��ת��Ϊ����
%L1=input('������B��ľ���ֵL1=');
%B1=input('������B���γ��ֵB1=');
%L2=input('������C��ľ���ֵL2=');
%B2=input('������C���γ��ֵB2=');

a= 6378245.0;%������
b = 6356863.019;%�������
format long,b;
e2=(a*a-b*b)/(b*b);%�ڶ�ƫ���ʵ�ƽ��
%��ʽϵ
A=sqrt(1+e2*cosd(B1)*cosd(B1)*cosd(B1)*cosd(B1));
B=sqrt(1+e2*cosd(B2)*cosd(B2));
C=sqrt(1+e2);
OB=B2-B1; %����A(L1,B1)���е㣨L2��B2��γ�Ȳ�
D=OB*(1+3*e2*OB*sind(2*B1+2*OB/3)/(4*B*B))/(2*B);
w=A*(L2-L1)/2;
E=sind(D)*cosd(w);
F=sind(w)*(B*cosd(B1)*cosd(D)-sind(B1)*sind(D))/A;
h=atand(F/E);
if (E<0)

     h=h+180;
      
end
O=2*asind(sqrt(E*E+F*F));
H=atand(tand(w)*(sind(B1)+B*cosd(B1)*tand(D))/A);
S=a*C*O/(1000*(B*B));%�����ľ��루ǧ�ף�

%fprintf('The distanse is %5.3f ', S);
%angle=h-H;%��λ��
%fprintf('The angle is %5.3f ', angle)
