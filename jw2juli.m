function [S]=jw2juli(L1,B1,L2,B2)
%由已知的经纬度转化为距离
%L1=input('请输入B点的经度值L1=');
%B1=input('请输入B点的纬度值B1=');
%L2=input('请输入C点的经度值L2=');
%B2=input('请输入C点的纬度值B2=');

a= 6378245.0;%地球长轴
b = 6356863.019;%地球短轴
format long,b;
e2=(a*a-b*b)/(b*b);%第二偏心率的平方
%公式系
A=sqrt(1+e2*cosd(B1)*cosd(B1)*cosd(B1)*cosd(B1));
B=sqrt(1+e2*cosd(B2)*cosd(B2));
C=sqrt(1+e2);
OB=B2-B1; %计算A(L1,B1)与中点（L2，B2）纬度差
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
S=a*C*O/(1000*(B*B));%两点间的距离（千米）

%fprintf('The distanse is %5.3f ', S);
%angle=h-H;%方位角
%fprintf('The angle is %5.3f ', angle)
