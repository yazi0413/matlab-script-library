function dec2bin(d,bytes)
%
% ����ʮ����С��ֵת��Ϊ��������ֵ
%����ʮ����С��������������
%d=input('������ʮ����С����');
%bytes=input('��������Ҫ���ֳ���');
dd=abs(d);
beq=[zeros(1,bytes)];
for k=1:bytes
    int=fix(2*dd);
    beq(k)=int;
    dd=2*dd-int;
end
if sign(d)==-1;
    bin=[1 beq];
else
    bin=[0 beq];
end 
disp('��Ӧ�Ķ�������Ϊ��');
disp(bin)