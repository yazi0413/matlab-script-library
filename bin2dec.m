function bin2dec(bin)
%
%
%������С��ת��Ϊʮ����
%bin=����Ķ�����С��
%bin����Ϊһ����01��ɵ�����
b=length(bin)-1;d=0;
for k=1:b
    d=d+bin(k+1)*2^(-k);
end
if sign(bin(1))==0;
    dec=d;
else
    dec=-d;
end
disp('ʮ���Ƶ�Чֵ��');
disp(dec)