function bin2dec(bin)
%
%
%二进制小数转化为十进制
%bin=输入的二进制小数
%bin必须为一个由01组成的向量
b=length(bin)-1;d=0;
for k=1:b
    d=d+bin(k+1)*2^(-k);
end
if sign(bin(1))==0;
    dec=d;
else
    dec=-d;
end
disp('十进制等效值：');
disp(dec)