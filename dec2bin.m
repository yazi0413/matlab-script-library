function dec2bin(d,bytes)
%
% 程序将十进制小数值转化为二进制数值
%包括十进制小数，正负数均可
%d=input('请输入十进制小数：');
%bytes=input('输入所需要的字长：');
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
disp('相应的二进制数为：');
disp(bin)