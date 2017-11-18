%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%做t1-t2的相关强度图
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;

deltaT=1;rmax=40000;ff=100;dr=1000;
Tau=560; %    delay time 范围，单位为分钟
%TT 起始时刻
Ft=12*360;
cd('过渡数据')
R=rmax;
%  ar=eye(561);
%%%生成各个时刻的变量，并平均
for t=0:6*deltaT:Ft;
    for k=t:1:t+5;
        name=['temcoh' num2str(k) 'f' num2str(ff)];
        load(name);
        eval(['p' num2str(k) '=pp(:,R/dr);'])
        
    end
    eval(['pp' num2str(t/6) '=mean([p' num2str(t) ' p' num2str(t+1) ' p' num2str(t+2) ' p' num2str(t+3) ' p'...
        num2str(t+4) ' p' num2str(t+5) '],2);'])
    clear(['p' num2str(t)],['p' num2str(t+1)],['p' num2str(t+2)],['p' num2str(t+3)],['p'...
        num2str(t+4)],['p' num2str(t+5)])
end
%%%%%%%%above 生成的ppnum对应为第num分钟内的声压平均值，以供下面的运算
i=1;
for TT=0:1:560;
    P0=0;
    for t=TT:10*deltaT:Ft/6;
        eval(['P0=P0+pp' num2str(t) '.^2;']) %<p(t1)p(t1)>
    end
    
    j=1;
    for tao=-TT:1:0; %tao 为负的情况
        P=0;Ptao=0;
        for t=TT:10*deltaT:Ft/6+tao;
            eval(['P=P+pp' num2str(t) '.*pp' num2str(t+tao) ';'])
            eval(['Ptao=Ptao+pp' num2str(t+tao) '.^2;'])
        end
        ar(i,j)=sum(P)/sqrt(sum(P0)*sum(Ptao));
        j=j+1;
    end
    
    for tao=1:Tau-TT;   %单位为分钟 tao为正的情况
        
        P=0;Ptao=0;
        for t=TT:10*deltaT:(Ft/6-tao);
            eval(['P=P+pp' num2str(t) '.*pp' num2str(t+tao) ';'])
            eval(['Ptao=Ptao+pp' num2str(t+tao) '.^2;'])
            
        end
        ar(i,j)=sum(P)/sqrt(sum(P0)*sum(Ptao));
        j=j+1;
    end
    i=i+1;
end

%%plot the curve

t1=linspace(0,Tau,Tau+1);
t2=t1;
pcolor(t1,t2,ar);shading interp;

cd ..