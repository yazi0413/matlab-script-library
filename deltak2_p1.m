%%%%%%%%%%%% deltak2
clear

ff=50;
load(['p1rsr_' num2str(ff)])
load(['p1rrr_' num2str(ff)])
load(['p1_' num2str(ff)])
load(['p0_' num2str(ff)])
deltac=0.1;c=(1500+deltac);
c0=1500;
w=2*pi*ff;
dk2=w^2/c0^2-w^2/c^2;
eval(['thesum=p1rsr_' num2str(ff) '.*p1rrr_' num2str(ff) '*dk2;'])


thesum=sum(thesum,1);
thesump0=sum(thesum,2)
eval(['dp=p1_' num2str(ff) '-p0_' num2str(ff) ';'])
dp(20)
% save thesump1 thesump1
