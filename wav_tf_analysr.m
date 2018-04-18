clc
close all;
clear
filepath='D:\MATLAB\work\tf_analyse_sound_quality';
cd(filepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s1 and s2 have same sampling rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the s1&s2 and plot
[Stemp,fs]=audioread('London90_DAI N6_Roger_Mic_Mute_DAI_Mute_.wav');
    %%%%%%%%%%  resample  %%%%%%%%%
if fs>12600;
    s1=zeros(floor(length(Stemp)/floor(fs/12600)),1);
   for ii=1:length(s1);
       s1(ii)=Stemp(ii*floor(fs/12600));
   end
   fs=12600;
end
    
N1=length(s1);
t1=(0:N1-1)/fs;

figure(1);
subplot(211);
plot(t1,s1);%time
xlabel('time(s)');
f=fs*(0:N1-1)/N1;
f=f';
subplot(212);
FFT_amp1=abs(fft(s1))*2/N1;
plot(f(1:floor(floor(N1/2))),FFT_amp1(1:floor(floor(N1/2))));%frequency spectrum
axis([0 fs/2 0 0.01])
xlabel('frequency(Hz)');
saveas(gcf,'signal1.jpg');

[Stemp,fs]=audioread('London90_DAI N6_GNR FM02_Mic_Mute_DAI_Mute_.wav');

    %%%%%%%%%%  resample  %%%%%%%%%
if fs>12600;
    s2=zeros(floor(length(Stemp)/floor(fs/12600)),1);
   for ii=1:length(s2);
       s2(ii)=Stemp(ii*floor(fs/12600));
   end
   fs=12600;
end

N2=length(s2);
if N2<N1;
   s2=[s2; zeros(N1-N2,1)];
elseif N2>N1;
    s2=s2(1:N1);
end
N2=length(s2);
t2=(0:N2-1)/fs;
% s2=s1+0.01*sin(2*pi*2000*t');

figure(2);
subplot(211);
plot(t2,s2);%time
xlabel('time(s)');
f=fs*(0:N2-1)/N2;
f=f';
subplot(212);
FFT_amp2=abs(fft(s2*2/N2));
plot(f(1:floor(N2/2)),FFT_amp2(1:floor(N2/2)));%freqency spectrum
axis([0 fs/2 0 0.01])
xlabel('frequency(Hz)');
saveas(gcf,'signal2.jpg');


figure(3);
subplot(211);
plot(t2,s2);%time
xlabel('time(s)');
subplot(212);
plot(t1,s1);%time
xlabel('time(s)');
saveas(gcf,'s1 vs s2.jpg');

figure(4)
subplot(211);
plot(f(1:floor(N2/2)),FFT_amp1(1:floor(N2/2)));%frequency spectrum
axis([0 fs/2 0 0.01])
xlabel('frequency(Hz)');
subplot(212);
plot(f(1:floor(N1/2)),FFT_amp2(1:floor(N1/2)));%freqency spectrum
axis([0 fs/2 0 0.01])
xlabel('frequency(Hz)');
saveas(gcf,'FFT s1 vs s2.jpg');

%% time-frequency analyze
window=128;
noverlap=64;
nfft=256;
mindb=20;
[S1,F,T]=spectrogram(s1,window,noverlap,nfft,fs);
[S2]=spectrogram(s2,window,noverlap,nfft,fs);
P1=10*log10(abs(S1));
P2=10*log10(abs(S2));

figure(5);
pcolor(T,F,P1);shading interp;
xlabel('time(s)');ylabel('frequency(Hz)');colorbar;
caxis([-50 10])
saveas(gcf,'time-freq_s1.jpg');
figure(6);
pcolor(T,F,P2);shading interp;
xlabel('time(s)');ylabel('frequency(Hz)');colorbar;
caxis([-50 10])
saveas(gcf,'time-freq_s2.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% STFT coherent analyse
coh_coef=abs(min(min(corrcoef(P1,P2))));
% den=sqrt(norm(P1))*sqrt(norm(P2));

% S=(P1-P2);
S = abs(S1)-abs(S2);


figure(7);
pcolor(T,F,S);shading interp;
xlabel('time(s)');
ylabel('frequency(Hz)');colorbar;
title(['corrcoef: ' num2str(coh_coef)], 'fontsize',12,'FontWeight','bold');
% caxis([-10 10]);
saveas(gcf,'difference_corr.jpg');
close all;