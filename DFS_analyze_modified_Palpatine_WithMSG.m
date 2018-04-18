close all;
clear;
clc;

mirror = -1; % set to -1 if gain curve is desired instead of attenuation curve
plot_calibrated_response = 1; % set to 1 if calibrations should be taken into account
sub_name = '\DFSout\'; % set subfolder name
production_calibration = 0; % output cal offset, just copy the number here

%% System cal in Palpatine
%

% palpatine2 Delhi80
cal =[-23 -23 -21 -21 -25 -27 -27 -27 -21 -24 -18 -13 -14 -17 -11 -5 -5];

% No system cal at all
%cal = [0	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

cal = cal + production_calibration -6;
% fre = [0 168 340 518 706 908 1131 1383 1674 2019 2442 2981 3693 4672 6059 7997 10417]; % High bandwitdth frequencies
fre = [0 163 329 501 682 876 1090 1328 1600 1919 2301 2773 3370 4140 5136 6383 7812]; % High bandwitdth frequencies

% frq_FOG=[250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000,8000]; %Dooku
frq_FOG=[250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000]; %palpatine

FOG=[57 59 59 61 62 58 56 51 51];  %Delhi80
%FOG=[51, 51, 51, 51, 51, 48, 48, 48, 48, 48];  %Luxor NP
%FOG=[50, 50, 50, 50, 54, 56, 53, 42, 42, 42]; %HP
%FOG=[43, 43, 43, 42, 43, 43, 42, 41, 41, 41]; %LP
color_array = {'-b','-r','-c','-m','-g','-y','-k','--b','--r','--c','--m','--g','--y','--k'};

% 
% %cal = [6	6	9	10	11	11	11	11	11	11	11	15	17	16	15	12	7]; %% System cal in Dooku1 BerlinRIE LP
% 
% %cal = [-3	-3 1 4 5 6 6 7 7 8 9 11 17 17 13 8 8]; %% System cal in Dooku1 BerlinRIE NP
% 
% %cal = [6	6	9	10	11	11	11	11	11	11	11	15	17	16	15	12	7];
% %cal = [-12 	-11	-10	-7	-6	-3	-4	-4	-4	-3	-2	0	7	8	6 5	5]; %p6 luxor NP System cal
% cal = [-0 	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0];
% cal = cal + production_calibration - 6;
% %fre = [0 168 340 518 706 908 1131 1383 1674 2019 2442 2981 3693 4672 6059 7997 10417]; % system cal frq at Dooku
% fre = [0 163 329 501 682 876 1090 1328 1600 1919 2301 2773 3370 4140 5136 6383 7812]; % P6
% frq_FOG=[250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000,8000]; %Dooku
% FOG=[51, 51, 51, 51, 51, 48, 48, 48, 48, 48];  %Luxor NP


%% plot more OMNI curves together
currentFolder = pwd;
currentFolder = strcat(currentFolder,'\');
string_path = strcat(currentFolder,sub_name);
D = dir([string_path, '\*.m']);
N = length(D(not([D.isdir])));
legend_labels_front = cell(N,1);
legend_labels_rear = cell(N,1);

%--- plot all logfiles in an overlay ---%
figfront = figure('Name','LuxorRHI','units','normalized','outerposition',[0 0.2 1/3 2/3]); % front
figrear = figure('Name','LuxorRHI','units','normalized','outerposition',[1/3 0.2 1/3 2/3]); % rear

h_front = zeros(N+1,1);
h_rear = zeros(N+1,1);

dont_come_here_again = 1;
checking_dual = 0;
checking_omni = 0;
mismatch = 0;

for R = 1:1:N 
  logPath = strcat(currentFolder,sub_name);
  logFile = D(R).name;

  if R > 1
     clear DFSinit 
  end

  run([logPath,logFile]);                                                            
  fs = DFSInit.Samplerate;
  [H_front,f] = freqz(DFSInit.init_resp_front,1,fs,fs);

  if (strcmp(DFSInit.DFSType,'TypeOmniAndDual') || strcmp(DFSInit.DFSType,'TypeDual'))
    type_dual = 1;
    [H_rear,f_rear] = freqz(DFSInit.init_resp_rear,1,fs,fs);
    checking_dual = 1;
  else
    type_dual = 0;
    figrear.delete;
    checking_omni = 1;
  end

  temp = strfind(logFile, 'ref');
  figfront.Parent.CurrentFigure = 1;
  hold on;
  if (temp ~= 0) % reference
    if plot_calibrated_response 
      front_response = 20*log10(abs(H_front));
      cal_interpolated = interp1(fre,cal,f(:,1)); 
      front_response_calibrated = front_response + cal_interpolated;
      h_front(R,1) = plot(f(:,1),mirror.*front_response_calibrated,'linewidth',2);
    else
      h_front(R,1) = plot(f(:,1),mirror.*20*log10(abs(H_front)),'linewidth',2);
    end
  else
    if plot_calibrated_response 
      front_response = 20*log10(abs(H_front));
      cal_interpolated = interp1(fre,cal,f(:,1)); 
      front_response_calibrated = front_response + cal_interpolated;
      h_front(R,1) = plot(f(:,1), mirror.*front_response_calibrated);
    else
      h_front(R,1) = plot(f(:,1),mirror.*20*log10(abs(H_front)),'-','linewidth',1);%,'Color',[1 0 0]);
    end
    
    
  end
  hold all
  legend_name = logFile;
  index = strfind(legend_name, '_fblog'); % ignore ALT log naming
  legend_name_modified = legend_name(1:index(1)-1);
  legend_labels_front{R,1} = sprintf(legend_name_modified);

  if type_dual % plot dual plots

    if checking_dual && checking_omni && dont_come_here_again
      figrear = figure('Name','LuxorRHI','units','normalized','outerposition',[1/3 0.2 1/3 2/3]); % rear
      dont_come_here_again = 0;
      mismatch = 1;
      fprintf('There has been a mix between OMNI and DUAL log files, overlay therefore not plotted\n');
    end
    
    figfront.Parent.CurrentFigure = 2;
    temp = strfind(logFile, 'ref');
    hold on
    if (temp ~= 0)
      if plot_calibrated_response 
        rear_response = 20*log10(abs(H_rear));
        cal_interpolated = interp1(fre,cal,f(:,1)); 
        rear_response_calibrated = rear_response + cal_interpolated;
        h_rear(R,1) = plot(f_rear(:,1),mirror.*rear_response_calibrated,'linewidth',2);
      else
        h_rear(R,1) = plot(f_rear(:,1),mirror.*20*log10(abs(H_rear)),'linewidth',2);
      end
    else
      if plot_calibrated_response 
        rear_response = 20*log10(abs(H_rear));
        cal_interpolated = interp1(fre,cal,f(:,1)); 
        rear_response_calibrated = rear_response + cal_interpolated;
        h_rear(R,1) = plot(f(:,1), mirror.*rear_response_calibrated);
      else
        h_rear(R,1) = plot(f_rear(:,1),mirror.*20*log10(abs(H_rear)),'-','linewidth',1);%,'Color',[0.5 1 0.5]);
      end
    end
    
    % plot front and rear figure in same plot (as many plots as we have log files)
    if ~mismatch
      figtogether = figure;
      figtogether.set('units','normalized');
      k = R;
      if ((R*0.2)>1)
        k = 1; % do not plot outside screen
      end
      figtogether.Position = [2/3 0.2*R 1/3 3/(5*N)]; % x,y,w,h
      plot(f_rear(:,1),mirror.*20*log10(abs(H_rear)),'-','linewidth',1)
      hold on
      plot(f(:,1),mirror.*20*log10(abs(H_front)),'-','linewidth',1)
      xlim([100,fs/2])
      hold off
      grid on  
      legend('rear response','front response')
      tit = title(legend_labels_front{R});    
      fig = gcf;
      set(tit,'Interpreter','none');
    end
    
    hold all
    legend_name_rear = logFile;
    index = strfind(legend_name_rear, '_fblog'); % ignore ALT log naming
    legend_name_modified_rear = legend_name_rear(1:index(1)-1);
    legend_labels_rear{R,1} = sprintf(legend_name_modified_rear);     
  end

end


figfront.Parent.CurrentFigure = 1;
h_front(N+1,1) = plot(frq_FOG,FOG);
legend_labels_front{N+1,1} = 'IG';
xlim([100,fs/2])
ylim([40,120])
xlabel('Frequency [Hz]')
ylabel('Magnitude[dB]')
title('Front Initialization')  
leg = legend(h_front,legend_labels_front{:}); 
set(leg,'Interpreter','none'); % do not treat underscore as a control character 
set(leg,'FontSize',10);
leg.Location = 'best';
set(gcf,'color','w');
hold off; grid on; box on;

if type_dual
  figrear.Parent.CurrentFigure = 2;
  legend_labels_rear{N+1,1} = 'IG';
  h_rear(N+1,1) = plot(frq_FOG,FOG);
  hold off; grid on; box on;
  xlim([100,fs/2])
  ylim([40,120])
  xlabel('Frequency [Hz]')
  ylabel('Magnitude[dB]')
  title('Rear Initialization')  
  skip_rest = 0;
  for k = 1:1:N
    if isempty(legend_labels_rear{k})
      skip_rest = 1;    
    end
  end
  
  if skip_rest ~= 1
    leg_rear = legend(h_rear,legend_labels_rear{:});
    set(leg_rear,'Interpreter','none'); % do not treat underscore as a control character 
    set(leg_rear,'FontSize',10);
    leg_rear.Location = 'best';
    set(gcf,'color','w');
  end
end  

%% do the MSG comparision between ALT and OPL
temp=importdata('validate this script.xlsx');
f_opl=temp.data(1:end-5,1);
gain0=temp.data(1:end-5,3);
MSG_alt=interp1(f,mirror.*front_response_calibrated,f_opl)+gain0;
MSG_opl=temp.data(1:end-5,5);
fig_msg = figure('Name','two MSG','units','normalized','outerposition',[0 0.2 1/3 2/3]); % msg
plot(f_opl,MSG_opl);
hold on;plot(f_opl,MSG_alt,'-r');
legend('MSG\_OPL','MSG\_ALT');
DIFF_MSG=MSG_opl-MSG_alt;
fig_diff = figure('Name','diff of msg','units','normalized','outerposition',[0 0.2 1/3 2/3]); % diff of msg
plot(f_opl,DIFF_MSG)
legend('diff of MSG');