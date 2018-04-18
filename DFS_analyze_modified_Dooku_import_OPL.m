close all;
clear;
clc;

mirror = -1; % set to -1 if gain curve is desired instead of attenuation curve
plot_calibrated_response = 1; % set to 1 if calibrations should be taken into account
sub_name = '\DFSout\';
production_calibration = 0; % output cal offset, just copy the number here

%% System cal in Palpatine
%

% palpatine2 Delhi80
cal =[-23 -23 -21 -21 -25 -27 -27 -27 -21 -24 -18 -13 -14 -17 -11 -5 -5];

% No system cal at all
%cal = [0	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

cal = cal + production_calibration;
% fre = [0 168 340 518 706 908 1131 1383 1674 2019 2442 2981 3693 4672 6059 7997 10417]; % High bandwitdth frequencies
fre = [0 163 329 501 682 876 1090 1328 1600 1919 2301 2773 3370 4140 5136 6383 7812]; % High bandwitdth frequencies

% frq_FOG=[250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000,8000]; %Dooku
frq_FOG=[250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000]; %palpatine

FOG=[57 59 59 61 62 58 56 51 51];  %Delhi80
%FOG=[51, 51, 51, 51, 51, 48, 48, 48, 48, 48];  %Luxor NP
%FOG=[50, 50, 50, 50, 54, 56, 53, 42, 42, 42]; %HP
%FOG=[43, 43, 43, 42, 43, 43, 42, 41, 41, 41]; %LP
color_array = {'-b','-r','-c','-m','-g','-y','-k','--b','--r','--c','--m','--g','--y','--k'};

%% System cal in Dooku1 BerlinRIE NP
% %
% cal = [-3	-3 2 6 6 8 8 9 9 10 10 13 19 17 13 7 7];
% 
% %Dooku1 Berlin RIE HP
% %cal = [2 2 6 6 6 6 5 2 3 4 8 11 16 16 10 7 7];
% 
% % No system cal at all
% %cal = [0	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% 
% %cal = [7 7 11 12 13 13 12 13 13 13 14 16 19 17 14 8 8];
%  
% 
% cal = cal + production_calibration;
% fre = [0 168 340 518 706 908 1131 1383 1674 2019 2442 2981 3693 4672 6059 7997 10417]; % High bandwitdth frequencies
% 
% frq_FOG=[250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000,8000]; %Dooku
% 
% FOG=[51, 51, 51, 51, 51, 48, 48, 48, 48, 48];  %Luxor NP
% %FOG=[50, 50, 50, 50, 54, 56, 53, 42, 42, 42]; %HP
% %FOG=[43, 43, 43, 42, 43, 43, 42, 41, 41, 41]; %LP
% color_array = {'-b','-r','-c','-m','-g','-y','-k','--b','--r','--c','--m','--g','--y','--k'};

%% plot more OMNI curves together
currentFolder = pwd;
currentFolder = strcat(currentFolder,'\');
string_path = strcat(currentFolder,sub_name);
D = dir([string_path, '\*.m']);
N = length(D(not([D.isdir])));
legend_labels_front = cell(N,1);
legend_labels_rear = cell(N,1);

%--- plot all logfiles in an overlay ---%
% model=input('please input the model as the figure title: ','s');
model = 'Delhi80';
figfront = figure('Name',model,'units','normalized','outerposition',[0 0.2 1/3 2/3]); % front
figrear = figure('Name',model,'units','normalized','outerposition',[1/3 0.2 1/3 2/3]); % rear

h_front = zeros(N+1,1);
h_rear = zeros(N+1,1);

dont_come_here_again = 1;
checking_dual = 0;
checking_omni = 0;
mismatch = 0;



%% import FOG data and gain0
cd D:\MATLAB\work\DFS\FOG
GAIN_D=dir([pwd,'\*.txt']);
GAIN_N = length(GAIN_D(not([GAIN_D.isdir])));
for data_i=1:GAIN_N
    temp_s=GAIN_D(data_i).name;
    index1=strfind(temp_s,'_I'); %text file name as the variable name
    eval([temp_s(1:3) temp_s(5:index1-1) '=importdata(temp_s);']) %import from text file
end
cd gain0
GAIN_D=dir([pwd,'\*.txt']);
GAIN_N = length(GAIN_D(not([GAIN_D.isdir])));
for data_i=1:GAIN_N
    temp_s=GAIN_D(data_i).name;
    index1=strfind(temp_s,'_I'); %text file name as the variable name
    eval([temp_s(1:index1-1) '=importdata(temp_s);']) %import from text file
end
cd ..
cd ..
%%
% import OPL data
cd D:\MATLAB\work\DFS\OPL
OPL_D=dir([pwd,'\*.txt']);
OPL_N = length(OPL_D(not([OPL_D.isdir])));
for data_i=1:OPL_N
    temp_s=OPL_D(data_i).name;
    index1=strfind(temp_s,'_I'); %text file name as the variable name
    index2=strfind(temp_s,'D80');
    eval([temp_s(1:index1-1) temp_s(index2+3:index2+6) '=importdata(temp_s);']) %import from text file
    eval(['temp=' temp_s(1:index1-1) temp_s(index2+3:index2+6) ';'])
    index_front=strfind(temp_s,'Front_');
    index_rear=strfind(temp_s,'Rear_');
    if ~isempty(index_front)
        figure(4);hold on;title('Front OPL');
        h_f=plot(temp(:,1),FOG20F(:,2)+20-temp(:,2),'k','linewidth',2);
    elseif ~isempty(index_rear)
        figure(5);hold on;title('Rear OPL');
        h_r=plot(temp(:,1),FOG20R(:,2)+20-temp(:,2),'k','linewidth',2);
    end
end
cd ..

%%
% start to handle with the DFS data
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
%     figrear.delete; %%%%%%%%%% error for old version, figrear is double type, not struct
    checking_omni = 1;
  end

  temp = strfind(logFile, 'ref');
  %figfront.Parent.CurrentFigure = 1;
  figure(1)
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
      h_front(R,1) = plot(f(:,1), mirror.*front_response_calibrated, color_array{R});
    else
      h_front(R,1) = plot(f(:,1),mirror.*20*log10(abs(H_front)),'-','linewidth',1);%,'Color',[1 0 0]);
    end
    
    
  end
  
  % add msg_dfs_off curves into the initialization curves
  hold on;
  plot(1:f(end)/128:f(end),DFSInit.front.msg_dfs_off,'-r','LineWidth',2);
  
  hold all
  legend_name = logFile;
  index = strfind(legend_name, '_fblog'); % ignore ALT log naming
  legend_name_modified = legend_name(1:index(1)-1);
  legend_labels_front{R,1} = sprintf(legend_name_modified);
  
  %%%%%%% add the attenu curve to OPL_front
  figure(4);hold on;
  h_f1=plot(f(:,1), mirror.*front_response_calibrated,'-r','LineWidth',2);

  if type_dual % plot dual plots

    if checking_dual && checking_omni && dont_come_here_again
      figrear = figure('Name',model,'units','normalized','outerposition',[1/3 0.2 1/3 2/3]); % rear
      dont_come_here_again = 0;
      mismatch = 1;
      fprintf('There has been a mix between OMNI and DUAL log files, overlay therefore not plotted\n');
    end
    
    %figfront.Parent.CurrentFigure = 2;
    figure(2)
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
      % add msg_dfs_off curves into the initialization curves
      hold on;
      plot(1:f(end)/128:f(end),DFSInit.rear.msg_dfs_off,'-r','LineWidth',2);
  
    else
      if plot_calibrated_response 
        rear_response = 20*log10(abs(H_rear));
        cal_interpolated = interp1(fre,cal,f(:,1)); 
        rear_response_calibrated = rear_response + cal_interpolated;
        h_rear(R,1) = plot(f(:,1), mirror.*rear_response_calibrated,color_array{R});
      else
        h_rear(R,1) = plot(f_rear(:,1),mirror.*20*log10(abs(H_rear)),'-','linewidth',1);%,'Color',[0.5 1 0.5]);
      end
      % add msg_dfs_off curves into the initialization curves
      hold on;
      plot(1:f(end)/128:f(end),DFSInit.rear.msg_dfs_off,'-r','LineWidth',2);
    end
    
    % plot front and rear figure in same plot (as many plots as we have log files)
    if ~mismatch
      %figtogether = figure;
        figure(3)
      %figtogether.set('units','normalized');
      k = R;
      if ((R*0.2)>1)
        k = 1; % do not plot outside screen
      end
      

      
      figtogether.Position = [2/3 0.2*R 1/3 3/(5*N)]; % x,y,w,h
      plot(f_rear(:,1),mirror.*20*log10(abs(H_rear)),'-b','linewidth',1)
      hold on
      plot(f(:,1),mirror.*20*log10(abs(H_front)),'-r','linewidth',1)
      hold on
      plot(1:f(end)/128:f(end),DFSInit.front.msg_dfs_off,'-k','LineWidth',2);%msg_dfs_off is the MSG without system calibration
      hold on
      plot(1:f(end)/128:f(end),DFSInit.rear.msg_dfs_off,'-c','LineWidth',2);
      xlim([100,fs/2])
      hold off
      grid on  
      legend('rear response','front response','msg\_dfs\_off\_front','msg\_dfs\_off\_rear')
%       tit = title(legend_labels_front{R});  
      tit = title('msg\_dfs\_off vs freq response without system cal');
      fig = gcf;
      set(tit,'Interpreter','none');
      % make the figure full screen
      set(gcf,'units','normalized','position',[0,0,1,1]);
    end
    
    hold all
    legend_name_rear = logFile;
    index = strfind(legend_name_rear, '_fblog'); % ignore ALT log naming
    legend_name_modified_rear = legend_name_rear(1:index(1)-1);
    legend_labels_rear{R,1} = sprintf(legend_name_modified_rear);   
      
    %%%%%%% add the attenu curve to OPL_rear
    figure(5);hold on;
    h_r1=plot(f(:,1),rear_response_calibrated,'-r','linewidth',2);
  end

end


%figfront.Parent.CurrentFigure = 1;
figure(1)
hold on
h_front(N+1,1) = plot(frq_FOG,FOG);
legend_labels_front{N+1,1} = 'IG';
xlim([100,fs/2])
ylim([20,100])
xlabel('Frequency [Hz]')
ylabel('Magnitude[dB]')
title('Front Initialization')  
leg = legend(h_front,legend_labels_front{:}); 
set(leg,'Interpreter','none'); % do not treat underscore as a control character 
set(leg,'FontSize',10);
leg.Location = 'best';
set(gcf,'color','w','position',[0,0,1,1]);
hold off; grid on; box on;


if type_dual
  %figrear.Parent.CurrentFigure = 2;
  figure(2)
  legend_labels_rear{N+1,1} = 'IG';
  h_rear(N+1,1) = plot(frq_FOG,FOG);
  hold off; grid on; box on;
  xlim([100,fs/2])
  ylim([20,100])
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
    set(gcf,'color','w','position',[0,0,1,1]);
  end
end  


%% save the figure
figure(1);
saveas(gcf,[model '_front_msg_dfs_off'],'bmp');
figure(2);
saveas(gcf,[model '_rear_msg_dfs_off'],'bmp');
figure(3);
saveas(gcf,[model '_all_freq_res_wo_cal'],'bmp');
% figure(4);
% ylim([-100,0])
legend([h_f,h_f1],'MSG_OPL','from feedback path');
% make the figure full-screen
set(gcf,'units','normalized','position',[0,0,1,1]);
grid on;
saveas(gcf,[model '_front_OPL_atten'],'bmp');

figure(5);
ylim([-100,0])
% legend([h_r,h_r1],'OPL','Attenu of feedback path');
grid on;
% make the figure full-screen
set(gcf,'units','normalized','position',[0,0,1,1]);
saveas(gcf,[model '_rear_OPL_atten'],'bmp');

% close all;