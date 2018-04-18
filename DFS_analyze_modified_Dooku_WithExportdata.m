close all;
clear;
clc;

mirror = -1; % set to -1 if gain curve is desired instead of attenuation curve
plot_calibrated_response = 1; % set to 1 if calibrations should be taken into account
sub_name = '\DFSout\';
production_calibration = 0; % output cal offset, just copy the number here


%% System cal in Dooku1 BerlinRIE NP
%
%cal = [-3	-3 2 6 6 8 8 9 9 10 10 13 19 17 13 7 7]; %nP before Dooku .14
%cal = [-3	-3 2 6 6 8 6 7 7 8 6 11 20 20 18 12 7];

%Dooku1 Berlin RIE HP
%cal = [2 2 6 6 6 6 5 2 3 4 8 11 16 16 10 7 7];

% No system cal at all
%cal = [0	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%% System cal in Dooku1 BerlinRIE LP
%cal = [7 7 11 12 13 13 12 13 13 13 14 16 19 17 14 8 8];


%% System cal in Dooku1 BerlinRIE UP
%cal = [-7	-7 -3 -3 -3 -3 -6 -10 -10 -4 4 11 14 8 6  2 2] %before D14
%cal = [-8	-7 -1 -2 -3 -4 -8 -9 -14 -10 2 14 10 7 11  5 5] %.19
% cal = [-6 -5 -1 -2 -2 -3 -6 -10 -10 -6 4 11 14 8 7 2 2]; %PA3 final
cal = [-10	-7	-3	-4	-3	-5	-6	-11	-10	-6	4	9	14	7	5	1	1]; % CAM13 UP

cal = cal + production_calibration;
fre = [0 168 340 518 706 908 1131 1383 1674 2019 2442 2981 3693 4672 6059 7997 10417]; % High bandwitdth frequencies
frq_FOG=[250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000,8000]; %Dooku

%FOG=[51, 51, 51, 51, 51, 48, 48, 48, 48, 48];  %Luxor NP
%FOG=[50, 50, 50, 50, 54, 56, 53, 42, 42, 42]; %HP
%FOG=[43, 43, 43, 42, 43, 43, 42, 41, 41, 41]; %LP
%FOG=[55, 55, 55, 57, 76, 79, 48, 54, 49, 46]; %UP 19.0 earlier
% FOG=[55, 55, 55, 68, 74, 76, 48, 57, 49, 40]; %UP PA3 final
FOG=[55	55	55	61	74	76	48	55	48	40]; % CAM13 UP

color_array = {'-b','-r','-c','-m','-g','-y','-k','--b','--r','--c','--m','--g','--y','--k'};
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
    else
      if plot_calibrated_response 
        rear_response = 20*log10(abs(H_rear));
        cal_interpolated = interp1(fre,cal,f(:,1)); 
        rear_response_calibrated = rear_response + cal_interpolated;
        h_rear(R,1) = plot(f(:,1), mirror.*rear_response_calibrated,color_array{R});
      else
        h_rear(R,1) = plot(f_rear(:,1),mirror.*20*log10(abs(H_rear)),'-','linewidth',1);%,'Color',[0.5 1 0.5]);
      end
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


%figfront.Parent.CurrentFigure = 1;
figure(1)
hold on
h_front(N+1,1) = plot(frq_FOG,FOG);
legend_labels_front{N+1,1} = 'IG';
xlim([100,4000])
ylim([50,120])
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
  %figrear.Parent.CurrentFigure = 2;
  figure(2)
  legend_labels_rear{N+1,1} = 'IG';
  h_rear(N+1,1) = plot(frq_FOG,FOG);
  hold off; grid on; box on;
  xlim([100,4000])
  ylim([50,120])
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

%% export the data to an excel file

figure(1);
hl=findall(gca,'type','line');
for R=1:size(hl)
    if R==1
        xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},1,'b1');
        xlswrite([pwd, '\data.xlsx'],hl(R).XData',1,'a2');
        xlswrite([pwd, '\data.xlsx'],hl(R).YData',1,'b2');
    elseif R==2
        xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},1,[char('c'+R-1) '1']);
        xlswrite([pwd, '\data.xlsx'],hl(R).XData',1,'c2');
        xlswrite([pwd, '\data.xlsx'],hl(R).YData',1,'d2');
    else
        xlswrite([pwd, '\data.xlsx'],hl(R).YData',1,[char('d'+R-2) '2']);
        xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},1,[char('d'+R-2) '1']);
    end
    
end

if type_dual==1
    figure(2);
    hl=findall(gca,'type','line');
    for R=1:size(hl)
        if R==1
            xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},2,'b1');
            xlswrite([pwd, '\data.xlsx'],hl(R).XData',2,'a2');
            xlswrite([pwd, '\data.xlsx'],hl(R).YData',2,'b2');
        elseif R==2
            xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},2,[char('c'+R-1) '1']);
            xlswrite([pwd, '\data.xlsx'],hl(R).XData',2,'c2');
            xlswrite([pwd, '\data.xlsx'],hl(R).YData',2,'d2');
        else
            xlswrite([pwd, '\data.xlsx'],hl(R).YData',2,[char('d'+R-2) '2']);
            xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},2,[char('d'+R-2) '1']);
        end

    end
end

% for R=1:size(hl)
%     xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},1,['a' num2str(4*(R-1)+1)]);
%     xlswrite([pwd, '\data.xlsx'],hl(R).XData,1,['a' num2str(4*(R-1)+2)]);
%     xlswrite([pwd, '\data.xlsx'],hl(R).YData,1,['a' num2str(4*(R-1)+3)]);
% end
% 
% if type_dual==1
%     figure(2);
%     hl=findall(gca,'type','line');
%     for R=1:size(hl)
%         xlswrite([pwd, '\data.xlsx'],{hl(R).DisplayName},2,['a' num2str(4*(R-1)+1)]);
%         xlswrite([pwd, '\data.xlsx'],hl(R).XData,2,['a' num2str(4*(R-1)+2)]);
%         xlswrite([pwd, '\data.xlsx'],hl(R).YData,2,['a' num2str(4*(R-1)+3)]);
%     end
% end

