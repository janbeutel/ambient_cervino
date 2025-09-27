function CERVINO_event_detection(year, station, channel)
   % year comes in as a number (e.g. 2018)
   year = string(year);

   close all
   clearvars -except year station channel  % don’t clear year

   % year    = "2018";
   network = "1I";
   % station = "MH44";
   % channel = "EHE.D";

   % fprintf("Running CERVINO_event_detection for year %s, station %s, channel %s\n", year, station, channel)

   data_directory = "../data/" + network + "/" + station + "/" + year + "/" + channel;
   % List all .miniseed files in the directory
   FileList = dir(fullfile(data_directory, '*.mat'));
   TOT = size(FileList,1);

   %PSD_all=zeros(16385,TOT);
   % DVec=zeros(TOT,6);
   % first=FileList(1).name;

   fprintf("Running CERVINO_event_detection for year %s, station %s, channel %s, and %d files\n", year, station, channel, TOT);

   savedir = "../results/events/" + network + "/"+ station + "/" + year + "/" + channel;
   if ~exist(savedir, 'dir')  % check if folder exists
      mkdir(savedir);         % create folder if it doesn't exist
   end

   % Define log file path
   logfile = "../results/events/" + network + "/" + station + "/" + year + "/" + channel + "_" + year + ".log";
   % Remove logfile if it exists
   if exist(logfile, "file")
      delete(logfile);
   end
   % Open in append mode ("a")
   fid = fopen(logfile, "a");
   if fid == -1
      error("Could not open log file.");
   end
   % Write header logfile
   fprintf(fid, "matfilename,numevents,eventfilename\n");

   for ite = 1:TOT
      % get the file name:
      folder = FileList(ite).folder;
      filename = FileList(ite).name;
      fprintf("Processing %s/%s ", folder, filename);

      % build full path
      fullpath = fullfile(folder, filename);
      load(fullpath); 

      L=length(sig);
      N=sig;
      Nm=mean(N);
      trace_v=N;

      year=str2num(filename(17:20));
      month=str2num(filename(21:22));
      day=str2num(filename(23:24));
      hour=str2num(filename(26:27));
            
      DVec(ite,1:4)=[year month day hour];
      DVec(ite,5:6)=0;
      
      %e=etime(DVec(ite,:),DV_ref);
      
      %time=T(end);
      x=linspace(datenum(DVec(ite,:)),datenum(DVec(ite,:))+(1/24),3600*FS);

      if L==3600*FS
         %%%%%%%%%%% STA/LTA    
         edp(1)=0.5;   %%%%%%Short time window (s)
         edp(2)=10;  %%%%%%Long time window (s)
         edp(3)=6; %%%%%% STA/LTA trigger on
         edp(4)=4; %%%%%%% STA/LTA trigger off
         edp(5)=1; %%%%%%%%% skip after end of the event
         edp(6)=0.5; %%%%%%%%%% min duration
         edp(7)=3; %rec time before trigger on
         edp(8)=10; %rec time after trigger off


         % type='wfa';

         v = trace_v(:,1);          % Waveform data
         Fs=FS;
         l_v = L; % Length of time series
         tv = x;  % Time vector of waveform
         abs_v = abs(v);                % Absolute value of time series

         %% Event Detection Parameters
         l_sta = edp(1)*Fs;   % STA window size
         l_lta = edp(2)*Fs;   % LTA window size
         ll_lta = l_lta;      % Current length of growing LTA window
         th_on = edp(3);      % Trigger on when sta_to_lta exceeds this theshold
         th_off = edp(4);     % Trigger off when sta_to_lta drops below threshold
         % for off_dur data points
         skip_int = edp(5)*Fs;% Skip ahead after end of event (seconds)
         min_dur = edp(6)*Fs; % Any triggers longer than min_dur become events

         %% Initialize flags and other variables
         lta_calc_flag = 0;       % has the full LTA window been calculated?
         ntrig = 0;               % number of triggers
         trig_array = zeros(1,2); % array of trigger times: [on,off;on,off;...]

         %% Loops over data
         % i is the primary reference point (right end of STA/LTA window)
         i = l_lta+1;
         while i <= l_v

            %% Skip data gaps (NaN values in LTA window)?
            if any(isnan(abs_v(i-l_lta:i)))
               gap = 1;
               lta_calc_flag = 0; % Force full claculations after gap
               while gap == 1
                  if i <= l_v, i = i+1; end
                  if ~any(isnan(abs_v(i-l_lta:i)))
                     gap = 0;
                  end
               end
            end

            %% Calculate STA & LTA Sum (Full)?
            if (lta_calc_flag == 0)
               lta_sum = 0;
               sta_sum = 0;
               for j = i-l_lta:i-1              % Loop to compute LTA & STA
                  lta_sum = lta_sum + abs_v(j); % Sum LTA window
                  if (i - j) <= l_sta           % Sum STA window (right side of LTA)
                     sta_sum = sta_sum + abs_v(j);
                  end
               end
               lta_calc_flag = 1;
            else
               %% Calculate STA & LTA Sum (Single new data point) if not Full
               lta_sum = lta_sum - abs_v(i-l_lta-1) + abs_v(i-1);
               sta_sum = sta_sum - abs_v(i-l_sta-1) + abs_v(i-1);
            end

            %% Calculate STA & LTA
            lta = lta_sum/l_lta;
            sta = sta_sum/l_sta;

            %% Calculate STA/LTA Ratio
            sta_to_lta = sta/lta;

            %% Trigger on? (Y/N)
            if sta_to_lta > th_on
               j = i;        % Set secondary reference point = primary
               while (sta_to_lta > th_off) % Track STA while LTA left edge static
                  j = j+1;
                  if j < l_v
                     sta_sum = sta_sum - abs_v(j-l_sta-1) + abs_v(j-1);
                     lta_sum = lta_sum + abs_v(j-1);
                     ll_lta = ll_lta+1;
                     sta = sta_sum/l_sta;
                     lta = lta_sum/ll_lta;
                     sta_to_lta = sta/lta;
                     
                     %% Skip data gaps (NaN values in STA window)?
                     if any(isnan(abs_v(j-l_sta:j)))
                        sta_to_lta = 0; % Force trigger off (data gap in STA window)
                     end
                     
                  else
                     sta_to_lta = 0; % Force trigger off (end of data)
                  end
               end
               duration = (j-i);

               %% Triggered period long enough? (Y/N)
               if duration > min_dur % If duration < min_dur then skip it
                  trig_t = i-l_sta;  % Beginning of STA window during trigger on
                  end_t  = j;        % End of STA window during trigger off
                  ntrig = ntrig + 1; % Event counter
                  trig_array(ntrig,:) = [trig_t, end_t];
               end
               ll_lta = l_lta;
               i = j + skip_int;  % Skip ahead
               lta_calc_flag = 0; % Reset LTA calc flag to force new computation
            end
            i = i + 1;
         end

         %% Return events
         if (trig_array(1,1)==0)&&(trig_array(1,2)==0)
            fprintf('\n No events detected\n')
         else
            fprintf("- %d events detected\n", size(trig_array,1))

            % switch lower(type)
            %    case {'ssd'}
            %       events = trig_array;
            %    case {'sst'}
            %       events = tv(trig_array);
            %    case {'wfa'}
            for n=1:size(trig_array,1)
               %fix = [2 5]; Fix event length to 2 seconds before trig on, 5 seconds
               %after trig on (All event are 7 seconds long), Change 'fix' as needed
               %events=[events;extract(wave,'INDEX',trig_array(n,1)-Fs*fix(1),...
               %                                    trig_array(n,1)+Fs)*fix(2)];
               
               if trig_array(n,1)-FS*edp(7)>0 & trig_array(n,2)+FS*edp(8)<3600*FS
                  X=x(trig_array(n,1)-FS*edp(7):trig_array(n,2)+FS*edp(8));
                  SIG=trace_v(trig_array(n,1)-FS*edp(7):trig_array(n,2)+FS*edp(8),:);
               else
                  X=x(trig_array(n,1):trig_array(n,2));
                  SIG=trace_v(trig_array(n,1):trig_array(n,2),:);   
               end
               sec=(trig_array(n,1)/FS);
               % sec_name=num2str(round(sec));
               sec_name = sprintf('%04d', round(sec));
               savename=[filename(1:31) + "_" + sec_name + ".mat"];
               % Write line to logfile
               fprintf(fid, "%s,%d,%s\n", filename, size(trig_array,1), savename);
               
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               save(savedir + "/" + savename,'SIG','X','sec', 'FS')            
            end
         end
      end
      clearvars -except FS FileList ite TOT network station year channel fid savedir
   end
   % Close logfile
   fclose(fid);
end
