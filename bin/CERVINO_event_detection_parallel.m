function CERVINO_event_detection_parallel(year, station, channel)
% CERVINO_event_detection detects events in seismic data using STA/LTA.
% Inputs:
%   year    - string or number, e.g. '2018' or 2018
%   station - string, e.g. 'MH44'
%   channel - string, e.g. 'EHE.D'

   year = string(year);
   close all
   clearvars -except year station channel

   network = "1I";
   data_directory = "../data/" + network + "/" + station + "/" + year + "/" + channel;
   FileList = dir(fullfile(data_directory,'*.mat'));
   TOT = numel(FileList);

   fprintf("Running CERVINO_event_detection for year %s, station %s, channel %s, with %d files\n", ...
      year, station, channel, TOT);

   % Output directories
   savedir = "../results/events_parallel/" + network + "/" + station + "/" + year + "/" + channel;
   if ~exist(savedir,'dir'), mkdir(savedir); end

   logfile = "../results/events_parallel/" + network + "/" + station + "/" + year + "/" + channel + "_" + year + ".csv";
   if exist(logfile,'file'), delete(logfile); end

   % Predefine edp (event detection parameters)
   edp = [0.5, 10, 6, 4, 1, 0.5, 3, 10];  % [STA(s), LTA(s), th_on, th_off, skip_after(s), min_duration(s), rec_before(s), rec_after(s)]

   % Preallocate DVec
   DVec = zeros(TOT,6);

   % Pre-extract folders and filenames for parfor safety
   folders = {FileList.folder};
   filenames = {FileList.name};

   parfor ite = 1:TOT
      folder = folders{ite};
      filename = filenames{ite};
      fullpath = fullfile(folder, filename);

      try
         out = process_event_file(fullpath, savedir, edp, filename);

         % Store timestamp
         DVec(ite,:) = out.DVec_row;

         % Write CSV via temporary files to avoid parfor conflicts
         if ~isempty(out.trig_info)
               tempFile = fullfile(savedir, "temp_log_" + ite + ".txt");
               fid_tmp = fopen(tempFile,'w');
               fprintf(fid_tmp,'%s\n',out.trig_info{:});
               fclose(fid_tmp);
         end
      catch ME
         fprintf("Error processing file %s: %s\n", filename, ME.message);
      end
   end

   % Merge temporary CSV files
   tempFiles = dir(fullfile(savedir,'temp_log_*.txt'));
   fid_final = fopen(logfile,'a');

   allLines = {};
   for k = 1:numel(tempFiles)
      fid_tmp = fopen(fullfile(savedir,tempFiles(k).name),'r');
      tline = fgetl(fid_tmp);
      while ischar(tline)
         allLines{end+1,1} = tline;  %#ok<SAGROW>
         tline = fgetl(fid_tmp);
      end
      fclose(fid_tmp);
      delete(fullfile(savedir,tempFiles(k).name)); % clean up
   end

   % Convert lines to table for sorting by 3rd column (eventfilename)
   numLines = numel(allLines);
   C = cell(numLines,8);  % 8 columns
   for k = 1:numLines
      C(k,:) = strsplit(allLines{k}, ',');
   end

   T = cell2table(C, 'VariableNames', ...
      {'matfilename','numevents','eventfilename','detectiontime','windowlength','peak_amp','rms_amp','triggerduration'});

   % Sort table by eventfilename
   T_sorted = sortrows(T, 'eventfilename');

   % Write final CSV
   fid_final = fopen(logfile,'w');
   fprintf(fid_final, "matfilename,numevents,eventfilename,detectiontime,windowlength,peak_amp,rms_amp,triggerduration\n");

   for k = 1:height(T_sorted)
      fprintf(fid_final, '%s,%s,%s,%s,%s,%s,%s,%s\n', ...
         T_sorted.matfilename{k}, T_sorted.numevents{k}, T_sorted.eventfilename{k}, ...
         T_sorted.detectiontime{k}, T_sorted.windowlength{k}, T_sorted.peak_amp{k}, ...
         T_sorted.rms_amp{k}, T_sorted.triggerduration{k});
   end

   fclose(fid_final);

   fprintf("Event detection finished for year %s, station %s, channel %s.\n", year, station, channel);

end



function out = process_event_file(fullpath, savedir, edp, filename)
% Returns a struct with:
%   .DVec_row   - timestamp from filename
%   .trig_info  - cell array of CSV lines
%   .sigseconds - vector of event start times (in seconds)

   fprintf("Processing %s ", fullpath)

   S = load(fullpath);
   sig = S.sig;
   t = S.t;
   FS = S.FS;
   L = length(sig);

   % Extract timestamp from filename
   year  = str2double(filename(17:20));
   month = str2double(filename(21:22));
   day   = str2double(filename(23:24));
   hour  = str2double(filename(26:27));
   DVec_row = [year month day hour 0 0];

   trig_info = {};
   sigseconds_list = [];

   if L ~= 3600*FS
      out = struct('DVec_row', DVec_row, 'trig_info', trig_info, 'sigseconds', sigseconds_list);
      return;
   end

   % Signal prep
   abs_v = abs(sig(:,1));
   l_v = L;

   % STA/LTA parameters
   l_sta = edp(1)*FS;
   l_lta = edp(2)*FS;
   th_on = edp(3);
   th_off = edp(4);
   skip_int = edp(5)*FS;
   min_dur = edp(6)*FS;
   rec_before = edp(7)*FS;
   rec_after  = edp(8)*FS;
   ll_lta = l_lta;

   lta_calc_flag = 0;
   trig_array = zeros(0,2);

   i = l_lta + 1;
   while i <= l_v
      if any(isnan(abs_v(i-l_lta:i)))
         lta_calc_flag = 0;
         while i<=l_v && any(isnan(abs_v(i-l_lta:i))), i=i+1; end
      end

      if lta_calc_flag==0
         lta_sum = sum(abs_v(i-l_lta:i-1));
         sta_sum = sum(abs_v(i-l_sta:i-1));
         lta_calc_flag=1;
      else
         lta_sum = lta_sum - abs_v(i-l_lta-1) + abs_v(i-1);
         sta_sum = sta_sum - abs_v(i-l_sta-1) + abs_v(i-1);
      end

      lta = lta_sum/l_lta;
      sta = sta_sum/l_sta;
      sta_to_lta = sta/lta;

      if sta_to_lta > th_on
         j = i;
         while j < l_v && sta_to_lta>th_off
               j=j+1;
               sta_sum = sta_sum - abs_v(j-l_sta-1)+abs_v(j-1);
               lta_sum = lta_sum + abs_v(j-1);
               ll_lta = ll_lta + 1;
               sta = sta_sum/l_sta;
               lta = lta_sum/ll_lta;
               sta_to_lta = sta/lta;
               if any(isnan(abs_v(max(1,j-l_sta):j))), sta_to_lta=0; end
         end

         duration = j - i;
         if duration>min_dur
               trig_array(end+1,:) = [i-l_sta,j];
         end
         ll_lta = l_lta;
         i = j + skip_int;
         lta_calc_flag = 0;
      end
      i=i+1;
   end

   % Save events & prepare CSV lines
   fprintf("- %d events detected\n", length(trig_array))

   for n=1:size(trig_array,1)
      trig_start = max(1,trig_array(n,1)-rec_before);
      trig_end   = min(l_v, trig_array(n,2)+rec_after);      

      % % compute relative seconds for signal segment
      % X = (trig_start:trig_end)/FS;  
      % SIG = sig(trig_start:trig_end,:);

      % % build base datetime from filename (year, month, day, hour)
      % base_dt = datetime(year, month, day, hour, 0, 0);

      % % actual detection time
      % dt = base_dt + seconds(trig_array(n,1)/FS);

      % X = (trig_start:trig_end)/FS
      % SIG = sig(trig_start:trig_end,:);
      % dt = datetime(X(1),'ConvertFrom','datenum')

      % --- Build 1-hour absolute time vector (datenum) ---
      x = linspace(datenum([year month day hour 0 0]), datenum([year month day hour+1 0 0]), 3600*FS);

      % --- Extract triggered segment using the same logic as original code ---
      if (trig_array(n,1) - rec_before > 0) && (trig_array(n,2) + rec_after < 3600*FS)
         X = x(trig_array(n,1)-rec_before : trig_array(n,2)+rec_after);
         SIG = sig(trig_array(n,1)-rec_before : trig_array(n,2)+rec_after, :);
      else
         X = x(trig_array(n,1) : trig_array(n,2));
         SIG = sig(trig_array(n,1) : trig_array(n,2), :);
      end

      % --- Detection time (datetime) ---
      dt = datetime(X(1), 'ConvertFrom', 'datenum');

      sigseconds = trig_array(n,1)/FS;
      windowlength = (trig_array(n,2)-trig_array(n,1))/FS;
      peak_amp = max(abs(SIG));
      rms_amp = rms(SIG);
      thr = 0.1*max(abs(SIG));
      idx = find(abs(SIG)>thr);
      if isempty(idx)
         triggerduration=0; 
      else 
         triggerduration=t(idx(end))-t(idx(1));
      end

      sec_name = sprintf('%04d', round(sigseconds));
      savename = filename(1:31) + "_" + sec_name + ".mat";
      save(fullfile(savedir,savename),"SIG","X","FS","sigseconds")

      % Append CSV line
      trig_info{end+1} = sprintf("%s,%d,%s,%s,%d,%d,%d,%d", ...
         filename, size(trig_array,1), savename, datestr(dt,'yyyy-mm-dd HH:MM:SS'), ...
         windowlength, peak_amp, rms_amp, triggerduration);

      % Collect sigseconds
      sigseconds_list(end+1) = sigseconds;
   end

   out = struct('DVec_row', DVec_row, 'trig_info', {trig_info}, 'sigseconds', sigseconds_list);
end
