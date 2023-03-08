function [Traces] = fetchData(E, S, pre, post, phase, tag, taper_fraction, sample_rate, ...
    channel_string, email, password, varargin)
%this function goes and fetches data from the IRIS and other DMC using irisFetch. 
%USAGE: [Traces] = fetchData(E,S,pre, post, phase, tag, taper_fraction, sample_rate, ...
%    channel_string, email, password)
%INPUT: pre and post event time in are in seconds
%       phase is either 'P' or 'S'
%       tag is a string that will be used to name the output file
%       taper_fraction is the fraction of the trace to taper at the beginning and end
%       sample_rate is the desired sample rate in Hz
%       channel_string is a string that will be used to select the channel
%       email and password are for IRIS DMC access
%       varargin is option for starting event number and datacenter list ({dacacenter1, datacenter2, ...},)
%       for datacenter {name,baselink}(e.g. {'IRIS','http://service.iris.edu/irisws/timeseries/1/query'})
% where E and S are event and station structures as would be produced by
% irisFetch.
% The function fetches, for each event, all traces recorded by stations in
% the station list, if any.
% The data is saved in a .mat file named tag_eventinfo_### where ### is sequential
% numbers.
% The user can optionally specify a file name root with an optional iuput
% The optional fourth argument is for re-starting the process at a certain
% event number if it was originally interrupted somehow. For example, if
% you're downloading 350 events, but your internet went out and the last
% event downloaded was 234 you can specify eventStart as 235 so it resumes
% downloading the missing events.
% This version does multiple traces for each station
% This version loops over data center as well. The station structure needs
% to have a field called DataCenter (this is not irisFetch standard). Right
% now that field can have values IRIS, GEOFON or ORFEUS and AUPASS

% Attention!!!!
% email and password are not in used now!!!
% !!!!!

myTrace = [];

% now double loop through events and stations

%seconds in fractional days (because that's what MATLAB uses for date/times)
pre  = pre/(24*60*60); 
post = post/(24*60*60);

kestart = 1;
datacenter_list={};

for i=1:numel(varargin)
    switch class(varargin{i})
        case 'double'
            kestart = varargin{i};
        case 'struct'
            datacenter_list = varargin{i};
    end
end

Sall=S; %these are the stations from all the datacenters

for ke=kestart:length(E)
    eventData = E(ke);
    eventtime_yyMMddHHmmss = datetime(eventData.PreferredTime);
    eventtime_yyMMddHHmmss.Format = 'yyMMddHHmmss';
    fname     = sprintf('%s_%s',tag,eventtime_yyMMddHHmmss);
    if isfile(strcat(fname,'.mat'))
        disp(['event ' fname ' already downloaded, skip'])
        continue
    end

    disp([ num2str(ke) ' of ' num2str(length(E)) ]);
    
    Traces=[];

    %loop over datacenters
    for kdc=1:length(datacenter_list)

        %loop over stations
        datacenter = datacenter_list(kdc);
        if datacenter.used == 0
            disp(['skipped download data from ' datacenter.name]);
            continue;
        else
            disp(['download data from ' datacenter.name]);
        end
        S = Sall(strcmpi({Sall.DataCenter},datacenter.name)); %these are the stations from this datacenter
        for ks=1:length(S)
            %first check to see if this station was active when the earthquake
            %happened
            sDate=datenum(S(ks).StartDate);
            if isempty(S(ks).EndDate); S(ks).EndDate='2500-01-01 00:00:00.000'; end
            eDate=datenum(S(ks).EndDate);
            eqDate=datenum(E(ke).PreferredTime);
            if ~(eqDate>sDate && eqDate<eDate) %if eq not between eDate and sDate
                %disp('skipped one')
                continue
            end
            %calculate Delta
            [D,Az]=distance(E(ke).PreferredLatitude,E(ke).PreferredLongitude,S(ks).Latitude,S(ks).Longitude);
            urlstr=sprintf('http://service.iris.edu/irisws/traveltime/1/query?distdeg=%f&evdepth=%f&noheader=true&mintimeonly=true',D,E(ke).PreferredDepth);
            urlstr2=sprintf('http://service.iris.edu/irisws/traveltime/1/query?distdeg=%f&evdepth=%f&mintimeonly=true',D,E(ke).PreferredDepth);
            tStr=urlread(urlstr);

            tCell=textscan(tStr,'%f %f %s %f %f %f %f %f %s %s');

            phases=tCell{3};
            times=tCell{4};
            
            pTime=times(strcmpi(phases,'P')); %this is time to P arrival.
            sTime=times(strcmpi(phases,'S'));
            
            pTime=min(pTime);
            sTime=min(sTime);
            
            if phase == 'P'
                phase_time = pTime;
            elseif phase == 'S'
                phase_time = sTime;     
            end
            if isempty(phase_time); phase_time=times(1); end
            originTimeStr  = E(ke).PreferredTime;
            originTimeNum  = datenum(originTimeStr);
            
            startTime = originTimeNum+phase_time/(24*60*60)-pre;
            startTime = datestr(startTime,'yyyy-mm-dd HH:MM:SS.FFF');
            
            endTime   = originTimeNum+phase_time/(24*60*60)+post;
            endTime   = datestr(endTime,'yyyy-mm-dd HH:MM:SS.FFF');

            station_starttime = datenum(S(ks).StartDate);
            station_endtime = datenum(S(ks).EndDate);
            if station_starttime > datenum(startTime) || station_endtime < datenum(endTime)
                fprintf('Event %s has no record on station %s\n',fname,S(ks).StationCode)
                continue
            end
            
            %fprintf('Trying event #%.0f station %s\n',ke,S(ks).StationCode)
            try
                myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                S(ks).StationCode,'*',channel_string,startTime,endTime,...
                datacenter.baselink);
                %if strcmp(email, '')
                
                %    myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                %            S(ks).StationCode,'*',channel_string,startTime, endTime);
                %else
                %    myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                %            S(ks).StationCode,'*',channel_string,startTime, endTime, ...
                %            { email, password });
                %end  
            catch ME
                continue         
            end
            %add phase names and times and add to the whole
            if ~isempty(myTrace)
                try
                    [myTrace.phaseNames]=deal({'P','S'});
                catch ME
                    keyboard
                end
                %note these times are relative to the start of the trace:
                if phase == 'P'
                    [myTrace.phaseTimes]=deal([0 (sTime-pTime)]+pre*(24*60*60)); %pre back to seconds
                elseif phase == 'S'
                    [myTrace.phaseTimes]=deal([ (pTime - sTime) 0]+pre*(24*60*60)); %pre back to seconds
                end
                %I know, this is growing inside  a loop.
                try
                    myTrace=getresp_fsdn(myTrace,datacenter);
                catch ME
                    disp(['getresp_fsdn failed with' ME.message]);
                end
                Traces = [Traces myTrace];
                
                assignin('base','T',Traces)
            end
            myTrace = [];
        end

    end
    %combine all saperate trace
    Traces = connect_trace_data(Traces,2,pre,post);
    %skip to next event if this event had no traces whatsover
    if isempty(Traces)
        disp([ 'event ' num2str(ke) ' is empty' ]);
        continue
    end
    
    % Get rid of empty structures:
    data = {Traces.data};
    tf_empty  = cellfun('isempty',data);
    Traces    = Traces(~tf_empty);
    for q = 1:length(Traces)
        %demean/detrend
        Traces(q).data = detrend(Traces(q).data);
        %taper
        Traces(q).data = tukeywin(length(Traces(q).data), taper_fraction).*Traces(q).data;
    end
        
    Traces   = wfResample(Traces, sample_rate);
    Traces   = wfRemInstResp_sacpz(Traces);
    
    save(fname,'Traces','eventData')
    
end
end

function sanity_check(Traces,netstacha)
%sanity_check checks the data in Traces
    for j = 1:length(Traces)
        if length(Traces(j).data) ~= Traces(j).sampleCount
            warning(['sampleCount is not equal to length of data ', netstacha{j}])
        end
        tem = round((Traces(j).endTime - Traces(j).startTime)*24*60*60*Traces(j).sampleRate+1);
        if Traces(j).sampleCount ~= tem
            warning(['sampleCount is not equal to duration ', netstacha{j} ' ' num2str(Traces(j).sampleCount) ' ' num2str(tem)])
        end
    end
end

function S_trace = connect_trace_data(Traces,flag,pre,post)
%connect_trace_data connects the seismic trace from one station if there is more than one in Traces
%   connect_trace_data(Traces,flag)
%   flag = 1: linear interplate
%   flag = 2: spling interplate
%   pre and post are in matlab datenum format

    S_trace = [];
    netstacha = {};
    pre = pre*24*60*60;
    post = post*24*60*60;
    
    for j = 1:length(Traces)
        netstacha(j) = {[Traces(j).network '_' Traces(j).station '_' Traces(j).channel]};
    end
    sanity_check(Traces,netstacha);
    unetstacha = unique(netstacha);
    for j = 1:length(unetstacha)
    % for each unique station
        trace_tem_list = Traces(strcmp(netstacha,unetstacha(j)));
        if length(trace_tem_list) == 1
            %only one trace
            sampleRate = trace_tem_list.sampleRate;
            total_count = (pre+post)*sampleRate;
            if trace_tem_list.sampleCount ~= total_count
                warning(['The sample count is not equal to the request time ',unetstacha{j}])
            else
                S_trace = [S_trace trace_tem_list];
            end
        else
            if length(unique([trace_tem_list.sampleRate])) ~= 1
                warning(['The sample rate of the same station is not the same ',unetstacha{j}])
                continue
            else
                sampleRate = unique([trace_tem_list.sampleRate]);
            end
            t_start=trace_tem_list(1).startTime;
            t_duration=[];
            data_tem=[];
            t_total_duration = 0:round((trace_tem_list(end).endTime-t_start)*24*60*60*sampleRate);
            total_count = (pre+post)*sampleRate;

            for k = 1:length(trace_tem_list)
            %for k = 1 : 5
                t_start_tem = round((trace_tem_list(k).startTime-t_start)*24*60*60*sampleRate);
                t_end_tem = round((trace_tem_list(k).endTime-t_start)*24*60*60*sampleRate); 
                t_duration_tem=t_start_tem:t_end_tem;
                t_duration=[t_duration t_duration_tem];
                data_tem=[data_tem trace_tem_list(k).data'];
            end
            [tt_duration,i_tt_duration,~]=unique(t_duration);
            data_tem_f = data_tem(i_tt_duration);
            if length(tt_duration) ~= length(t_duration)
                % check if there is any overlaped record
                warning(['The potential sample overlap exist ',unetstacha{j}])
            end
            if flag == 1
                data_final = interp1(tt_duration,data_tem_f,t_total_duration,'linear');
            elseif flag == 2
                data_final = interp1(tt_duration,data_tem_f,t_total_duration,'spline');
            end
            trace_tem = trace_tem_list(1);
            trace_tem.data = data_final';
            trace_tem.startTime = t_start;
            trace_tem.endTime = trace_tem_list(end).endTime;
            trace_tem.sampleCount = length(data_final);
            if trace_tem.sampleCount <= total_count
                warning(['The sample count is less than to the request time ',unetstacha{j}])
            elseif trace_tem.sampleCount > total_count
                warning(['The sample count is more than to the request time, cut off the later data point ',unetstacha{j}])
                trace_tem.data = trace_tem.data(1:total_count);
                trace_tem.sampleCount = total_count;
                trace_tem.endTim=trace_tem.startTime+(total_count-1)/sampleRate/24/60/60;
                S_trace = [S_trace trace_tem];
            else
                % sample count is equal to the request time
                S_trace = [S_trace trace_tem];
            end
        end
    end
end
