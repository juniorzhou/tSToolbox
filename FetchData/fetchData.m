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
% now that field can have values IRIS, GEOFON or ORFEUS. I will implement
% AUSPASS later.

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
    
    disp([ num2str(ke) ' of ' num2str(length(E)) ]);
    
    Traces=[];

    %loop over datacenters
    for kdc=1:length(datacenter_list)

        %loop over stations
        datacenter = datacenter_list(kdc);
        if datacenter.used == 0
            disp(['skipped download data from ' datacenter.name]);
            continue;
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
                Traces = [Traces myTrace];
                
                assignin('base','T',Traces)
            end
            myTrace = [];
        end

    end

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
        
    %Traces   = wfResample(Traces, sample_rate);
    %Traces   = wfRemInstResp(Traces);
    
    eventData = E(ke);
    eventtime_yyMMddHHmmss = datetime(eventData.PreferredTime);
    eventtime_yyMMddHHmmss.Format = 'yyMMddHHmmss';
    %fname     = sprintf('%s_%s_%03.0f',tag,eventtime_yyMMddHHmmss,ke);
    fname     = sprintf('%s_%s',tag,eventtime_yyMMddHHmmss);
    save(fname,'Traces','eventData')
    
end