function [S] = findStations(srchParams,datacenter_list)

%function [] = findStations(latLonbox,sTime,eTime)
%latLonbox = [minLat maxLat minLon maxLon]
%sTime, eTime format:'2010-01-01 00:00:00'
%get stations

%unpack search parameters and use defaults as needed.
sTime=srchParams.sTime;
eTime=srchParams.eTime;
latLonBox=srchParams.latLonBox;
% 'detail' parameter
if isfield(srchParams,'detl')
    detl=srchParams.detl;
else
    detl='';
end
% 'network' parameter
if isfield(srchParams,'net')
    net=srchParams.net;
else
    net='';
end
% 'station' parameter
if isfield(srchParams,'sta')
    sta=srchParams.sta;
else
    sta='';
end
% 'location' parameter
if isfield(srchParams,'loc')
    loc=srchParams.loc;
else
    loc='';
end
% 'channel' parameter
if isfield(srchParams,'chan')
    chan=srchParams.chan;
else
    chan='BHZ,HHZ';
end


for i=1:length(datacenter_list)
    if strcmp(datacenter_list(i).name,'IRIS')
        Siris=irisFetch.Stations(detl,net,sta,loc,chan,'boxcoordinates'...
            ,latLonBox,'startTime',sTime,'endTime',eTime,'StartAfter','1950-01-01');
        [Siris.DataCenter]=deal('IRIS');
    elseif datacenter_list(i).used == true
        stationurl = strcat(datacenter_list(i).baselink,'/fdsnws/station/1/');
        temS = irisFetch.Stations(detl,net,sta,loc,chan,'boxcoordinates'...
            ,latLonBox,'startTime',sTime,'endTime',eTime,'StartAfter','1950-01-01','BASEURL',stationurl);
        [temS.DataCenter]=deal(datacenter_list(i).name);
        if isempty(temS)
            disp(['no data from datacenter: ', datacenter_list(i).name])
            continue
        end 
        if exist('allothers','var')
            allothers = [allothers, temS];
        else
            allothers = temS;
        end
    else
        %skip this datacenter
        continue
    end
end
%note: the startafter parameter is to avoid picking up synthetic stations
%that have startTimes of 1900-01-01.

%SS=[Siris, Sgeofon, Sorfeus, Sauspass];
if exist('allothers','var')
    SS=[Siris, allothers];
else
    SS=Siris;
end
cde={SS.StationCode};
[ucd ix]= unique(cde);

S=SS(ix);
