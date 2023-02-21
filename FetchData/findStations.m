function [S] = findStations(srchParams)

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

%S=irisFetch.Stations(detl,net,sta,loc,chan,'boxcoordinates'...
%    ,latLonBox,'startTime',sTime,'endTime',eTime,'StartAfter','1950-01-01');

Sgeofon=irisFetch.Stations(detl,net,sta,loc,chan,'boxcoordinates'...
    ,latLonBox,'startTime',sTime,'endTime',eTime,'StartAfter','1950-01-01','BASEURL','https://geofon.gfz-potsdam.de/fdsnws/station/1/');
[Sgeofon.DataCenter]=deal('GEOFON');

Sorfeus=irisFetch.Stations(detl,net,sta,loc,chan,'boxcoordinates'...
    ,latLonBox,'startTime',sTime,'endTime',eTime,'StartAfter','1950-01-01','BASEURL','http://www.orfeus-eu.org/fdsnws/station/1/');
[Sorfeus.DataCenter]=deal('ORFEUS');

Siris=irisFetch.Stations(detl,net,sta,loc,chan,'boxcoordinates'...
    ,latLonBox,'startTime',sTime,'endTime',eTime,'StartAfter','1950-01-01');
[Siris.DataCenter]=deal('IRIS');

%note: the startafter parameter is to avoid picking up synthetic stations
%that have startTimes of 1900-01-01.

SS=[Siris, Sgeofon, Sorfeus];
cde={SS.StationCode};
[ucd ix]= unique(cde);

S=SS(ix);
