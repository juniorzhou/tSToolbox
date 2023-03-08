function [trace_resp] = getresp_fsdn(Trace,datacenter)
%this function fetches resp from the datacenter use fsdn request with level resp. 
%USAGE: [Trace] = fetchData(Trace,datacenter)
%INPUT: Trace: a structure array with the following fields
%           network: network code
%           station: station code
%           location: location code
%           channel: channel code
%           startTime: start time in matlab datenum format
%           endTime: end time in matlab datenum format
%           sacpz: a structure with the following fields
%               units: units of the Response
%               constant: constant of the Response
%               poles: poles of the Response
%               zeros: zeros of the Response
%       datacenter: a structure with the following fields
%           baselink: the base link to the data center
%           example: https://auspass.edu.au/
%OUTPUT: trace_resp is the same as the input trace but with the sacpz field
%example of the http request for xmlfile
%       https://auspass.edu.au/fdsnws/station/1/query?net=6F&sta=BL01&
%       cha=BHZ&start=2009-09-29T17:50:00&end=2009-09-29T18:50:00&level=resp

    trace_resp=Trace;
    for i=1:length(Trace)
        %get the response for each trace
        url=datacenter.baselink;
        network=Trace(i).network;
        station=Trace(i).station;
        location=Trace(i).location;
        channel=Trace(i).channel;
        starttime_m=Trace(i).startTime;
        endtime_m=Trace(i).endTime;
        starttime=datetime(starttime_m,'ConvertFrom','datenum');
        endtime=datetime(endtime_m,'ConvertFrom','datenum');
        starttime_str=datestr(starttime,'yyyy-mm-ddTHH:MM:SS');
        endtime_str=datestr(endtime,'yyyy-mm-ddTHH:MM:SS');
        url=[url '/fdsnws/station/1/query?net=' network '&sta=' station '&cha=' channel...
         '&start=' starttime_str '&end=' endtime_str  '&location=' location '&level=resp'];
        %get the response
        try
            xmlresp_str=webread(url);
        catch
            warning('could not get response for station %s_%s\n%s', ...
                network,station,url);
            trace_resp(i).sacpz.used=0;
            continue
        end
        if ~strcmp(xmlresp_str(3:5),'xml')
            warning("wrong return value for link %s, try again",url)
            try
                xmlresp_str=webread(url);
            catch
                warning('could not get response for station %s_%s\n%s', ...
                    network,station,url);
                trace_resp(i).sacpz.used=0;
                continue
            end
        end
        if ~strcmp(xmlresp_str(3:5),'xml')
            warning("wrong return value for link %s, skip",url)
            continue
        end    
        resp=parse_xmlresp(xmlresp_str,url);
        trace_resp(i)=Trace(i);
        if ~strcmp(resp.units, trace_resp(i).sensitivityUnits)
            warning(['response units are not same as the Trace, %s vs %s, this will likely ' ...
                'cause problems station %s_%s'],resp.units,trace_resp(i).sensitivityUnits,network,station);
        end
        trace_resp(i).sacpz.units=resp.units;
        trace_resp(i).sacpz.constant=resp.constant;
        trace_resp(i).sacpz.poles=resp.poles;
        trace_resp(i).sacpz.zeros=resp.zeros;
        trace_resp(i).sacpz.used=1;
    end
end


function resp = parse_xmlresp(xmlresp_str,url)
    resp={};
    X=xml2struct_own(xmlresp_str,url);
    %subset to network
    N=X.Children(strcmp({X.Children.Name},'Network'));
    %subset to station
    S=N.Children(strcmp({N.Children.Name},'Station'));
    %subset to channel
    C=S.Children(strcmp({S.Children.Name},'Channel'));
    %subset to response
    R=C.Children(strcmp({C.Children.Name},'Response'));
    %subset ot instrument Sensitivity
    IR=R.Children(strcmp({R.Children.Name},'InstrumentSensitivity'));
    unitir=IR.Children(strcmp({IR.Children.Name},'InputUnits'));
    unitir_str=unitir.Children(strcmp({unitir.Children.Name},'Name')).Children.Data;
    resp.units=unitir_str;
    %subset to appropriate Stage
    %first find the stage that has the poles and zeros
    Rc=R.Children;
    for k=1:length(Rc)
        try
            tf(k)=any(strcmp({Rc(k).Children.Name},'PolesZeros'));
        catch
            tf(k)=false;
        end
    end
    R1=R.Children(tf); 
    %check, just in case
    if sum(tf)>1 
        warning('more than one polezeros detected, use first one default');
        R1=R1(1);
    end
    %now subset
    %subset to poles and zeros
    PZ=R1.Children(strcmp({R1.Children.Name},'PolesZeros'));
    
    zix=strcmp({PZ.Children.Name},'Zero');
    pix=strcmp({PZ.Children.Name},'Pole');
    
    %extract values for each pole
    PolStruc=PZ.Children(pix);
    for k=1:length(PolStruc)
        P=PolStruc(k);
        Pr=P.Children(strcmp({P.Children.Name},'Real'));
        rl=str2double(Pr.Children.Data);
        Pi=P.Children(strcmp({P.Children.Name},'Imaginary'));
        im=str2double(Pi.Children.Data);
        Poles(k)=rl+1i*im;
    end
    
    %extract values for each zero
    ZerStruc=PZ.Children(zix);
    for k=1:length(ZerStruc)
        Z=ZerStruc(k);
        Zr=Z.Children(strcmp({Z.Children.Name},'Real'));
        rl=str2double(Zr.Children.Data);
        Zi=Z.Children(strcmp({Z.Children.Name},'Imaginary'));
        im=str2double(Zi.Children.Data);
        Zeros(k)=rl+1i*im;
    end
    
    %other things
    nZeros=sum(zix); 
    nPoles=sum(pix);
    
    %normalization factor
    NF=PZ.Children(strcmp({PZ.Children.Name},'NormalizationFactor'));
    NormFactor=str2double(NF.Children.Data);
    resp.poles=Poles;
    resp.zeros=Zeros;
    resp.constant=NormFactor;
end

function theStruct = xml2struct_own(xmlString,url)
% xml2struct_own Convert XML string to a MATLAB structure.
    try
        % The following avoids the need for file I/O:
        inputObject = java.io.StringBufferInputStream(xmlString);  % or: org.xml.sax.InputSource(java.io.StringReader(xmlString))
        try
            % Parse the input data directly using xmlread's core functionality
            parserFactory = javaMethod('newInstance','javax.xml.parsers.DocumentBuilderFactory');
            p = javaMethod('newDocumentBuilder',parserFactory);
            xmlTreeObject = p.parse(inputObject);
        catch
            % Use xmlread's semi-documented inputObject input feature
            xmlTreeObject = xmlread(inputObject);
        end
    catch
        % Fallback to standard xmlread usage, using a temporary XML file:
        % Store the XML data in a temp *.xml file
        filename = append('tem_',string(datetime("now","Format","mmssSSS")),'.xml');
        fid = fopen(filename,'Wt');
        fwrite(fid,xmlString);
        fclose(fid);
        % Read the file into an XML model object
        xmlTreeObject = xmlread(filename);
        % Delete the temp file
        delete(filename);
    end
    try
        theStruct = parseChildNodes(xmlTreeObject);
    catch
        error(message('bioinfo:xml2struct:XMLParseError'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodeStruct = makeStructFromNode(theNode)

    nodeStruct = struct('Name',char(theNode.getNodeName),...
        'Attributes',parseAttributes(theNode),'Data','',...
        'Children',parseChildNodes(theNode));
    
    if any(strcmp(methods(theNode),'getData'))
       nodeStruct.Data = char(theNode.getData); 
    else
        nodeStruct.Data = '';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function attributes = parseAttributes(theNode)
    % Create attributes struct
    attributes = [];
    if theNode.hasAttributes
        theAttributes = theNode.getAttributes;
        numAttributes = theAttributes.getLength;
        allocCell = cell(1,numAttributes);
        attributes = struct('Name',allocCell,'Value',allocCell);
        for count = 1:numAttributes
            attrib = theAttributes.item(count-1);
            attributes(count).Name = char(attrib.getName);
            attributes(count).Value = char(attrib.getValue);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function children = parseChildNodes(theNode)
% Recurse over node children
    children = [];
    if theNode.hasChildNodes
        childNodes = theNode.getChildNodes;
        numChildNodes = childNodes.getLength;
        allocCell = cell(1,numChildNodes);
        children = struct('Name',allocCell,'Attributes',allocCell,...
                                        'Data',allocCell,'Children',allocCell);
        for count = 1:numChildNodes
            theChild = childNodes.item(count-1);
            children(count) = makeStructFromNode(theChild);
        end
    end
end