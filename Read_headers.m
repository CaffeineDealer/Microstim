function acquisition_events=Read_headers(name)
%% Info about the headers can be found at:
% https://rippleneuro.s3-us-west-2.amazonaws.com/sites/5817b02bbc9d752748a1409d/assets/58559c85bc9d756e810d825d/NEVspec2_2_v07.pdf

% Read the nsX file. Contains electrode information and data.
% It saves a binary file with the raw signal per electrode.
% The headers are saved on a separate file (headers.mat)

tic
path = 'C:\research\data\saul';
data_file = [name,'.nev.ns5'];

% Get the Headers: Header, extendedHeader
run = 1;

signals_and_headers = struct;
more_runs_in_this_file = 'true';  %If you press pause, while recording, it will add headers again after each run, so the whole reading process needs to be repeated. See the diagram on the pdf


if ~exist([path '/' data_file(1:end-4)])
    mkdir([path '/' data_file(1:end-4)])
end

while more_runs_in_this_file

    fid = fopen([path '/' data_file],'rb');
    signals_and_headers(run).header.FileTypeID               = fread(fid, [1,8],  'uint8=>char');
    signals_and_headers(run).header.FileSpec                 = fread(fid, [1,2],  'uchar'); % Major.Minor revision version (e.g. 2.2)
    signals_and_headers(run).header.BytesInHeaders           = fread(fid, [1,1],  'uint32');
    signals_and_headers(run).header.Label                    = fread(fid, [1,16], 'uint8=>char');
    signals_and_headers(run).header.Comments                 = fread(fid, [1,200],'uint8=>char');
    signals_and_headers(run).header.ApplicationToCreateFile  = fread(fid, [1,52], 'uint8=>char');
    signals_and_headers(run).header.ProcessorTimestamp       = fread(fid, [1,1],  'uint32');
    signals_and_headers(run).header.Period                   = fread(fid, [1,1],  'uint32');
    
    fs_binary_file_index = ftell(fid); % This is used so we can 
    
    signals_and_headers(run).header.TimeResolutionTimeStamps = fread(fid, [1,1],  'uint32');
    signals_and_headers(run).header.TimeOrigin               = fread(fid, [1,8],  'uint16'); % When the data was recorded: Year,Month,DayOfWeek,Day,Hour,Minute,Second,Millisecond
    signals_and_headers(run).header.ChannelCount             = fread(fid, [1,1],  'uint32'); %

    signals_and_headers(run).extendedheader = struct;

    % Get Extended Header
    for ielectrode = 1:signals_and_headers(run).header.ChannelCount
        signals_and_headers(run).extendedheader(ielectrode).Type                = fread(fid, [1,2], 'uint8=>char'); %% CC = Continuous Channels
        signals_and_headers(run).extendedheader(ielectrode).ElectrodeID         = fread(fid, [1,1], 'uint16');
        signals_and_headers(run).extendedheader(ielectrode).ElectrodeLabel      = fread(fid, [1,16],'uint8=>char');
        signals_and_headers(run).extendedheader(ielectrode).FrontEndID          = fread(fid, [1,1], 'uchar');
        signals_and_headers(run).extendedheader(ielectrode).FrontEndPin         = fread(fid, [1,1], 'uchar');
        signals_and_headers(run).extendedheader(ielectrode).MinDigitalValue     = fread(fid, [1,1], 'int16');
        signals_and_headers(run).extendedheader(ielectrode).MaxDigitalValue     = fread(fid, [1,1], 'int16');
        signals_and_headers(run).extendedheader(ielectrode).MinAnalogValue      = fread(fid, [1,1], 'int16');
        signals_and_headers(run).extendedheader(ielectrode).MaxAnalogValue      = fread(fid, [1,1], 'int16');
        signals_and_headers(run).extendedheader(ielectrode).Units               = fread(fid, [1,16],'uint8=>char'); 
        signals_and_headers(run).extendedheader(ielectrode).HighPassCornerFreq  = fread(fid, [1,1], 'uint32'); %in mHz
        signals_and_headers(run).extendedheader(ielectrode).HighPassFilterOrder = fread(fid, [1,1], 'uint32'); %0 = NONE
        signals_and_headers(run).extendedheader(ielectrode).HighPassFilterType  = fread(fid, [1,1], 'uint16'); %0=None,1=Butterworth,2=Chebyshev
        signals_and_headers(run).extendedheader(ielectrode).LowPassCornerFreq   = fread(fid, [1,1], 'uint32'); %in mHz
        signals_and_headers(run).extendedheader(ielectrode).LowPassFilterOrder  = fread(fid, [1,1], 'uint32'); %0 = NONE
        signals_and_headers(run).extendedheader(ielectrode).LowPassFilterType   = fread(fid, [1,1], 'uint16'); %0=None,1=Butterworth,2=Chebyshev
    end

%     position = ftell(fid)
% stop
    % Get the data

    % This section of the NSx file contains an open-ended number of packets consisting of a
    % header, a timestamp, the total number of sampled data points, and a variable number of
    % data points.
    % The timestamp reported in each data packet indicates the sampling time for the first data
    % point stored in that packet. If a pause occurs during data storage (see Grapevine User
    % Manual), a new data packet will be created when data storage resumes at the end of the
    % pause. The timestamp field is followed by a field specifying the number of data points
    % stored in the packet. Data points correspond to a single point in time and the data points
    % are in order of increasing time. Each data point consist of samples from one or more
    % channels (see Diagram below table). 

    signals_and_headers(run).Data.Header             = fread(fid, [1,1], 'int8');
    signals_and_headers(run).Data.Timestamp          = fread(fid, [1,1], 'uint32');
    signals_and_headers(run).Data.NumberOfDataPoints = fread(fid, [1,1], 'uint32');
%     position = ftell(fid)
% stop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signals_and_headers(run).Data.NumberOfDataPoints = 6000000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    data_binary_file_index = ftell(fid);

    signals_and_headers(run).Data.F = fread(fid, [signals_and_headers(run).header.ChannelCount, signals_and_headers(run).Data.NumberOfDataPoints], 'int16=>int16'); % This saves all the electrodes on a single matrix. Consider saving per electrode.
    signals_and_headers(run).Data.F = signals_and_headers(run).Data.F.*(signals_and_headers(run).extendedheader(1).MaxAnalogValue/signals_and_headers(run).extendedheader(1).MaxDigitalValue); % Convert to uV - This assumes that all 
                                                                                                                                                                % electrodes have the same MaxDigitalValue and MaxAnalogValue
    
                                                                                                                                                                
    pause_pressed = fread(fid, [1,2], 'uint8=>char'); % If you press pause, while recording it will add headers again after each run, so the whole process needs to be repeated. See the diagram on the pdf
    more_runs_in_this_file = strcmp(pause_pressed,'CC');
    run = run+1;
end

fclose(fid);

header = signals_and_headers.header;
extendedheader = signals_and_headers.extendedheader;
NumberOfDataPoints = signals_and_headers.Data.NumberOfDataPoints;
ElectrodeIDs = [extendedheader.ElectrodeID];

save([path '/' data_file(1:end-4) '/headers.mat'],'header','extendedheader','NumberOfDataPoints','ElectrodeIDs', 'fs_binary_file_index', 'data_binary_file_index')



for ielectrode = 1:signals_and_headers.header.ChannelCount
        if ~exist([path '/' data_file(1:end-4)])
            mkdir([path '/' data_file(1:end-4)])
        end
        
        ftemp = fopen([path '/' data_file(1:end-4) '/raw_elec' num2str(signals_and_headers.extendedheader(ielectrode).ElectrodeID) '.bin'], 'W');  % uppercase W
        fwrite(ftemp, signals_and_headers(1).Data.F(ielectrode,:),'int16'); % I have signals_and_headers(1)... here since the pause button was not pressed. GENERALIZE
        fclose(ftemp);    
end           


clear data_file fid isegment ielectrode ans more_runs_in_this_file run pause_pressed ielectrode ftemp path single_electrode_signal



% path = 'F:/Ripple Files/m120baseline';
% data_file = 'raw_elec1.bin'; % Reading takes ~
% 
% fid = fopen([path '/' data_file],'rb');
% aa = fread(fid, [1, 18000420], 'int16=>int16');
% fclose(fid);



toc

%% Read the NEV file (Contains events - Spikes - Stimulation Waveforms etc.)
% path = 'F:/Ripple Files/';
path = 'C:\research\data\saul\';
% file = 'F:/Ripple Files/m120baseline.nev';
% file = 'F:/Ripple Files/m020b.nev';
file = [name,'.nev.nev'];   

info_file = [path file];     


fid = fopen(info_file,'rb');
header.FileTypeID                 = fread(fid, [1,8],  'uint8=>char');
header.FileSpec                   = fread(fid, [1,2],  'uchar'); % Major.Minor revision version (e.g. 2.2)
header.AdditionalFlags            = fread(fid, [1,1],  'uint16');
header.BytesInHeaders             = fread(fid, [1,1],  'uint32');
header.BytesinDataPackets         = fread(fid, [1,1],  'uint32');
header.TimeResolutionOfTimestamps = fread(fid, [1,1],  'uint32');
header.TimeResolutionOfSamples    = fread(fid, [1,1],  'uint32');
header.TimeOrigin                 = fread(fid, [1,8],  'int16'); % When the data was recorded: Year,Month,DayOfWeek,Day,Hour,Minute,Second,Millisecond
header.ApplicationToCreateFile    = fread(fid, [1,32], 'uint8=>char');
header.CommentField               = fread(fid, [1,200],'uint8=>char'); 
header.Reserved                   = fread(fid, [1,52], 'uint8'); %
header.ProcessorTimestamp         = fread(fid, [1,1],  'uint32'); %
header.NumberOfExtendedHeaders    = fread(fid, [1,1],  'uint32'); %


clear extendedheader

% Get Extended Header
for iExtendedHeader = 1:header.NumberOfExtendedHeaders
    extendedheader(iExtendedHeader).PacketID = fread(fid, [1,8], 'uint8=>char'); %% CC = Continuous Channels
    
    switch extendedheader(iExtendedHeader).PacketID
        
        case 'NEUEVWAV'
            extendedheader(iExtendedHeader).ElectrodeID                 = fread(fid, [1,1], 'uint16');
            extendedheader(iExtendedHeader).FrontEndID                  = fread(fid, [1,1], 'uchar');
            extendedheader(iExtendedHeader).FrontEndConnectorPin        = fread(fid, [1,1], 'uchar');
            extendedheader(iExtendedHeader).NeuralAmpDigitizationFactor = fread(fid, [1,1], 'uint16');
            extendedheader(iExtendedHeader).EnergyThreshold             = fread(fid, [1,1], 'uint16');
            extendedheader(iExtendedHeader).HighThreshold               = fread(fid, [1,1], 'int16');
            extendedheader(iExtendedHeader).LowThreshold                = fread(fid, [1,1], 'int16');
            extendedheader(iExtendedHeader).NumberOfSortedUnits         = fread(fid, [1,1], 'uchar');
            extendedheader(iExtendedHeader).BytesPerSample              = fread(fid, [1,1], 'uchar');
            extendedheader(iExtendedHeader).StimAmplDigitizationFactor  = fread(fid, [1,1], 'float');
            extendedheader(iExtendedHeader).Reserved                    = fread(fid, [1,6], 'int8'); %CHECK IF THIS IS ZERO
            bytesPerSample = extendedheader(iExtendedHeader).BytesPerSample; % This will be used in the packets later

        case 'NEUEVFLT'
            extendedheader(iExtendedHeader).ElectrodeID             = fread(fid, [1,1], 'uint16');
            extendedheader(iExtendedHeader).HighPassCornerFrequency = fread(fid, [1,1], 'uint32');
            extendedheader(iExtendedHeader).HighPassFilterOrder     = fread(fid, [1,1], 'uint32');
            extendedheader(iExtendedHeader).HighPassFilterType      = fread(fid, [1,1], 'uint16');
            extendedheader(iExtendedHeader).LowPassCornerFrequency  = fread(fid, [1,1], 'uint32');
            extendedheader(iExtendedHeader).LowPassFilterOrder      = fread(fid, [1,1], 'uint32');
            extendedheader(iExtendedHeader).LowPassFilterType       = fread(fid, [1,1], 'uint16');
            extendedheader(iExtendedHeader).Reserved                = fread(fid, [1,2], 'int8');  %CHECK IF THIS IS ZERO
            
              
        case 'NEUEVLBL'
            extendedheader(iExtendedHeader).ElectrodeID = fread(fid, [1,1], 'uint16');
            extendedheader(iExtendedHeader).Label       = fread(fid, [1,16],'uint8=>char');
            extendedheader(iExtendedHeader).Reserved    = fread(fid, [1,6], 'int8');  %CHECK IF THIS IS ZERO

            
        case 'DIGLABEL'
            extendedheader(iExtendedHeader).Label       = fread(fid, [1,16], 'uint8=>char');
            extendedheader(iExtendedHeader).Mode        = fread(fid, [1,1], 'int8'); % 0:serial, 1:parallel
            extendedheader(iExtendedHeader).Reserved    = fread(fid, [1,7], 'int8'); % CHECK IF THIS IS ZERO
    end
    
end

%

% CREATE STIMULATION EVENTS - SPIKING EVENTS - PARALLEL PORT EVENTS
% The struct packets contains everything that was captured from the
% acquisition system. It also saves the waveforms of each spike
% (redundant). 
% The struct: events, contains only what is really needed, in Brainstorm
% friendly format

clear packets events
packets = struct;

ipacket = 1;
are_there_more_packets = true;
tic
while are_there_more_packets %ipacket<10000
    packets(ipacket).Timestamp = fread(fid, [1,1], 'uint32');
    packets(ipacket).PacketID  = fread(fid, [1,1], 'uint16');
    
    if isempty(packets(ipacket).Timestamp)
        warning('No more packets')
        break
    end
    
    if packets(ipacket).PacketID ==0 % These are events from the parallel port
        
        packets(ipacket).PacketInsertionReason = fread(fid, [1,1], 'uchar');
        packets(ipacket).Reserved              = fread(fid, [1,1], 'uchar');
        packets(ipacket).ParallelInput         = fread(fid, [1,1], 'uint16');
        packets(ipacket).SAMinput1             = int16(fread(fid, [1,1], '*int16'));
        packets(ipacket).SAMinput2             = int16(fread(fid, [1,1], '*int16'));
        packets(ipacket).SAMinput3             = int16(fread(fid, [1,1], '*int16'));
        packets(ipacket).SAMinput4             = int16(fread(fid, [1,1], '*int16'));
        packets(ipacket).ReservedEnding        = fread(fid, [1,header.BytesinDataPackets-18],'*char');
        
        
    elseif packets(ipacket).PacketID >0 && packets(ipacket).PacketID <513  % These packets indicate spikes
        packets(ipacket).UnitClassificationNumber = fread(fid, [1,1], 'uchar');
        packets(ipacket).Reserved                 = fread(fid, [1,1], 'uchar');
        
        if bytesPerSample == 2
            precision = '*int16';
            packets(ipacket).Waveform = fread(fid, [1,header.BytesinDataPackets/bytesPerSample - 4], precision);
        elseif bytesPerSample == 1
            precision = '*int8';
            packets(ipacket).Waveform = int8(fread(fid, [1,header.BytesinDataPackets/bytesPerSample - 8], precision));
        elseif bytesPerSample == 4
            precision = '*int32';
            packets(ipacket).Waveform = int32(fread(fid, [1,header.BytesinDataPackets/bytesPerSample - 2], precision));
        end
            

    elseif packets(ipacket).PacketID >5120 && packets(ipacket).PacketID <5633 % The packets with these IDs indicate stimulation
        packets(ipacket).Reserved = fread(fid, [1,2], 'uchar');
        
        if bytesPerSample == 2
            precision = '*int16';
            packets(ipacket).Waveform = int16(fread(fid, [1,52], precision));
        elseif bytesPerSample == 1
            precision = '*int8';
            packets(ipacket).Waveform = int8(fread(fid, [1,52], precision));
        elseif bytesPerSample == 4
            precision = '*int32';
            packets(ipacket).Waveform = int32(fread(fid, [1,52], precision));
        end
   
        
    elseif packets(ipacket).PacketID > 4294967290 % 0xFFFFFFFF this means that the packet is a continuation of the previous packet
                                                  % This can probably occur
                                                  % if there is a very long
                                                  % stimulation waveform
        disp ('CONTINUATION OF THE PREVIOUS PACKET') 
        STOP : HAMMERTIME
    end

    are_there_more_packets = ~isempty(packets(ipacket).Timestamp);
    ipacket = ipacket +1;
end
toc

% save('packets.mat','packets','-v7.3') % This thing is huge. Improve (18GB)






% Convert the packets to events in Brainstorm format, keeping only what is necessary
tic
events = struct;
events.label   = [];

ID = [packets.PacketID];
already_checked = false(length(ID),1);

ipacket = 1;
while sum(already_checked) ~= length(ID) % Don't go through all of them if they are already taken care of
    disp(num2str(length(ID) - sum(already_checked)))
    if ~already_checked(ipacket)

        all_packets_with_similar_ID = find(ID == ID(ipacket)); % This will be used to concatenate the similar events
        
        if ID(ipacket) == 0            

            % Write the packet to events
            index = find(strcmp({events.label},['Event ' num2str(packets(ipacket).ParallelInput)])); % This is the index for the events

            if isempty(index)
                index = length(events)+1;
                events(index).label = ['Event ' num2str(packets(ipacket).ParallelInput)];
                all_packets_with_similar_ID=find(cellfun(@(x)isequal(x,packets(ipacket).ParallelInput),{packets.ParallelInput})==1);
            end
            events(index).color                  = rand(1,3);
            events(index).epochs                 = ones(1,length(all_packets_with_similar_ID));
            events(index).samples                = [packets(all_packets_with_similar_ID).Timestamp];
            events(index).times                  = [packets(all_packets_with_similar_ID).Timestamp]./header.TimeResolutionOfTimestamps;
            events(index).reactTimes             = [];
            events(index).select                 = 1;
            events(index).PacketInsertionReason  = [packets(all_packets_with_similar_ID).PacketInsertionReason];
            events(index).indicesOnPacketsStruct = all_packets_with_similar_ID;

        elseif ID(ipacket)>0  &&  ID(ipacket)<513
            % Write the packet to events
            index = find(strcmp({events.label},['Spikes Channel raw ' num2str(packets(ipacket).PacketID)]))

            if isempty(index)
                index = length(events)+1;
                events(index).label = ['Spikes Channel raw ' num2str(packets(ipacket).PacketID)];
            end
            events(index).color                    = rand(1,3);
            events(index).epochs                   = ones(1,length(all_packets_with_similar_ID));
            events(index).samples                  = [packets(all_packets_with_similar_ID).Timestamp];
            events(index).times                    = [packets(all_packets_with_similar_ID).Timestamp]./header.TimeResolutionOfTimestamps;
            events(index).reactTimes               = [];
            events(index).select                   = 1;
            events(index).UnitClassificationNumber = [packets(all_packets_with_similar_ID).UnitClassificationNumber];   
            events(index).indicesOnPacketsStruct   = all_packets_with_similar_ID;

        elseif ID(ipacket)>5120  &&  ID(ipacket)<5633
            % Write the packet to events
            index = find(strcmp({events.label},['Stimulation Channel raw ' num2str(packets(ipacket).PacketID-5120)]));

            if isempty(index)
                index = length(events)+1;
                events(index).label = ['Stimulation Channel raw' num2str(packets(ipacket).PacketID-5120)];
            end
            events(index).color                    = rand(1,3);
            events(index).epochs                   = ones(1,length(all_packets_with_similar_ID));
            events(index).samples                  = [packets(all_packets_with_similar_ID).Timestamp];
            events(index).times                    = [packets(all_packets_with_similar_ID).Timestamp]./header.TimeResolutionOfTimestamps;
            events(index).reactTimes               = [];
            events(index).select                   = 1;
            events(index).UnitClassificationNumber = [packets(all_packets_with_similar_ID).UnitClassificationNumber]; 
            events(index).indicesOnPacketsStruct   = all_packets_with_similar_ID;
            
        end

        already_checked(all_packets_with_similar_ID) = true;
    end
    ipacket = find(~already_checked,1); % This makes sure that every index for the packet struct will be for a non already accessed entry. This takes the conversion down to 9seconds
end    
toc

%
% Get rid of the empty entries and order the events alphabetically based on
% their label.
select_events_that_have_values = false(1,length(events));
for i = 1:length(events)
    select_events_that_have_values(i) = ~isempty(events(i).label);
end
events = events(select_events_that_have_values);

[~, order_alphabetically_indices] = sort(string({events.label}));
events = events(order_alphabetically_indices);

acquisition_events.header         = header;
acquisition_events.extendedheader = extendedheader;
acquisition_events.events         = events;
acquisition_events.packets        = packets;

fclose(fid);
clear info_file fid ans bytesPerSample are_there_more_packets i index ipacket order_alphabetically_indices precision select_events_that_have_values iExtendedHeader
clear header extendedheader packets ID already_checked all_packets_with_similar_ID


save ([path name 'events.mat'],'events')
% save ('events.mat','events')
end
