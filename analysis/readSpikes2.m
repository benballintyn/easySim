function [data] = readSpikes2(fname,cells2record)
if (size(cells2record,1) > size(cells2record,2))
    cells2record = cells2record';
end
fid = fopen(fname);
raw=fread(fid,inf,'int32');
if (isempty(raw))
    data = [];
    return
end
nRaw = length(raw);
timestepStartInds = find(raw == -1);
timestepInds = timestepStartInds + 1;
nTimeSteps = length(timestepStartInds);
nspikes = nRaw - 2*nTimeSteps;

data = zeros(nspikes,2);
curDataSize = 0;
for i=1:nTimeSteps
    if (i ~= nTimeSteps)
        curNspikes = timestepStartInds(i+1) - timestepInds(i) - 1;
        %disp([num2str(timestepInds(i)+1) ' ' num2str(timestepStartInds(i+1)-1)])
        curCellIDs = cells2record(raw((timestepInds(i)+1):(timestepStartInds(i+1)-1)))';
        %disp(['curNspikes = ' num2str(curNspikes) ' nCurCellIDs = ' num2str(length(curCellIDs))])
        oldDataSize = curDataSize;
        curDataSize = curDataSize + curNspikes;
        data((oldDataSize + 1):curDataSize,:) = [ones(curNspikes,1)*raw(timestepInds(i)) curCellIDs];
    else
        curNspikes = nRaw - timestepInds(end);
        curCellIDs = cells2record(raw((timestepInds(i)+1):end))';
        %disp(['curNspikes = ' num2str(curNspikes) ' nCurCellIDs = ' num2str(length(curCellIDs))])
        oldDataSize = curDataSize;
        curDataSize = curDataSize + curNspikes;
        data((oldDataSize + 1):curDataSize,:) = [ones(curNspikes,1)*raw(timestepInds(i)) curCellIDs];
    end
end
end

