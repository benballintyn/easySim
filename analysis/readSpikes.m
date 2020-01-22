function [data] = readSpikes(fname,cells2record)
raw=fread(fopen(fname),inf,'int32');
if (isempty(raw))
    data = [];
    return
end
count = 1;
doSkip = false;
for i=1:length(raw)
    if (doSkip)
        doSkip = false;
        continue;
    end
    if (raw(i) == -1)
        curTime = raw(i+1);
        doSkip = true;
    else
        data(count,1) = curTime;
        data(count,2) = cells2record(raw(i));
        count=count+1;
    end
end
end

