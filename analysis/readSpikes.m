function [data] = readSpikes(fname,cells2record)
% This function is used to read the spike time data from a .bin file
% readSpikes(fname,cells2record)
%   INPUTS:
%       fname        - filename containing the spikes
%       cells2record - vector containing the id #'s of all the neurons that
%                      were recorded
%
%   OUTPUTS:
%       data - ? x 2 matrix where in each row the 1st value is the timestep
%              index (NOT TIME IN SECONDS) of a spike and the 2nd value is
%              the id # of the neuron that spiked. The size of the 1st
%              dimension of this matrix will be proportional to the # of
%              spikes in the simulation
fid = fopen(fname);
raw=fread(fid,inf,'int32');
if (isempty(raw))
    data = [];
    return
end
nRaw = length(raw);
modVal = floor(nRaw/100);
count = 1;
doSkip = false;
for i=1:nRaw
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
    if (mod(i,modVal) == 0)
        disp([num2str((i/nRaw)*100) '% done'])
    end
end
fclose('all');
end

