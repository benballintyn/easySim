function [] = makeFakeSpikeData(nTimeSteps,cells,maxFrac)
ncells = length(cells);

fid = fopen('results/test_spikes.bin','W');

for i=1:nTimeSteps
    fwrite(fid,-1,'int32');
    fwrite(fid,i,'int32');
    nSpikes = ceil(rand*maxFrac*ncells);
    spikeIDs = randperm(ncells)';
    spikeIDs = spikeIDs(1:nSpikes);
    fwrite(fid,spikeIDs,'int32');
end
end

