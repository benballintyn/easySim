function [frs] = getFiringRates(spikeData,nNeurons,nTimesteps,dt,downsampleFactor,window)

frs = zeros(nNeurons,nTimesteps);
if (isempty(spikeData))
    frs = downsample(frs',downsampleFactor);
    return
end
alpha=1/window;
tau=-1000:1:1000;
w=(alpha^2)*tau.*exp(-alpha*tau); %alpha function
w(w<0) = 0;

for i=1:nNeurons
    curSpikeTimes = spikeData((spikeData(:,2) == i),1);
    spikes = zeros(1,nTimesteps);
    spikes(curSpikeTimes) = 1;
    frs(i,:) = conv(spikes,w,'same');
end
frs = downsample(frs',downsampleFactor);
frs = frs./dt;

end

