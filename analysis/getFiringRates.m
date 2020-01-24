function [frs] = getFiringRates(spikeData,cells2record,nTimesteps,dt,downsampleFactor,window)
% This function takes in the spike data form readSpikes as well as
% additional parameters to convert the spike times into continuous firing
% rate traces. It does this by using the alpha function method from Dayan &
% Abbott Chapter 1. This is a causal filter that only uses past spikes to
% compute the firing rate.
% getFiringRates(spikeData,nNeurons,nTimesteps,dt,downsampleFactor,window)
%   INPUTS:
%       spikeData        - matrix that is output from readSpikes
%
%       cells2record     - vector containing the id #'s of all the neurons that
%                          were recorded
%
%       nTimesteps       - # of timesteps in the simulation. Note that
%                          max(spikeData(:,1)) should be <= nTimesteps
%
%       dt               - timestep size used in the simulation
%
%       downsampleFactor - factor by which to downsample time. e.g. if
%                          nTimesteps is 10,000 and downsampleFactor is 10,
%                          the output frs will have 1,000 time points
%
%       window           - factor controlling how much spikes are smoothed
%
%   OUTPUTS:
%       frs - (nTimesteps/downsampleFactor) x length(cells2record) matrix
%             of firing rates. Rows give timepoints and columns give
%             individual neurons
nNeurons = length(cells2record);
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

