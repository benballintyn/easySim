function [spikeGenProbs] = setSpikeGenProbs(net,spikeGenProbs,spikeGenIDs,newProbs)
N = net.nNeurons;
for i=1:spikeGenIDs
    s = net.spikeGeneratorInfo(i).start_ind - N;
    e = net.spikeGeneratorInfo(i).end_ind - N;
    spikeGenProbs(s:e) = newProbs(i);
end
end

