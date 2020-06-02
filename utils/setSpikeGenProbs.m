function [spikeGenProbs] = setSpikeGenProbs(net,spikeGenProbs,spikeGenIDs,newProbs)
N = net.nNeurons;
for i=1:length(spikeGenIDs)
    s = net.spikeGeneratorInfo(spikeGenIDs(i)).start_ind - N;
    e = net.spikeGeneratorInfo(spikeGenIDs(i)).end_ind - N;
    spikeGenProbs(s:e) = newProbs(i);
end
end

