function [spikeGenProbs] = getSpikeGenProbs(net,dt)
N = net.nNeurons;
for i=1:net.nSpikeGenerators
    s = net.spikeGeneratorInfo(i).start_ind - N;
    e = net.spikeGeneratorInfo(i).end_ind - N;
    spikeGenProbs(s:e) = net.spikeGeneratorInfo(i).firingRate*dt;
end
end

