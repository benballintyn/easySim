% v1_test
clear all;

net = AEVLIFnetwork();

net.addGroup('1',800,'excitatory',1,'std_noise',100e-12,'depressed_synapses',true);

net.addGroup('2',200,'inhibitory',1,'std_noise',100e-12);

weightRange = 1e-9:1e-11:1.8e-7;
%uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
meanWt = .5e-7; stdWt = .5*meanWt;
px = (1/(stdWt.*sqrt(2*pi)))*exp(-.5.*((weightRange - meanWt)./stdWt).^2);
gaussWeightDist = weightDistribution(weightRange,px);
clusterConnParams.nClusters = 8;
clusterConnParams.intraConnProb = .3;
clusterConnParams.interConnProb = .01;
clusterConnParams.intraWeightDist = gaussWeightDist;
clusterConnParams.interWeightDist = gaussWeightDist;
clusterConnParams.is_plastic = true;

net.connect(1,1,'clustered',clusterConnParams);

randomConnParams.connProb = .1;
randomConnParams.weightDistribution = gaussWeightDist;
net.connect(1,2,'random',randomConnParams);
net.connect(2,1,'random',randomConnParams);

net.addSpikeGenerator('input',500,'excitatory',10)

net.connect(-1,1,'random',randomConnParams);

ntimesteps = 100000;
useGpu = true;
spikefile = 'spikes.bin';
sim_dir = 'results/v1_test';
[outputs] = easysim(net,ntimesteps,useGpu,'recompile',true,...
                                          'sim_dir',sim_dir,...
                                          'spikefile',spikefile);

[spikeData] = readSpikes([sim_dir '/' spikefile],outputs.cells2record);
window = .1/outputs.dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(outputs.cells2record),ntimesteps,outputs.dt,downsampleFactor,window);


      