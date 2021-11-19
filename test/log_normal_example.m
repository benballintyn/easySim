% log_normal_example
clear all; close all;
net = AEVLIFnetwork();

net.addGroup('1',800,'excitatory',1,'std_noise',1000e-12);
net.addGroup('2',200,'inhibitory',1,'std_noise',1000e-12);

x = 1e-12:1e-12:1e-6;
g = 1.785e-8;
lnd_px = lognpdf(x,log(g)-.5,1);
lnd_weightDist = weightDistribution(x,lnd_px);

randConnParamsE2E.connProb = .1;
randConnParamsE2E.weightDistribution = lnd_weightDist;

randConnParamsEI.connProb = .5;
randConnParamsEI.weightDistribution = lnd_weightDist;

randConnParamsI2I.connProb = .1;
randConnParamsI2I.weightDistribution = lnd_weightDist;

net.connect(1,1,'random',randConnParamsE2E);
net.connect(1,2,'random',randConnParamsEI);
net.connect(2,1,'random',randConnParamsEI);
net.connect(2,2,'random',randConnParamsI2I);

ntimesteps = 100000;
useGpu = true;
spikefile = 'spikes.bin';
sim_dir = 'results/log_normal_example';
[outputs] = easysim(net,ntimesteps,useGpu,...
                                    'sim_dir',sim_dir,...
                                    'spikefile',spikefile,...
                                    'recompile',true,...
                                    'recordVars',false,...
                                    'timeStepSize',.0001);

% retrieve spike data and compute firing rates
[spikeData] = readSpikes([sim_dir '/' spikefile],outputs.cells2record);
window = .1/outputs.dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(outputs.cells2record),ntimesteps,outputs.dt,downsampleFactor,window);
