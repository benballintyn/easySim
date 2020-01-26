% noise_example
clear all;

net = EVLIFnetwork();

net.addGroup('1',10,'excitatory',1,'std_noise',50e-12);

ntimesteps = 1000000;
useGpu = false;
spikefile = 'spikes.bin';
sim_dir = 'results/noise_example';
[dt,cells2record,sim_dir,recordV,recordVth,iappRecord] = easysim(net,ntimesteps,useGpu,...
                                    'sim_dir',sim_dir,...
                                    'spikefile',spikefile,...
                                    'recompile',true,...
                                    'recordVars',true);

% retrieve spike data and compute firing rates
[spikeData] = readSpikes([sim_dir '/' spikefile],cells2record);
window = .1/dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(cells2record),ntimesteps,dt,downsampleFactor,window);
