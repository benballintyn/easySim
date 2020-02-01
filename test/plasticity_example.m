% plasticity_example
clear all;

net = AEVLIFnetwork();

net.addGroup('1',1000,'excitatory',1);
net.addGroup('2',1000,'excitatory',1);
net.addGroup('3',500,'inhibitory',1);

weightRange = 1e-9:1e-11:1e-8;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
plastic_conn_params.connProb = .1;
plastic_conn_params.weightDistribution = uniformWeightDist;
plastic_conn_params.is_plastic = true; % key line
net.connect(1,2,'random',plastic_conn_params);

static_conn_params = plastic_conn_params;
static_conn_params.is_plastic = false;
net.connect(3,1,'random',static_conn_params);
net.connect(3,2,'random',static_conn_params);

net.print(true)

useGpu = true;
nT = 100000;
spikefile = 'spikes.bin';
sim_dir = 'results/plasticity_example';
[dt,cells2record] = easysim(net,nT,useGpu,'sim_dir',sim_dir,...
                                          'spikefile',spikefile,...
                                          'recompile',true);

% retrieve spike data and compute firing rates
[spikeData] = readSpikes([sim_dir '/' spikefile],cells2record);
window = .1/dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(cells2record),ntimesteps,dt,downsampleFactor,window);
