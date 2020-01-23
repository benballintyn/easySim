% EVLIF_example

% Initialize an empty EVLIFnetwork object
net = EVLIFnetwork();

% Add an excitatory group with 1000 neurons with default parameters to
% coordinate frame 1
net.addGroup('1',1000,'excitatory',1);

% Add a second excitatory group with 1000 neurons and custom parameters to
% coordinate frame 1
net.addGroup('2',1000,'excitatory',1,'mean_tau_synE',10e-3,'std_tau_synE',1e-3);

% Add an inhibitory group with 500 neurons and default parameters to
% coordinate frame 1
net.addGroup('3',500,'inhibitory',1);

% Connect group 1 to group 2 using a gaussian connection
% type. 1 is group 1's ID and does not necessarily need to match its name
% (given in the addGroup call)
weightFunction = @(x) 1.1e-10*rand(size(x)); % function handle that returns a synaptic weight (for dGsyn matrix) given a distance
gaussConnProbFunction = @(x)(sqrt(2)/(.05*sqrt(pi)))*exp(-x.^2/(2*.05^2)); % function handle that returns a connection probability given a distance
conn_1_2_params.connProbFunction = gaussConnProbFunction; % assign connection probability function to connProbFunction field
conn_1_2_params.weightFunction = weightFunction; % assign weight function to weightFunction field
conn_1_2_params.useWrap = true; % set useWrap to true to allow distances to be wrapped around the coordinate frame (avoids boundary effects)
% call the net.connect(src_id,tgt_id,connType,params) method with
% parameters
net.connect(1,2,'gaussian',conn_1_2_params);

% now connect group 2 to group 1 with the same connection type and
% parameters
conn_2_1_params = conn_1_2_params;
net.connect(2,1,'gaussian',conn_2_1_params);

% connect the inhibitory group (group 3) to both group 1 (with 10% connection probability) 
% and group 2 (with 5% connection probability) using a random connectivity 
%(with synaptic weights drawn from a uniform distribution)
weightRange = 0:1e-15:1.1e-10;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
conn_3_1_params.connProb = .1;
conn_3_1_params.weightDistribution = uniformWeightDist;
conn_3_2_params.connProb = .05;
conn_3_2_params.weightDistribution = uniformWeightDist;
net.connect(3,1,'random',conn_3_1_params);
net.connect(3,2,'random',conn_3_2_params);

% Last, add an excitatory Poisson spike generator group of 100 neurons that
% spike at 10Hz and connect it to the inhibitory group (group 3). NOTE THAT
% SPIKE GENERATOR GROUP IDS ARE GIVEN NEGATIVE VALUES TO DISTINGUISH THEM
% FROM SIMULATED NEURON GROUPS. THEREFORE THIS GROUP WILL RECEIVE AN ID OF
% -1
net.addSpikeGenerator('spikegen1',100, 'excitatory',10);
conn_spkgen_3_params = conn_3_1_params;
net.connect(-1,3,'random',conn_spkgen_3_params); 

% Simulate this network on the GPU for 10,000 timesteps and record the spikes
% from groups 1-3 in a file 'spikes.bin'. Adding the optional input
% 'recompile' will cause the easysim function to recompile the simulation
% code into CUDA code for this particular network
ntimesteps = 10000;
useGpu = true;
spikefile = 'spikes.bin';
easysim(net,ntimesteps,useGpu,'spikefile',spikefile,'recompile',true)