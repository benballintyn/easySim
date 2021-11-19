% Parameters that work to have noise-induced switches (no input)
% E --> E
%   weightRange = 5e-9:1e-12:5e-8
%   connProb = .5, meanWt = mean(weightRange), stdWt = .3*meanWt
% E --> I
%   weightRange = 5e-9:1e-12:5e-8
%   connProb = .51, meanWt = mean(weightRange), stdWt = .3*meanWt
% I --> E
%   weightRange = 8e-8:1e-12:8e-7;
%   connProb = .5, meanWt = mean(weightRange), stdWt = .3*meanWt
clear all;

net = AEVLIFnetwork();

% Add excitatory stay and switch groups and inhibitory stay and switch
% groups. Excitatory 
net.addGroup('stayE',400,'excitatory',1,'std_noise',1000e-12,'depressed_synapses',true);
net.addGroup('switchE',400,'excitatory',1,'std_noise',1000e-12,'depressed_synapses',true);
net.addGroup('stayI',100,'inhibitory',1);
net.addGroup('switchI',100,'inhibitory',1);

% Add excitatory Poisson inputs to both stayE and switchE
net.addSpikeGenerator('stayInput',100,'excitatory',10);
net.addSpikeGenerator('switchInput',100,'excitatory',10);

% Define weight distribution and connection parameters for random
% self-connections
weightRange = 5e-9:1e-12:5e-8;
meanWt = mean(weightRange); stdWt = .3*meanWt;
px = (1/(stdWt.*sqrt(2*pi)))*exp(-.5.*((weightRange - meanWt)./stdWt).^2);
gaussWeightDist = weightDistribution(weightRange,px);
randomConnParams.connProb = .5;
randomConnParams.weightDistribution = gaussWeightDist;

% Self-connections for stayE and switchE
net.connect(1,1,'random',randomConnParams);
net.connect(2,2,'random',randomConnParams);

% Define connection parameters for random E-I connections
weightRange = 5e-9:1e-12:5e-8;
meanWt = mean(weightRange); stdWt = .3*meanWt;
px = (1/(stdWt.*sqrt(2*pi)))*exp(-.5.*((weightRange - meanWt)./stdWt).^2);
gaussWeightDist = weightDistribution(weightRange,px);
randomConnParams.connProb = .508;
randomConnParams.weightDistribution = gaussWeightDist;

% Make E-I
net.connect(1,4,'random',randomConnParams);
net.connect(2,3,'random',randomConnParams);

% Define connection parameters for random I-E connections
weightRange = 8e-8:1e-12:8e-7;
meanWt = mean(weightRange); stdWt = .3*meanWt;
px = (1/(stdWt.*sqrt(2*pi)))*exp(-.5.*((weightRange - meanWt)./stdWt).^2);
gaussWeightDist = weightDistribution(weightRange,px);
randomConnParams.connProb = .5;
randomConnParams.weightDistribution = gaussWeightDist;

% Make I-E connections
net.connect(4,2,'random',randomConnParams);
net.connect(3,1,'random',randomConnParams);

% Connect spikeGenerators to the appropriate pools with own parameters
weightRange = 1e-8:1e-11:1e-7;
meanWt = mean(weightRange); stdWt = .3*meanWt;
px = (1/(stdWt.*sqrt(2*pi)))*exp(-.5.*((weightRange - meanWt)./stdWt).^2);
gaussWeightDist = weightDistribution(weightRange,px);
randomConnParams.connProb = 0;
randomConnParams.weightDistribution = gaussWeightDist;
net.connect(-1,1,'random',randomConnParams);
%net.connect(-2,2,'random',randomConnParams);

% Set useGPU flag
useGPU = true;

% Setup network variables
% Initialize the network variables
[V,Vreset,Cm,Gl,El,Vth,Vth0,Vth_max,tau_ref,dth,p0,GsynE,GsynI,VsynE,VsynI,tau_synE,tau_synI,...
          Iapp,std_noise,GsynMax,Isra,tau_sra,a,b,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
          ecells,icells,spikeGenProbs,cells2record,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
          tau_plus,tau_x,tau_minus,tau_y,is_plastic,C,dt] = ...
          setupNet(net,useGPU);

% Compile simulation code
compileSimulator(net,useGPU,length(cells2record));

% set plasticity type even though it is not being used
plasticity_type='';

% explicitly set dt
if (useGPU)
    dt = single(1e-4);
else
    dt = 1e-4;
end

% Make directory to store data in
datadir = ['~/phd/easySim/results/stay_switch_test/'];
if (~exist(datadir,'dir'))
    mkdir(datadir)
end

% Set spikeGenerator spike probabilities for ON and OFF periods
spikeGenProbsOFF = setSpikeGenProbs(net,spikeGenProbs,[1 2],[0 0]);
spikeGenProbsON = setSpikeGenProbs(net,spikeGenProbsOFF,[1 2],[20*dt 20*dt]);

% Reinitialize the network variables
[V,Vth,Isra,GsynE,GsynI,D,F,r1,r2,o1,o2,Iapp] = resetVars(net,useGPU);

% Set number of timesteps to simulate
nT_stim = double(2/dt);

% Open file in which spikes will be stored for the ON period
spikeONfname = [datadir 'stay_switch_test_spikesON.bin'];
spkfid = fopen(spikeONfname,'W');

tic;
% Run the simulation for 2 seconds with the spikeGenerators on
[~,V,Vth,Isra,GsynE,GsynI,D,F,r1,r2,o1,o2] = runAEVLIFNetGPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
              Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbsON,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT_stim,spkfid);
timeON = toc;
fclose(spkfid); % close the spike file


% Run the simulation for 98 seconds with the spikeGenerators off
nT_switches = double(500/dt);
spikeOFFfname = [datadir 'stay_switch_test_spikesOFF.bin'];
spkfid = fopen(spikeOFFfname,'W'); % create a new spike file for the off period

tic;
[~,V,Vth,Isra,GsynE,GsynI,D,F,r1,r2,o1,o2] = runAEVLIFNetGPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
              Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbsOFF,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT_switches,spkfid);
timeOFF = toc;
fclose(spkfid); % close the spike filed

totalTime = timeON + timeOFF;
simTime = (nT_stim + nT_switches)*dt;
disp(['Simulated time = ' num2str(simTime) ' actual time = ' num2str(totalTime) ' : ' num2str(totalTime/simTime) 'x real time'])


% make sure all files are closed
fclose('all');

% Save the network and connectivity matrix
save([datadir 'net.mat'],'net','-mat')
save([datadir 'GsynMax.mat'],'GsynMax','-mat')
nT = nT_stim + nT_switches;
save([datadir 'nT.mat'],'nT','-mat')
save([datadir 'dt.mat'],'dt','-mat')

% Read the spikes from both spike files and concatenate them. NOTE, if there are many spikes this
% will produce a very tall matrix (Nspikes x 2)
cells2record = gather(cells2record);
save([datadir 'cells2record.mat'],'cells2record','-mat')
[dataON] = readSpikes2(spikeONfname,cells2record);
[dataOFF] = readSpikes2(spikeOFFfname,cells2record);
dataOFF(:,1) = dataOFF(:,1) + nT_stim; % offset spike timesteps for OFF period by nT_stim
allSpikes = [dataON; dataOFF];

% Turn spike data into firing rates
downsampleFactor = 20;
window = 100;
save([datadir 'downsampleFactor.mat'],'downsampleFactor','-mat')
save([datadir 'window.mat'],'window','-mat')
[frs] = getFiringRates(allSpikes,length(cells2record),nT_stim + nT_switches,dt,downsampleFactor,window);
tvec = dt:dt*downsampleFactor:(nT_stim + nT_switches)*dt;

figure;
for i=1:4
    plot(tvec,mean(frs(:,net.groupInfo(i).start_ind:net.groupInfo(i).end_ind),2))
    hold on;
end