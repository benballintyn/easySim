% ELIFnetwork_test
clear all;
net = EVLIFnetwork();
net.addGroup('group1',1000,'excitatory',1);
net.addGroup('group2',1000,'excitatory',2);
net.addGroup('group3',1000,'excitatory',3);
net.addGroup('group4',1000,'inhibitory',4);

weightRange = 0:1e-12:10e-8;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
randConnProb = .1;
weightFunction = @(x) 1e-8*rand(size(x));
gaussConnProbFunction = @(x)(sqrt(2)/(.05*sqrt(pi)))*exp(-x.^2/(2*.05^2));

gaussConnParams.connProbFunction = gaussConnProbFunction;
gaussConnParams.weightFunction = weightFunction;
gaussConnParams.useWrap = true;

randConnParams.connProb = randConnProb;
randConnParams.weightDistribution = uniformWeightDist;

randConnParams2.connProb = .1;
randConnParams2.weightDistribution = uniformWeightDist;

net.connect(1,1,'gaussian',gaussConnParams);
net.connect(2,2,'gaussian',gaussConnParams);
net.connect(3,3,'gaussian',gaussConnParams);
net.connect(1,2,'random',randConnParams);
net.connect(1,3,'random',randConnParams);
net.connect(2,3,'random',randConnParams);
net.connect(1,4,'random',randConnParams);
net.connect(4,1,'random',randConnParams2);
net.connect(4,2,'random',randConnParams2);
net.connect(4,3,'random',randConnParams2);

weightRange = 0:1e-12:4.5e-8;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
spkInputParams.connProb = .05;
spkInputParams.weightDistribution = uniformWeightDist;
net.addSpikeGenerator('spk1',1000,'excitatory',10);
net.connect(-1,1,'random',spkInputParams);

useGpu = 0;
disp(['useGPU = ' num2str(useGpu)])

[V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
          setupEVLIFNet(net,useGpu);
simTime = 1;
nT = ceil((simTime/dt));

spkfid = fopen('test.bin','W');
if (~useGpu)
    disp('compiling')
    compile_loopUpdateEVLIFNetCPU();
    disp('starting simulation')
    tic;
    loopUpdateEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
    sim_dur = toc;
    fclose(spkfid);
    disp(['loopUpdateEVLIFNetCPU: Total sim time: ' num2str(sim_dur) '. Time per timestep = ' num2str(sim_dur/(simTime/dt)) ' --> ' num2str((sim_dur/(simTime/dt))/dt) 'x real time'])
else
    compile_loopUpdateEVLIFNetGPU_fast(net.nNeurons,length(spikeGenProbs),length(cells2record));
    tic;
    loopUpdateEVLIFNetGPU_fast(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
    sim_dur = toc;
    fclose(spkfid);
    disp(['loopUpdateEVLIFNetGPU_fast: Total sim time: ' num2str(sim_dur) '. Time per timestep = ' num2str(sim_dur/(simTime/simobj.dt)) ' --> ' num2str((sim_dur/(simTime/simobj.dt))/simobj.dt) 'x real time'])
end

% retrieve spike data and compute firing rates
[spikeData] = readSpikes('test.bin',cells2record);
window = .1/dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(cells2record),nT,dt,downsampleFactor,window);