% speed_test
netParams.meanWeight = 1.785e-8; % from calculation for 1mV EPSP
netParams.connProbE2E = .1; % (citation)
netParams.connProbI2I = .8; % (Fino, Packer, and Yuste, 2014)
netParams.connProbGC2E = .3;
netParams.connProbGC2I = .2;

sigma = 200e-6;
f = @(D) (sqrt(2)/(sigma*sqrt(pi)))*(exp((-D.^2)./(2*sigma^2)));
netParams.connProbFunctionE2I = @(D) f(D)./max(f(D)); % (Fino, Packer, and Yuste, 2014)
netParams.gaussWeightFunction = @(D) lognrnd(log(netParams.meanWeight)-.5,1,size(D)); % (citation: Koulakov)

net = EVLIFnetwork();

net.addGroup('E',800,'excitatory',1,'std_noise',100e-12,'depressed_synapses',true)
net.addGroup('I',200,'inhibitory',1,'std_noise',100e-12,'depressed_synapses',true)

% Predefined parameters
maxWeight = 10e-7;

% Set parameters for E --> E connection
weightRange = 1e-12:1e-12:maxWeight;
px = lognpdf(weightRange,log(netParams.meanWeight)-.5,1);
lognWeightDist = weightDistribution(weightRange, px);
randomConnParamsE2E.connProb = netParams.connProbE2E;
randomConnParamsE2E.weightDistribution = lognWeightDist;

% Set parameters for E --> I connection
gaussConnParamsE2I.connProbFunction = netParams.connProbFunctionE2I;
gaussConnParamsE2I.weightFunction = netParams.gaussWeightFunction;
gaussConnParamsE2I.useWrap = true;

% Set parameters for I --> E connection
gaussConnParamsI2E.connProbFunction = netParams.connProbFunctionE2I;
gaussConnParamsI2E.weightFunction = netParams.gaussWeightFunction;
gaussConnParamsI2E.useWrap = true;

% Set parameters for I --> I connection
randomConnParamsI2I.connProb = netParams.connProbI2I;
randomConnParamsI2I.weightDistribution = lognWeightDist;

% add connections to network object
net.connect(1,1,'random',randomConnParamsE2E);
net.connect(1,2,'gaussian',gaussConnParamsE2I);
net.connect(2,1,'gaussian',gaussConnParamsI2E);
net.connect(2,2,'random',randomConnParamsI2I);

% Add GC input via a Poisson spike generator
net.addSpikeGenerator('GC',500,'excitatory',10)
% Connect the GC input to E group
randomConnParamsGC2E.connProb = netParams.connProbGC2E;
randomConnParamsGC2E.weightDistribution = lognWeightDist;
randomConnParamsGC2I.connProb = netParams.connProbGC2I;
randomConnParamsGC2I.weightDistribution = lognWeightDist;
net.connect(-1,1,'random',randomConnParamsGC2E);
net.connect(-1,2,'random',randomConnParamsGC2I);

% use GPU
useGpu = 1;

% Initialize the network variables
[V,Vreset,Cm,Gl,El,Vth,Vth0,Vth_max,tau_ref,dth,p0,GsynE,GsynI,VsynE,VsynI,tau_synE,tau_synI,...
          Iapp,std_noise,GsynMax,Isra,tau_sra,a,b,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
          ecells,icells,spikeGenProbs,cells2record,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
          tau_plus,tau_x,tau_minus,tau_y,is_plastic,C,dt] = ...
          setupNet(net,useGpu);
      

compile_runAEVLIFNetGPU_basic(net,length(cells2record));
compileSimulator(net,useGpu,length(cells2record));

% set plasticity type even though it is not being used
plasticity_type='';

% explicitly set dt
if (useGpu)
    dt = single(1e-4);
else
    dt = 1e-4;
end

datadir = ['~/phd/easySim/test/results/'];
if (~exist(datadir,'dir'))
    mkdir(datadir)
end
nT = double(10/dt);
spkfid = fopen([datadir 'speedtest_normal.bin'],'W');
tic;
runAEVLIFNetGPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
      Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
      p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
      is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
      tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
fclose(spkfid);
time = toc;
disp(num2str(time))

spkfid = fopen([datadir 'speedtest_basic.bin'],'W');
tic;
runAEVLIFNetGPU_basic_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
              Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
fclose(spkfid);
time = toc;
disp(num2str(time))