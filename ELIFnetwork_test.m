% ELIFnetwork_test
clear all;
net = ELIFnetwork();
net.addGroup('group1',1000,'excitatory',1);
net.addGroup('group2',1000,'excitatory',2);
net.addGroup('group3',1000,'excitatory',3);

weightRange = 0:.001:1;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
randConnProb = .01;
weightFunction = @() rand();
gaussConnProbFunction = @(x)(sqrt(2)/(.05*sqrt(pi)))*exp(-x.^2/(2*.05^2));

gaussConnParams.connProbFunction = gaussConnProbFunction;
gaussConnParams.weightFunction = weightFunction;
gaussConnParams.useWrap = true;

randConnParams.connProb = randConnProb;
randConnParams.weightDistribution = uniformWeightDist;

net.connect(1,1,'gaussian',gaussConnParams);
net.connect(2,2,'gaussian',gaussConnParams);
net.connect(3,3,'gaussian',gaussConnParams);
net.connect(1,2,'random',randConnParams);
net.connect(1,3,'random',randConnParams);
net.connect(2,3,'random',randConnParams);

simobj.dt = 1e-5;
useGpu = false;
[V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells] = setupNet(net,simobj,useGpu);

tic;
for i=1:(1/simobj.dt)
    [V,Gref,GsynE,GsynI,spiked] = updateNet_mex(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells);
end
sim_dur = toc;
disp(['Total sim time: ' num2str(sim_dur) '. Time per timestep = ' num2str(sim_dur/(1/simobj.dt)) ' --> ' num2str((sim_dur/(1/simobj.dt))/simobj.dt) 'x real time'])

