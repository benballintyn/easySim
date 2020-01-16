% ELIFnetwork_test
clear all;
net = ELIFnetwork();
net.addGroup('group1',2000,'excitatory',1);
net.addGroup('group2',2000,'excitatory',2);
net.addGroup('group3',2000,'excitatory',3);
net.addGroup('group4',2000,'inhibitory',4);

weightRange = 0:1e-15:1.1e-10;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
randConnProb = .01;
weightFunction = @(x) 1.1e-10*rand(size(x));
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
net.connect(1,4,'random',randConnParams);
net.connect(4,1,'random',randConnParams);
net.connect(4,2,'random',randConnParams);
net.connect(4,3,'random',randConnParams);

useGpu = 0;
disp(['useGPU = ' num2str(useGpu)])
dt=single(1e-4);
simTime = 1;
if (useGpu)
    %allSpikes = gpuArray(zeros(net.nNeurons,nT,'single'));
    simobj.dt = dt;
else
    %allSpikes = zeros(net.nNeurons,nT);
    simobj.dt = dt;
end
nT = ceil((simTime/simobj.dt));
if (useGpu)
    allSpikes = gpuArray(zeros(net.nNeurons,nT,'single'));
    allVs = gpuArray(zeros(net.nNeurons,nT,'single'));
else
    allSpikes = zeros(net.nNeurons,nT);
    allVs = zeros(net.nNeurons,nT);
end

[V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells] = setupNet(net,simobj,useGpu);

disp('starting simulation')
tic;
%{
for i=1:nT
    [V,Gref,GsynE,GsynI,spiked] = updateNet_mex(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells);
    allVs(:,i) = V;
    allSpikes(:,i) = spiked;
    if (mod(i,1000)==0)
        disp(num2str(i))
    end
end
%}
if (~useGpu)
    compile_easySim_noGpuArray(net.nNeurons);
    [allVs,allSpikes] = loopUpdateNet_mex(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                      dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells,nT);
else
    [allVs,allSpikes] = loopUpdateNet(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells,nT);
end
sim_dur = toc;
disp(['Total sim time: ' num2str(sim_dur) '. Time per timestep = ' num2str(sim_dur/(simTime/simobj.dt)) ' --> ' num2str((sim_dur/(simTime/simobj.dt))/simobj.dt) 'x real time'])
figure;
imagesc(allSpikes)
