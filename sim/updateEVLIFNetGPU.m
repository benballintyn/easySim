function [V,Vth,GsynE,GsynI,spiked] = updateEVLIFNetGPU(V,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,dt,ecells,icells,spikeGenProbs) %#codegen

coder.gpu.kernelfun;

nSpikeGen = length(spikeGenProbs);
useSpikeGen = (nSpikeGen > 0);

% update spike thresholds
vth1 = arrayfun(@minus,Vth0,Vth);
dVthdt = arrayfun(@rdivide,vth1,tau_ref);
Vth = arrayfun(@plus,Vth,dVthdt*dt);
    
spiked = (V > Vth);
if (useSpikeGen)
    spikeGenSpikes = (rand(nSpikeGen,1) < spikeGenProbs);
    allSpikes = [spiked; spikeGenSpikes];
else
    allSpikes = spiked;
end

areSpikes = any(spiked);
if (areSpikes)
    V(spiked) = -.08;
    Vth(spiked) = Vth_max(spiked);
end

e_spiked = logical(allSpikes.*ecells);
i_spiked = logical(allSpikes.*icells);
dGsynEdt = arrayfun(@rdivide,-GsynE,tau_synE);
dGsynIdt = arrayfun(@rdivide,-GsynI,tau_synI);
GsynE = arrayfun(@plus,GsynE,dGsynEdt*dt);
GsynI = arrayfun(@plus,GsynI,dGsynIdt*dt);

if (areSpikes)
    dGsynE_sum = sum(dGsyn(:,e_spiked),2);
    dGsynI_sum = sum(dGsyn(:,i_spiked),2);
    GsynE = arrayfun(@plus,GsynE,dGsynE_sum);
    GsynE = arrayfun(@min,GsynE,maxGsynE);
    GsynI = arrayfun(@plus,GsynI,dGsynI_sum);
    GsynI = arrayfun(@min,GsynI,maxGsynI);
end
Isyn = arrayfun(@times,GsynE,arrayfun(@minus,VsynE,V)) + arrayfun(@times,GsynI,arrayfun(@minus,VsynI,V));

f1 = (1./Cm);
f2 = arrayfun(@minus,El,V);
f3 = arrayfun(@minus,V,Vth);
f4 = arrayfun(@rdivide,f3,dth);
f5 = arrayfun(@exp,f4);
f6 = arrayfun(@times,dth,f5);
f7 = arrayfun(@plus,f2,f6);
f8 = arrayfun(@times,Gl,f7);
f9 = arrayfun(@plus,f8,Isyn);
f10 = arrayfun(@plus,f9,Iapp);
dVdt = arrayfun(@times,f1,f10);

V = arrayfun(@plus,V,dVdt*dt);
V = arrayfun(@max,Vreset,V); % bound membrane potentials to be >= than the reset value
end