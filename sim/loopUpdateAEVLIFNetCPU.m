function [] = loopUpdateAEVLIFNetCPU(V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid) %#codegen

N = size(V,1);
nSpikeGen = length(spikeGenProbs);
useSpikeGen = (nSpikeGen > 0);
n2record = length(cells2record);
useRecord = (n2record > 0);
nSimulatedSpikes = 0;
nGeneratedSpikes = 0;

% if no spike file was given, don't record any spikes
if (spkfid < 0)
    useRecord = false;
end

% Loop through nT timepoints
for i=1:nT
    
    % Update spike thresholds
    vth1 = arrayfun(@minus,Vth0,Vth);
    dVthdt = arrayfun(@rdivide,vth1,tau_ref);
    Vth = arrayfun(@plus,Vth,dVthdt*dt);
    
    % Update adaptation currents
    sra1 = arrayfun(@minus,V,El);
    sra2 = arrayfun(@times,a,sra1);
    sra3 = arrayfun(@minus,sra2,Isra);
    dIsradt = arrayfun(@rdivide,sra3,tau_sra);
    Isra = arrayfun(@plus,Isra,dIsradt*dt);
    
    spiked = (V > Vth);
    
    % check if any spike generators spiked
    if (useSpikeGen)
        spikeGenSpikes = (rand(nSpikeGen,1) < spikeGenProbs);
        allSpikes = [spiked; spikeGenSpikes];
        nGeneratedSpikes = nGeneratedSpikes + sum(spikeGenSpikes);
    else
        allSpikes = spiked;
    end
    nSimulatedSpikes = nSimulatedSpikes + sum(spiked);
    
    areSimSpikes = any(spiked);
    areAnySpikes = any(allSpikes);
    if (areSimSpikes)
        V(spiked) = Vreset(spiked);
        Vth(spiked) = Vth_max(spiked);
        Isra(spiked) = arrayfun(@plus,Isra(spiked),b(spiked));
    end
    
    e_spiked = logical(allSpikes.*ecells);
    i_spiked = logical(allSpikes.*icells);
    dGsynEdt = arrayfun(@rdivide,-GsynE,tau_synE);
    dGsynIdt = arrayfun(@rdivide,-GsynI,tau_synI);
    
    GsynE = arrayfun(@plus,GsynE,dGsynEdt*dt);
    GsynI = arrayfun(@plus,GsynI,dGsynIdt*dt);

    if (areAnySpikes)
        dGsynE_sum = sum(dGsyn(:,e_spiked),2);
        dGsynI_sum = sum(dGsyn(:,i_spiked),2);
        GsynE = arrayfun(@plus,GsynE,dGsynE_sum);
        GsynE = arrayfun(@min,GsynE,maxGsynE);
        GsynI = arrayfun(@plus,GsynI,dGsynI_sum);
        GsynI = arrayfun(@min,GsynI,maxGsynI);
    end
    
    Isyn = arrayfun(@plus,arrayfun(@times,GsynE,arrayfun(@minus,VsynE,V)),arrayfun(@times,GsynI,arrayfun(@minus,VsynI,V)));
    
    curIapp = Iapp;
    curIapp = arrayfun(@plus,curIapp,arrayfun(@times,std_noise,randn(N,1)));
    
    f1 = (1./Cm);
    f2 = arrayfun(@minus,El,V);
    f3 = arrayfun(@minus,V,Vth);
    f4 = arrayfun(@rdivide,f3,dth);
    f5 = arrayfun(@exp,f4);
    f6 = arrayfun(@times,dth,f5);
    f7 = arrayfun(@plus,f2,f6);
    f8 = arrayfun(@times,Gl,f7);
    f9 = arrayfun(@plus,f8,Isyn);
    f10 = arrayfun(@plus,f9,curIapp);
    f11 = arrayfun(@minus,f10,Isra);
    dVdt = arrayfun(@times,f1,f11);
    
    V = arrayfun(@plus,V,dVdt*dt);
    V = arrayfun(@max,Vreset,V); % bound membrane potentials to be >= than the reset value
    if (areSimSpikes)
        if (useRecord)
            fwrite(spkfid,-1,'int32');
            fwrite(spkfid,i,'int32');
            fwrite(spkfid,find(spiked(cells2record)),'int32');
        end
    end
    
end
fprintf('%i simulated spikes\n',int32(nSimulatedSpikes));
fprintf('%i generated spikes\n',int32(nGeneratedSpikes));
end