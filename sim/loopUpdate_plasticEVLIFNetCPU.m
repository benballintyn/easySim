function [dGsyn] = loopUpdate_plasticEVLIFNetCPU(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
              VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,tau_synI,...
              Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid) %#codegen
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
    
    % decay plasticity variables
    dr1dt = arrayfun(@rdivide,-r1,tau_plus);
    dr2dt = arrayfun(@rdivide,-r2,tau_x);
    do1dt = arrayfun(@rdivide,-o1,tau_minus);
    do2dt = arrayfun(@rdivide,-o2,tau_y);
    r1 = arrayfun(@plus,r1,dr1dt*dt);
    r2 = arrayfun(@plus,r2,dr2dt*dt);
    o1 = arrayfun(@plus,o1,do1dt*dt);
    o2 = arrayfun(@plus,o2,do2dt*dt);
    
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
    
    % update synaptic weights
    % LTD first (includes simulated neuron and spikeGenerator spikes)
    if (areAnySpikes)
        ltd1 = bsxfun(@times,A2minus,C(:,allSpikes));
        ltd2 = bsxfun(@times,ltd1,o1);
        ltd3 = bsxfun(@times,C(:,allSpikes),A3minus);
        ltd4 = bsxfun(@times,ltd3,o1);
        ltd5 = bsxfun(@times,ltd4,r2(allSpikes)');
        ltd6 = ltd2 + ltd5;
        ltd7 = bsxfun(@times,ltd6,is_plastic(:,allSpikes));
        ltd8 = bsxfun(@times,ltd7,dGsyn(:,allSpikes));
        dGsyn(:,allSpikes) = dGsyn(:,allSpikes) - ltd8;
    end
    % LTP next (includes simulated neuron spikes only)
    if (areSimSpikes)
        ltp1 = bsxfun(@times,C(spiked,:),A2plus');
        ltp2 = bsxfun(@times,ltp1,r1');
        ltp3 = bsxfun(@times,C(spiked,:),A3plus');
        ltp4 = bsxfun(@times,ltp3,r1');
        ltp5 = bsxfun(@times,ltp4,o2(spiked));
        ltp6 = ltp2 + ltp5;
        ltp7 = bsxfun(@times,ltp6,is_plastic(spiked,:));
        ltp8 = bsxfun(@times,ltp7,dGsyn(spiked,:));
        dGsyn(spiked,:) = dGsyn(spiked,:) + ltp8;
    end
    
    if (areSimSpikes)
        V(spiked) = Vreset(spiked);
        Vth(spiked) = Vth_max(spiked);
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
    
    % Compute total synaptic current for each neuron
    Isyn = arrayfun(@plus,arrayfun(@times,GsynE,arrayfun(@minus,VsynE,V)),arrayfun(@times,GsynI,arrayfun(@minus,VsynI,V)));
    
    % add noise to any input currents
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
    dVdt = arrayfun(@times,f1,f10);
    
    V = arrayfun(@plus,V,dVdt*dt);
    V = arrayfun(@max,Vreset,V); % bound membrane potentials to be >= than the reset value
    
    % update synapse strengths and plasticity variables
    if (strcmp(plasticity_type,'all-to-all'))
        r1(spiked) = r1(spiked) + 1;
        r2(spiked) = r2(spiked) + 1;
        o1(spiked) = o1(spiked) + 1;
        o2(spiked) = o2(spiked) + 1;
    elseif (strcmp(plasticity_type,'nearest'))
        r1(spiked) = 1;
        r2(spiked) = 1;
        o1(spiked) = 1;
        o2(spiked) = 1;
    end
    
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