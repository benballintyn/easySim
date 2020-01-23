function [] = loopUpdateEVLIFNetGPU_fast(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid) %#codegen

coder.gpu.kernelfun; % for code generation

N = size(V,1); % # of simulated neurons
nSpikeGen = length(spikeGenProbs); % # of poisson spike generator neurons
useSpikeGen = (nSpikeGen > 0); % determine if there are any spike generators
n2record = length(cells2record); % # of neurons to record
useRecord = (n2record > 0); % determine if any neurons should be recorded

% if no spike file was given, don't record any spikes
if (spkfid < 0)
    useRecord = false;
end

% Loop through nT timepoints
for i=2:(nT+1)
    
    % Update spike thresholds
    vth1 = arrayfun(@minus,Vth0,Vth);
    dVthdt = arrayfun(@rdivide,vth1,tau_ref);
    Vth = arrayfun(@plus,Vth,dVthdt*dt);
    
    spiked = (V > Vth); % determine simulated neurons that spiked
    
    % check if any spike generators spiked
    if (useSpikeGen)
        spikeGenSpikes = (rand(nSpikeGen,1) < spikeGenProbs);
        allSpikes = [spiked; spikeGenSpikes];
    else
        allSpikes = spiked;
    end
    
    areSpikes = any(spiked); % determine if there were any spikes
    if (areSpikes)
        V(spiked) = Vreset(spiked); % reset membrane voltages
        Vth(spiked) = Vth_max(spiked); % set spike threshold to max
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
    
    Isyn = arrayfun(@plus,arrayfun(@times,GsynE,arrayfun(@minus,VsynE,V)),arrayfun(@times,GsynI,arrayfun(@minus,VsynI,V)));
    Iapp = arrayfun(@plus,Iapp,arrayfun(@times,std_noise,randn(N,1)));
    %Iapp = arrayfun(@plus,Iapp,normrnd(0,std_noise));
    
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
    
    if (areSpikes)
        if (useRecord)
            fwrite(spkfid,-1,'int32');
            fwrite(spkfid,i,'int32');
            fwrite(spkfid,find(spiked(cells2record)),'int32');
        end
    end
    
end
end