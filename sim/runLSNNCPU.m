function [W,recordV,recordVth,recordD,recordF,recordIsyn,recordIapp] = runLSNNCPU(V_temp,Vreset,Vth_temp,Vth0,...
              t_refrac,tau_m,beta,tau_a,W_temp,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,nT,spkfid,recordVars) %#codegen

W = W_temp;
V = V_temp;
Vth = Vth_temp;

N = size(V,1); % # of simulated neurons
timeSinceSpike = nan(N,1);
rho = exp(-dt./tau_a);
alpha = exp(-dt./tau_m);

nSpikeGen = length(spikeGenProbs); % # of poisson spike generator neurons
useSpikeGen = (nSpikeGen > 0); % determine if there are any spike generators
n2record = length(cells2record); % # of neurons to record
useRecord = (n2record > 0); % determine if any neurons should be recorded
nSimulatedSpikes = 0;
nGeneratedSpikes = 0;
if (~isnan(C))
    usePlasticity = true;
else
    usePlasticity = false;
end
if (recordVars)
    recordV = zeros(N,nT);
    recordVth = zeros(N,nT);
    recordIsyn = zeros(N,nT);
    recordIapp = zeros(N,nT);
    if (usePlasticity) % NEED TO CHANGE THIS FOR E-PROP
        recordr1 = zeros(N,nT);
        recordr2 = zeros(N,nT);
        recordo1 = zeros(N,nT);
        recordo2 = zeros(N,nT);
    else
        recordr1 = nan;
        recordr2 = nan;
        recordo1 = nan;
        recordo2 = nan;
    end
else
    recordV = zeros(N,1);
    recordVth = zeros(N,1);
    recordIsyn = zeros(N,1);
    recordIapp = zeros(N,1);
    if (usePlasticity) % NEED TO CHANGE THIS FOR E-PROP
        recordr1 = zeros(N,1);
        recordr2 = zeros(N,1);
        recordo1 = zeros(N,1);
        recordo2 = zeros(N,1);
    else
        recordr1 = nan;
        recordr2 = nan;
        recordo1 = nan;
        recordo2 = nan;
    end
end
% if no spike file was given, don't record any spikes
if (spkfid < 0)
    useRecord = false;
end

% Loop through nT timepoints
for i=1:nT
    
    % Update spike thresholds
    Vth = Vth0 + beta.*a; % Threshold update (Bellec et al., 2020)
    
    if (usePlasticity) % Need to change for LSNNs
        % decay plasticity variables
        r1 = r1 + dt*(-r1./tau_plus);
        r2 = r2 + dt*(-r2./tau_x);
        o1 = o1 + dt*(-o1./tau_minus);
        o2 = o2 + dt*(-o2./tau_y);
    end
    
    spiked = (V > Vth); % determine simulated neurons that spiked
    
    % check if any spike generators spiked
    if (useSpikeGen)
        spikeGenSpikes = (rand(nSpikeGen,1) < spikeGenProbs);
        allSpikes = [spiked; spikeGenSpikes];
        nGeneratedSpikes = nGeneratedSpikes + sum(spikeGenSpikes);
    else
        allSpikes = spiked;
    end
    nSimulatedSpikes = nSimulatedSpikes + sum(spiked);
    
    areSimSpikes = any(spiked); % determine if there were any spikes
    areAnySpikes = any(allSpikes);
    
    if (areSimSpikes) % If there are any simulated spikes
        V(spiked) = Vreset(spiked); % reset membrane voltages of spiking neurons
    end
    
    % Compute total synaptic current for each neuron
    Isyn = W*allSpikes;
    
    % add noise to any input currents
    curIapp = Iapp;
    curIapp = curIapp + std_noise.*randn(N,1);
    
    % update membrane voltages
    V = alpha.*V + Isyn + curIapp - spiked.*Vth;
    V = max(V,Vreset);

    % update adaptation variables
    a = rho.*a + spiked;
    
    if (usePlasticity)
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
            ltd8 = bsxfun(@times,ltd7,GsynMax(:,allSpikes));
            GsynMax(:,allSpikes) = GsynMax(:,allSpikes) - ltd8;
            %ltd = (A2minus.*C(:,allSpikes).*o1  + C(:,allSpikes).*A3minus.*o1.*r2(allSpikes)').*is_plastic(:,allSpikes).*GsynMax(:,allSpikes);
            %GsynMax(:,allSpikes) = GsynMax(:,allSpikes) - ltd;
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
            ltp8 = bsxfun(@times,ltp7,GsynMax(spiked,:));
            GsynMax(spiked,:) = GsynMax(spiked,:) + ltp8;
            %ltp = (C(spiked,:).*A2plus'.*r1' + C(spiked,:).*A3plus'.*r1'.*o2(spiked)).*is_plastic(spiked,:).*GsynMax(spiked,:);
            %GsynMax(spiked,:) = GsynMax(spiked,:) + ltp;
        end

        % update plasticity variables
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
    end
    
    % record variables if requested
    if (recordVars)
        recordV(:,i) = V;
        recordVth(:,i) = Vth;
        recordIsyn(:,i) = Isyn;
        recordIapp(:,i) = curIapp;
        if (usePlasticity)
            recordr1(:,i) = r1;
            recordr2(:,i) = r2;
            recordo1(:,i) = o1;
            recordo2(:,i) = o2;
        end
    else
        recordV = V;
        recordVth = Vth;
        recordIsyn = Isyn;
        recordIapp = curIapp;
        if (usePlasticity)
            recordr1 = r1;
            recordr2 = r2;
            recordo1 = o1;
            recordo2 = o2;
        end
    end
    
    % if there are any spikes from SIMULATED neurons, write them to file
    if (areSimSpikes)
        if (useRecord)
            fwrite(spkfid,-1,'int32'); % each timestep with spikes starts with -1
            fwrite(spkfid,i,'int32'); % then the integer # of the timestep
            % then the INDICES FROM THE cells2record VECTOR ARE WRITTEN. NOT THE NEURON IDS THEMSELVES
            fwrite(spkfid,find(spiked(cells2record)),'int32');
        end
    end
end
fprintf('%i simulated spikes ',int32(nSimulatedSpikes));
fprintf('%i generated spikes\n',int32(nGeneratedSpikes));
end
