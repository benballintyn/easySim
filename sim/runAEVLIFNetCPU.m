function [GsynMax,recordV,recordVth,recordIsra,recordGsynE,recordGsynI,recordD,recordF,recordr1,recordr2,recordo1,recordo2,recordIsyn,recordIapp] = runAEVLIFNetCPU(V_temp,Vreset,tau_ref,Vth_temp,Vth0,Vth_max,...
              Isra_temp,tau_sra,a,b,VsynE,VsynI,GsynE_temp,GsynI_temp,GsynMax_temp,tau_D,tau_F,f_fac,D_temp,F_temp,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1_temp,r2_temp,o1_temp,o2_temp,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid,recordVars) %#codegen

GsynMax = GsynMax_temp;
V = V_temp;
Vth = Vth_temp;
Isra = Isra_temp;
GsynE = GsynE_temp;
GsynI = GsynI_temp;
D = D_temp;
F = F_temp;
r1 = r1_temp;
r2 = r2_temp;
o1 = o1_temp;
o2 = o2_temp;          

N = size(V,1); % # of simulated neurons
nSpikeGen = length(spikeGenProbs); % # of poisson spike generator neurons
useSpikeGen = (nSpikeGen > 0); % determine if there are any spike generators
n2record = length(cells2record); % # of neurons to record
useRecord = (n2record > 0); % determine if any neurons should be recorded
nSimulatedSpikes = 0;
nGeneratedSpikes = 0;
Fmax = 1./p0;
if (~isnan(D))
    useSynDynamics = true;
else
    useSynDynamics = false;
end
if (~isnan(C))
    usePlasticity = true;
else
    usePlasticity = false;
end
if (recordVars)
    recordV = zeros(N,nT);
    recordVth = zeros(N,nT);
    recordIsra = zeros(N,nT);
    recordGsynE = zeros(N,nT);
    recordGsynI = zeros(N,nT);
    recordIsyn = zeros(N,nT);
    recordIapp = zeros(N,nT);
    if (useSynDynamics)
        recordD = zeros(N,nT);
        recordF = zeros(N,nT);
    else
        recordD = nan;
        recordF = nan;
    end
    if (usePlasticity)
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
    recordIsra = zeros(N,1);
    recordGsynE = zeros(N,1);
    recordGsynI = zeros(N,1);
    recordIsyn = zeros(N,1);
    recordIapp = zeros(N,1);
    if (useSynDynamics)
        recordD = zeros(N,1);
        recordF = zeros(N,1);
    else
        recordD = nan;
        recordF = nan;
    end
    if (usePlasticity)
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
    Vth = Vth + dt*((Vth0 - Vth)./tau_ref);
    
    % Update adaptation currents
    Isra = Isra + dt*((a.*(V - El) - Isra)./tau_sra);
    
    if (useSynDynamics)
        % update depression/facilitation variables
        D = D + dt*((1-D)./tau_D);
        F = F + dt*((1-F)./tau_F);
    end
    
    if (usePlasticity)
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
        Vth(spiked) = Vth_max(spiked); % set spike threshold to max
        Isra(spiked) = Isra(spiked) + b(spiked);
    end
    
    e_spiked = logical(allSpikes.*ecells); % all excitatory simulated and spike generating neurons that spiked
    i_spiked = logical(allSpikes.*icells); % all inhibitory simulated and spike generating neurons that spiked
    
    GsynE = GsynE + dt*(-GsynE./tau_synE);
    GsynI = GsynI + dt*(-GsynI./tau_synI);

    if (areAnySpikes)
        if (useSynDynamics)
            % multiply release probability with maximum synaptic conductance
            dGsynE1 = p0(e_spiked)'.*F(e_spiked)'.*D(e_spiked)';
            dGsynE2 = bsxfun(@times,GsynMax(:,e_spiked),dGsynE1);
            
            dGsynI1 = p0(i_spiked)'.*F(i_spiked)'.*D(i_spiked)';
            dGsynI2 = bsxfun(@times,GsynMax(:,i_spiked),dGsynI1);
            
            % sum across all inputs
            dGsynE_sum = sum(dGsynE2,2);
            dGsynI_sum = sum(dGsynI2,2);
            
            % update depression/facilitation variables for neurons that spiked
            D(allSpikes) = D(allSpikes) - p0(allSpikes).*F(allSpikes).*D(allSpikes).*has_depression(allSpikes);
            F(allSpikes) = F(allSpikes) + (Fmax(allSpikes) - F(allSpikes)).*f_fac(allSpikes).*has_facilitation(allSpikes);
        else
            dGsynE = bsxfun(@times,GsynMax(:,e_spiked),p0(e_spiked)');
            dGsynI = bsxfun(@times,GsynMax(:,i_spiked),p0(i_spiked)');
            dGsynE_sum = sum(dGsynE,2);
            dGsynI_sum = sum(dGsynI,2);
        end
        GsynE = GsynE + dGsynE_sum;
        GsynI = GsynI + dGsynI_sum;
    end
    
    % Compute total synaptic current for each neuron
    Isyn = GsynE.*(VsynE - V) + GsynI.*(VsynI - V);
    
    % add noise to any input currents
    curIapp = Iapp;
    curIapp = curIapp + std_noise.*randn(N,1);
    
    % update membrane voltages
    V = V + dt*(1./Cm).*(Gl.*((El - V) + dth.*exp((V - Vth)./dth)) + Isyn + curIapp - Isra);
    V = max(V,Vreset);
    
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
        recordIsra(:,i) = Isra;
        recordGsynE(:,i) = GsynE;
        recordGsynI(:,i) = GsynI;
        recordIsyn(:,i) = Isyn;
        recordIapp(:,i) = curIapp;
        if (useSynDynamics)
            recordD(:,i) = D;
            recordF(:,i) = F;
        end
        if (usePlasticity)
            recordr1(:,i) = r1;
            recordr2(:,i) = r2;
            recordo1(:,i) = o1;
            recordo2(:,i) = o2;
        end
    else
        recordV = V;
        recordVth = Vth;
        recordIsra = Isra;
        recordGsynE = GsynE;
        recordGsynI = GsynI;
        recordIsyn = Isyn;
        recordIapp = curIapp;
        if (useSynDynamics)
            recordD = D;
            recordF = F;
        end
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