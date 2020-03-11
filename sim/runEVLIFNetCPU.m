function [GsynMax,recordV,recordVth,recordIsyn,recordIapp,recordD,recordF] = runEVLIFNetCPU(V,Vreset,...
              tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid,recordVars) %#codegen

N = size(V,1); % # of simulated neurons
nSpikeGen = length(spikeGenProbs); % # of poisson spike generator neurons
useSpikeGen = (nSpikeGen > 0); % determine if there are any spike generators
n2record = length(cells2record); % # of neurons to record
useRecord = (n2record > 0); % determine if any neurons should be recorded
nSimulatedSpikes = 0;
nGeneratedSpikes = 0;
if (~isempty(D))
    useSynDynamics = true;
    Fmax = 1./p0;
else
    useSynDynamics = false;
end
if (~isempty(C))
    usePlasticity = true;
else
    usePlasticity = false;
end
if (recordVars)
    recordV = zeros(N,nT);
    recordVth = zeros(N,nT);
    recordIsyn = zeros(N,nT);
    recordIapp = zeros(N,nT);
    if (useSynDynamics)
        recordD = zeros(N,nT);
        recordF = zeros(N,nT);
    else
        recordD = [];
        recordF = [];
    end
else
    recordV = [];
    recordVth = [];
    recordIsyn = [];
    recordIapp = [];
    recordD = [];
    recordF = [];
end

% if no spike file was given, don't record any spikes
if (spkfid < 0)
    useRecord = false;
end

% Loop through nT timepoints
for i=1:nT
    
    % Update spike thresholds
    vth1 = arrayfun(@minus,Vth0,Vth); % (Vth0 - Vth)
    dVthdt = arrayfun(@rdivide,vth1,tau_ref); % (Vth0 - Vth)./tau_ref
    Vth = arrayfun(@plus,Vth,dVthdt*dt); % Vth = Vth + dt*[(Vth0 - Vth)./tau_ref]
    
    if (useSynDynamics)
        % update depression/facilitation variables
        dDdt = arrayfun(@rdivide,(1 - D),tau_D);
        dFdt = arrayfun(@rdivide,(1 - F),tau_F);
        D = arrayfun(@plus,D,dDdt*dt);
        F = arrayfun(@plus,F,dFdt*dt);
    end
    
    if (usePlasticity)
        % decay plasticity variables
        dr1dt = arrayfun(@rdivide,-r1,tau_plus);
        dr2dt = arrayfun(@rdivide,-r2,tau_x);
        do1dt = arrayfun(@rdivide,-o1,tau_minus);
        do2dt = arrayfun(@rdivide,-o2,tau_y);
        r1 = arrayfun(@plus,r1,dr1dt*dt);
        r2 = arrayfun(@plus,r2,dr2dt*dt);
        o1 = arrayfun(@plus,o1,do1dt*dt);
        o2 = arrayfun(@plus,o2,do2dt*dt);
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
    end
    
    e_spiked = logical(allSpikes.*ecells); % all excitatory simulated and spike generating neurons that spiked
    i_spiked = logical(allSpikes.*icells); % all inhibitory simulated and spike generating neurons that spiked
    dGsynEdt = arrayfun(@rdivide,-GsynE,tau_synE); % decay in excitatory synaptic conductances
    dGsynIdt = arrayfun(@rdivide,-GsynI,tau_synI); % decay in inhibitory synaptic conductances
    
    GsynE = arrayfun(@plus,GsynE,dGsynEdt*dt);
    GsynI = arrayfun(@plus,GsynI,dGsynIdt*dt);

    if (areAnySpikes)
        if (useSynDynamics)
            % multiply release probability with maximum synaptic conductance
            dGsynE1 = bsxfun(@times,GsynMax(:,e_spiked),p0(e_spiked)');
            dGsynE2 = bsxfun(@times,dGsynE1,F(e_spiked)');
            dGsynE3 = bsxfun(@times,dGsynE2,D(e_spiked)');
            dGsynI1 = bsxfun(@times,GsynMax(:,i_spiked),p0(i_spiked)');
            dGsynI2 = bsxfun(@times,dGsynI1,F(i_spiked)');
            dGsynI3 = bsxfun(@times,dGsynI2,D(i_spiked)');
            % sum across all inputs
            dGsynE_sum = sum(dGsynE3,2);
            dGsynI_sum = sum(dGsynI3,2);
            
            % update depression/facilitation variables for neurons that spiked
            d1 = arrayfun(@times,p0(allSpikes),F(allSpikes));
            d2 = arrayfun(@times,d1,D(allSpikes));
            d3 = arrayfun(@times,d2,has_depression);
            D(allSpikes) = arrayfun(@minus,D(allSpikes),d3);
            
            f1 = arrayfun(@minus,Fmax(allSpikes),F(allSpikes));
            f2 = arrayfun(@times,f_fac(allSpikes),f1);
            f3 = arrayfun(@times,f2,has_facilitation);
            F(allSpikes) = arrayfun(@plus,F(allSpikes),f3);
        else
            dGsynE = bsxfun(@times,GsynMax(:,e_spiked),p0(e_spiked)');
            dGsynI = bsxfun(@times,GsynMax(:,i_spiked),p0(i_spiked)');
            dGsynE_sum = sum(dGsynE,2);
            dGsynI_sum = sum(dGsynI,2);
        end
        GsynE = arrayfun(@plus,GsynE,dGsynE_sum);
        GsynI = arrayfun(@plus,GsynI,dGsynI_sum);
    end
    
    % Compute total synaptic current for each neuron
    Isyn = arrayfun(@plus,arrayfun(@times,GsynE,arrayfun(@minus,VsynE,V)),arrayfun(@times,GsynI,arrayfun(@minus,VsynI,V)));
    
    % add noise to any input currents
    curIapp = Iapp;
    curIapp = arrayfun(@plus,curIapp,arrayfun(@times,std_noise,randn(N,1)));
    
    % update membrane voltages
    f1 = (1./Cm); % (1./Cm)
    f2 = arrayfun(@minus,El,V); % El - V
    f3 = arrayfun(@minus,V,Vth); % V - Vth
    f4 = arrayfun(@rdivide,f3,dth); % (V - Vth)/dth
    f5 = arrayfun(@exp,f4); % exp((V - Vth)/dth)
    f6 = arrayfun(@times,dth,f5); % dth.*exp((V - Vth)/dth)
    f7 = arrayfun(@plus,f2,f6); % (El - V) + dth.*exp((V - Vth)/dth)
    f8 = arrayfun(@times,Gl,f7); % Gl.*[(El - V) + dth.*exp((V - Vth)/dth)]
    f9 = arrayfun(@plus,f8,Isyn); % Gl.*[(El - V) + dth.*exp((V - Vth)/dth)] + Isyn
    f10 = arrayfun(@plus,f9,curIapp); % Gl.*[(El - V) + dth.*exp((V - Vth)/dth)] + Isyn + Iapp
    dVdt = arrayfun(@times,f1,f10); % (1./Cm).*(Gl.*[(El - V) + dth.*exp((V - Vth)/dth)] + Isyn + Iapp)
    
    V = arrayfun(@plus,V,dVdt*dt); % V = V + dVdt*dt
    V = arrayfun(@max,Vreset,V); % bound membrane potentials to be >= than the reset value
    
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
        recordIapp(:,i) = Iapp;
        if (useSynDynamics)
            recordD(:,i) = D;
            recordF(:,i) = F;
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
fprintf('%i simulated spikes\n',int32(nSimulatedSpikes));
fprintf('%i generated spikes\n',int32(nGeneratedSpikes));
end