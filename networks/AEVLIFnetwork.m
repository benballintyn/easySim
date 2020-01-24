classdef AEVLIFnetwork < handle
    
    properties
       groupInfo
       nGroups
       nSpikeGenerators
       nNeurons
       coordinateFrames
       spikeGeneratorInfo
    end
    
    methods
        function obj = AEVLIFnetwork()
            obj.nGroups = 0;
            obj.nSpikeGenerators = 0;
            %obj.groupInfo = struct(
            obj.nNeurons = 0;
            obj.groupInfo = struct('id',{},'name',{},'N',{},'neuronType',{},'isExcitatory',{},'isInhibitory',{},...
                'coordinateFrame',{},'start_ind',{},'end_ind',{},'targets',{},'connections',{},'connectionParams',{},...
                'std_noise',{},'mean_V0',{},'std_V0',{},'mean_Vreset',{},'std_Vreset',{},'mean_Vth0',{},'std_Vth0',{},...
                'mean_Vth_max',{},'std_Vth_max',{},'mean_tau_ref',{},'std_tau_ref',{},...
                'mean_tau_sra',{},'std_tau_sra',{},'mean_a',{},'std_a',{},'mean_b',{},'std_b',{},...
                'mean_VsynE',{},'std_VsynE',{},'mean_VsynI',{},'std_VsynI',{},...
                'mean_max_GsynE',{},'std_max_GsynE',{},'mean_max_GsynI',{},'std_max_GsynI',{},...
                'mean_tau_synE',{},'std_tau_synE',{},'mean_tau_synI',{},'std_tau_synI',{},...
                'mean_Cm',{},'std_Cm',{},'mean_Gl',{},'std_Gl',{},'mean_El',{},...
                'std_El',{},'mean_dth',{},'std_dth',{},'xcoords',{},'ycoords',{},'record',{});
            obj.spikeGeneratorInfo = struct('id',{},'name',{},'N',{},'neuronType',{},'isExcitatory',{},'isInhibitory',{},...
                'firingRate',{},'start_ind',{},'end_ind',{},'targets',{},'connections',{},'connectionParams',{});
            obj.coordinateFrames = struct('ID',{},'xmin',{},'xmax',{},'ymin',{},'ymax',{});
        end
        
        function addGroup(obj,name,N,neuronType,coordinateFrame,varargin)
            % ordered field names
            orderedFieldNames = {'id','name','N','neuronType','isExcitatory','isInhibitory',...
                'coordinateFrame','start_ind','end_ind','targets','connections','connectionParams',...
                'std_noise','mean_V0','std_V0','mean_Vreset','std_Vreset','mean_Vth0','std_Vth0',...
                'mean_Vth_max','std_Vth_max','mean_tau_ref','std_tau_ref',...
                'mean_tau_sra','std_tau_sra','mean_a','std_a','mean_b','std_b',...
                'mean_VsynE','std_VsynE','mean_VsynI','std_VsynI',...
                'mean_max_GsynE','std_max_GsynE','mean_max_GsynI','std_max_GsynI',...
                'mean_tau_synE','std_tau_synE','mean_tau_synI','std_tau_synI',...
                'mean_Cm','std_Cm','mean_Gl','std_Gl','mean_El',...
                'std_El','mean_dth','std_dth','xcoords','ycoords','record'};
            % Create inputParser and assign default values and checks
            p = inputParser;
            validNeuronTypes = {'excitatory','inhibitory'};
            % default values
            defaultXmin = 0;
            defaultXmax = 1;
            defaultYmin = 0;
            defaultYmax = 1;
            default_std_noise = 1e-12;%10e-12; % 50pA*s
            default_mean_V0 = -.07; % -70mV
            default_std_V0 = 0;
            default_mean_Vreset = -.08; % -80mV
            default_std_Vreset = 0;
            default_mean_Vth0 = -.05; % -50mV
            default_std_Vth0 = 0;
            default_mean_Vth_max = .2; % 200mV
            default_std_Vth_max = 0;
            default_mean_tau_ref = 1e-3; % 1 ms
            default_std_tau_ref = 0;
            default_mean_tau_sra = 200e-3; % 200ms
            default_std_tau_sra = 0;
            default_mean_a = 2e-9; % 2nS
            default_std_a = 0;
            default_mean_b = 20e-12; % 20pA
            default_std_b = 0;
            default_mean_VsynE = 0; % 0mV
            default_std_VsynE = 0;
            default_mean_VsynI = -.08; % -80mV
            default_std_VsynI = 0;
            default_mean_max_GsynE = 10e-6; % 10uS
            default_std_max_GsynE = 0;
            default_mean_max_GsynI = 15e-6; % 15uS
            default_std_max_GsynI = 0;
            default_mean_tau_synE = 20e-3; % 10ms
            default_std_tau_synE = 0;
            default_mean_tau_synI = 10e-3; % 1ms
            default_std_tau_synI = 0;
            default_mean_Cm = 10e-9; % 10nF/mm^2%.1e-9; % 100pF
            default_std_Cm = 0;
            default_mean_Gl = 1e-6; % 1/(1Mohm*mm^2) %10e-9; % 10nS
            default_std_Gl = 0;
            default_mean_El = -.07; % -70mV
            default_std_El = 0;
            default_mean_dth = .002; %2mV
            default_std_dth = 0;
            default_record = true;
            checkNeuronType = @(x) any(validatestring(neuronType,validNeuronTypes));
            validNumCheck = @(x) isnumeric(x) && ~isinf(x) && ~isnan(x);
            positiveNoInfCheck = @(x) x > 0 && ~isinf(x);
            nonNegativeNoInfCheck = @(x) x>= 0 && ~isinf(x);
            addRequired(p,'name',@ischar)
            addRequired(p,'N',positiveNoInfCheck);
            addRequired(p,'neuronType',checkNeuronType);
            addRequired(p,'coordinateFrame',positiveNoInfCheck);
            addParameter(p,'xmin',defaultXmin,validNumCheck);
            addParameter(p,'xmax',defaultXmax,validNumCheck);
            addParameter(p,'ymin',defaultYmin,validNumCheck);
            addParameter(p,'ymax',defaultYmax,validNumCheck);
            addParameter(p,'std_noise',default_std_noise,nonNegativeNoInfCheck);
            addParameter(p,'mean_V0',default_mean_V0,validNumCheck);
            addParameter(p,'std_V0',default_std_V0,validNumCheck);
            addParameter(p,'mean_Vreset',default_mean_Vreset,validNumCheck);
            addParameter(p,'std_Vreset',default_std_Vreset',validNumCheck);
            addParameter(p,'mean_Vth0',default_mean_Vth0,validNumCheck);
            addParameter(p,'std_Vth0',default_std_Vth0,validNumCheck);
            addParameter(p,'mean_Vth_max',default_mean_Vth_max,validNumCheck);
            addParameter(p,'std_Vth_max',default_std_Vth_max,validNumCheck);
            addParameter(p,'mean_tau_ref',default_mean_tau_ref,validNumCheck);
            addParameter(p,'std_tau_ref',default_std_tau_ref,validNumCheck);
            addParameter(p,'mean_tau_sra',default_mean_tau_sra,positiveNoInfCheck);
            addParameter(p,'std_tau_sra',default_std_tau_sra,positiveNoInfCheck);
            addParameter(p,'mean_a',default_mean_a,positiveNoInfCheck);
            addParameter(p,'std_a',default_std_a,validNumCheck);
            addParameter(p,'mean_b',default_mean_b,positiveNoInfCheck);
            addParameter(p,'std_b',default_std_b,validNumCheck);
            addParameter(p,'mean_VsynE',default_mean_VsynE,validNumCheck);
            addParameter(p,'std_VsynE',default_std_VsynE,validNumCheck);
            addParameter(p,'mean_VsynI',default_mean_VsynI,validNumCheck);
            addParameter(p,'std_VsynI',default_std_VsynI,validNumCheck);
            addParameter(p,'mean_max_GsynE',default_mean_max_GsynE,validNumCheck);
            addParameter(p,'std_max_GsynE',default_std_max_GsynE,validNumCheck);
            addParameter(p,'mean_max_GsynI',default_mean_max_GsynI,validNumCheck);
            addParameter(p,'std_max_GsynI',default_std_max_GsynI,validNumCheck);
            addParameter(p,'mean_tau_synE',default_mean_tau_synE,validNumCheck);
            addParameter(p,'std_tau_synE',default_std_tau_synE,validNumCheck);
            addParameter(p,'mean_tau_synI',default_mean_tau_synI,validNumCheck);
            addParameter(p,'std_tau_synI',default_std_tau_synI,validNumCheck);
            addParameter(p,'mean_Cm',default_mean_Cm,validNumCheck);
            addParameter(p,'std_Cm',default_std_Cm,validNumCheck);
            addParameter(p,'mean_Gl',default_mean_Gl,validNumCheck);
            addParameter(p,'std_Gl',default_std_Gl,validNumCheck);
            addParameter(p,'mean_El',default_mean_El,validNumCheck);
            addParameter(p,'std_El',default_std_El,validNumCheck);
            addParameter(p,'mean_dth',default_mean_dth,validNumCheck);
            addParameter(p,'std_dth',default_std_dth,validNumCheck);
            addParameter(p,'record',default_record,@islogical)
            parse(p,name,N,neuronType,coordinateFrame,varargin{:})
            
            % Housekeeping to keep track of # of groups and neurons
            obj.nGroups = obj.nGroups + 1;
            groupInfo.id = obj.nGroups;
            groupInfo.start_ind = obj.nNeurons+1;
            groupInfo.end_ind = obj.nNeurons + p.Results.N;
            obj.nNeurons = obj.nNeurons + p.Results.N;
            
            % Funnel parser results into net.groupInfo
            names = fieldnames(p.Results);
            for i=1:length(names)
                if (any(strcmp(names{i},{'xmin','xmax','ymin','ymax'})))
                    continue;
                else
                    groupInfo.(names{i}) = p.Results.(names{i});
                end
            end
            
            % Determine whether the group is excitatory or inhibitory
            groupInfo.isExcitatory = strcmp(p.Results.neuronType,'excitatory');
            groupInfo.isInhibitory = strcmp(p.Results.neuronType,'inhibitory');
            
            % Set coordinate frame of new group and add to network list of
            % coordinate frame if it is new
            cf.ID = p.Results.coordinateFrame;
            cf.xmin = p.Results.xmin;
            cf.xmax = p.Results.xmax;
            cf.ymin = p.Results.ymin;
            cf.ymax = p.Results.ymax;
            if (~ismember(p.Results.coordinateFrame,[obj.coordinateFrames.ID]))
                obj.coordinateFrames = [obj.coordinateFrames cf];
            end
            groupInfo.coordinateFrame = cf;
            
            % assign x and y coordinates to each new neuron
            groupInfo.xcoords = rand(1,p.Results.N)*(p.Results.xmax - p.Results.xmin) + p.Results.xmin;
            groupInfo.ycoords = rand(1,p.Results.N)*(p.Results.ymax - p.Results.ymin) + p.Results.ymin;
            [~,inds] = sort(groupInfo.xcoords);
            groupInfo.xcoords = groupInfo.xcoords(inds);
            groupInfo.ycoords = groupInfo.ycoords(inds);
            
            % empty holders for connections
            groupInfo.targets = [];
            groupInfo.connections = [];
            groupInfo.connectionParams = {};
            
            % reorder field names in the groupInfo for this group
            s = orderfields(groupInfo,orderedFieldNames);
            obj.groupInfo(obj.nGroups) = s;
        end
        
        function addSpikeGenerator(obj,name,N,neuronType,firingRate)
            obj.nSpikeGenerators = obj.nSpikeGenerators + 1;
            info.id = -obj.nSpikeGenerators;
            info.name = name;
            info.N = N;
            info.neuronType = neuronType;
            info.isExcitatory = strcmp(neuronType,'excitatory');
            info.isInhibitory = strcmp(neuronType,'inhibitory');
            info.firingRate = firingRate;
            info.start_ind = []; % to be assigned later
            info.end_ind = []; % to be assigned later
            
            % empty holders for connections
            info.targets = [];
            info.connections = [];
            info.connectionParams = {};
            obj.spikeGeneratorInfo(obj.nSpikeGenerators) = info;
        end
        
        function connect(obj,src_id,tgt_id,connType,connParams)
            if (tgt_id < 0)
                error('Cannot have a SpikeGenerator group as a connection target')
            end
            if (src_id < 0)
                if (AEVLIFnetwork.checkConnInputs(connType,connParams))
                    if (strcmp(connType,'random'))
                        conn=randomConnector(src_id,tgt_id,connParams.connProb,connParams.weightDistribution,obj.groupInfo,obj.spikeGeneratorInfo);
                    elseif (strcmp(connType,'clustered'))
                        error('SpikeGenerators do not support clustered connections')
                    elseif (strcmp(connType,'gaussian'))
                        error('SpikeGenerators do not support gaussian connections')
                    elseif (strcmp(connType,'gradient'))
                        conn=gradientConnector(src_id,tgt_id,connParams.connProbFunction,connParams.weightFunction,obj.groupInfo);
                    end
                    obj.spikeGeneratorInfo(-src_id).targets = [obj.spikeGeneratorInfo(-src_id).targets tgt_id];
                    obj.spikeGeneratorInfo(-src_id).connections = [obj.spikeGeneratorInfo(-src_id).connections conn];
                    obj.spikeGeneratorInfo(-src_id).connectionParams{end+1} = connParams;
                end
            else
                if (AEVLIFnetwork.checkConnInputs(connType,connParams))
                    if (strcmp(connType,'random'))
                        conn=randomConnector(src_id,tgt_id,connParams.connProb,connParams.weightDistribution,obj.groupInfo);
                    elseif (strcmp(connType,'clustered'))
                        if (src_id ~= tgt_id)
                            error('Source and target must be the same for clustered connections')
                        end
                        conn=clusterConnector(src_id,connParams.nClusters,connParams.intraConnProb,connParams.interConnProb,connParams.intraWeightDist,connParams.interWeightDist,obj.groupInfo);
                    elseif (strcmp(connType,'gaussian'))
                        if (obj.groupInfo(src_id).coordinateFrame.ID ~= obj.groupInfo(tgt_id).coordinateFrame.ID)
                            error('Source and target groups must be in the same coordinate frame for gaussian connection')
                        end
                        conn=gaussianConnector(src_id,tgt_id,connParams.connProbFunction,connParams.weightFunction,connParams.useWrap,obj.groupInfo);
                    elseif (strcmp(connType,'gradient'))
                        if (obj.groupInfo(src_id).coordinateFrame.ID == obj.groupInfo(tgt_id).coordinateFrame.ID)
                            error('Source and target groups must be in different coordinate frames for gradient connection')
                        end
                        conn=gradientConnector(src_id,tgt_id,connParams.connProbFunction,connParams.weightFunction,obj.groupInfo);
                    end
                    obj.groupInfo(src_id).targets = [obj.groupInfo(src_id).targets tgt_id];
                    obj.groupInfo(src_id).connections = [obj.groupInfo(src_id).connections conn];
                    obj.groupInfo(src_id).connectionParams{end+1} = connParams;
                end
            end
        end
        
        function print(obj,verbose)
            if (nargin < 2)
                verbose = false;
            end
            if (~verbose)
                fprintf('======================================== %s ======================================== \n',class(obj))
                fprintf('Total # of neurons: %i \n',obj.nNeurons)
                fprintf('Excitatory        : %1$i across groups %2$s \n',sum([obj.groupInfo([obj.groupInfo.isExcitatory]).N]),sprintf('%d ',[obj.groupInfo([obj.groupInfo.isExcitatory]).id]))
                fprintf('Inhibitory        : %1$i across groups %2$s \n',sum([obj.groupInfo([obj.groupInfo.isInhibitory]).N]),sprintf('%d ',[obj.groupInfo([obj.groupInfo.isInhibitory]).id]))
                fprintf('\n')
                fprintf('CONNECTION SUMMARY: \n')
                for i=1:obj.nGroups
                    for j=1:length(obj.groupInfo(i).targets)
                        fprintf('%1$i ====> %2$i %3$s \n',i,obj.groupInfo(i).targets(j),class(obj.groupInfo(i).connections(j)))
                    end
                end
            else
                fprintf('======================================== %s ======================================== \n',class(obj))
                fprintf('Total # of neurons: %i \n',obj.nNeurons)
                fprintf('Excitatory        : %1$i across groups %2$s \n',sum([obj.groupInfo([obj.groupInfo.isExcitatory]).N]),sprintf('%d ',[obj.groupInfo([obj.groupInfo.isExcitatory]).id]))
                fprintf('Inhibitory        : %1$i across groups %2$s \n',sum([obj.groupInfo([obj.groupInfo.isInhibitory]).N]),sprintf('%d ',[obj.groupInfo([obj.groupInfo.isInhibitory]).id]))
                fprintf('\n')
                fprintf('__________________________________________________________________________________________________\n')
                fprintf('GROUP INFO: \n')
                for i=1:obj.nGroups
                    fprintf('GROUP ID: %i \n',obj.groupInfo(i).id)
                    fprintf('GROUP NAME: %s \n', obj.groupInfo(i).name)
                    fprintf('NEURON TYPE: %s \n',obj.groupInfo(i).neuronType)
                    fprintf('NEURON #: %i \n',obj.groupInfo(i).N);
                    fprintf('COORDINATE FRAME: %i \n',obj.groupInfo(i).coordinateFrame.ID)
                    fprintf('     xmin: %1$i     xmax: %2$i \n',obj.groupInfo(i).coordinateFrame.xmin,obj.groupInfo(i).coordinateFrame.xmax)
                    fprintf('     ymin: %1$i     ymax: %2$i \n',obj.groupInfo(i).coordinateFrame.ymin,obj.groupInfo(i).coordinateFrame.ymax)
                    fnames = fieldnames(obj.groupInfo(i));
                    for j=1:length(fnames)
                        if (contains(fnames{j},'mean'))
                            fprintf('%1$15s = %2$10g    %3$15s = %4$10g \n',fnames{j},obj.groupInfo(i).(fnames{j}),fnames{j+1},obj.groupInfo(i).(fnames{j+1}))
                        end
                    end
                    fprintf('____________________________________________\n')
                    fprintf('\n')
                end
                fprintf('CONNECTION SUMMARY: \n')
                for i=1:obj.nGroups
                    for j=1:length(obj.groupInfo(i).targets)
                        fprintf('%1$i ====> %2$i \n',i,obj.groupInfo(i).targets(j))
                        fprintf('    %s \n',class(obj.groupInfo(i).connections(j)))
                        fnames = fieldnames(obj.groupInfo(i).connectionParams{j});
                        for k=1:length(fnames)
                            if (isa(obj.groupInfo(i).connectionParams{j}.(fnames{k}),'function_handle'))
                                fprintf('    %1$s : %2$s \n',fnames{k},func2str(obj.groupInfo(i).connectionParams{j}.(fnames{k})))
                            elseif (isa(obj.groupInfo(i).connectionParams{j}.(fnames{k}),'weightDistribution'))
                                w = obj.groupInfo(i).connectionParams{j}.(fnames{k});
                                fprintf('    %s : \n',fnames{k})
                                fprintf('    min : %1$10g       max : %2$10g \n',w.xrange(1),w.xrange(end))
                            end
                        end
                    end
                end
            end
        end
    end
    
    methods (Static)
        function all_ok=checkConnInputs(connType,connParams)
            if (strcmp(connType,'random'))
                reqFieldNames = {'connProb','weightDistribution'};
            elseif (strcmp(connType,'clustered'))
                reqFieldNames = {'nClusters','intraConnProb','interConnProb','intraWeightDist','interWeightDist'};
            elseif (strcmp(connType,'gaussian'))
                reqFieldNames = {'connProbFunction','weightFunction','useWrap'};
            elseif (strcmp(connType,'gradient'))
                reqFieldNames = {'connProbFunction','weightFunction'};
            end
            presentFields = isfield(connParams,reqFieldNames);
            nMissingFields = sum(~presentFields);
            missingFieldInds = find(presentFields == 0);
            if (nMissingFields > 0)
                all_ok=0;
                if (nMissingFields == 1)
                    error([reqFieldNames{~presentFields} ' is missing from the connectivity parameters'])
                else
                    msg='';
                    for i=1:nMissingFields
                        if (i < nMissingFields)
                            msg = [msg reqFieldNames{missingFieldInds(i)} ', and '];
                        else
                            msg = [msg reqFieldNames{missingFieldInds(i)} ' are missing'];
                        end
                    end
                    error(msg)
                end
            else
                all_ok=1;
            end     
        end
    end
end