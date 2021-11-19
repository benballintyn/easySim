classdef AEVLIFnetwork < handle
    
    properties
       groupInfo
       nGroups
       nSpikeGenerators
       nNeurons
       nSpikeGeneratorNeurons
       coordinateFrames
       spikeGeneratorInfo
       is_plastic
       is_dynamic
    end
    
    methods
        function obj = LSNN()
            % This class holds all of the necessary information about an
            % LSNN(Bellec et al., 2020). These networks consist of a mixture
	    % of standard LIF neurons and LIF neurons with an adaptive threshold
	    % (ALIF neurons). This type of network is meant to be paired with the
	    % e-prop learning algorithm.
            %
            % METHODS:
            %   addGroup(name,N,neuronType,coordinateFrame,varargin)
            %       - add a neuron group to be simulated.
            %         Type help EVLIFnetwork.addGroup for more info
            %
            %   addSpikeGenerator(name,N,neuronType,firingRate)
            %       - add a Poisson spike generator group.
            %         Type help EVLIFnetwork.addSpikeGenerator for more
            %         info
            %
            %   connect(src_id,tgt_id,connType,connParams)
            %       - connect 2 neuron groups using a connType connection.
            %         This can be used to connect a neuron group to itself
            %         or to another group. It can also be used to connect a
            %         spikeGenerator group to a neuron group. Check the
            %         different connectionTypes to see which can be used in
            %         which situation.
            %         Type help EVLIFnetwork.connect for more info
            %
            %   print(verbose)
            %       - This prints out information about the network.
            %         Calling net.print() or net.print(false) will print a
            %         short description of the network. Calling
            %         net.print(true) will print a much longer description
            %         of the network with each group's parameter sets.
            
            obj.nGroups = 0;
            obj.nSpikeGenerators = 0;
            obj.nNeurons = 0;
            obj.nSpikeGeneratorNeurons = 0;
            obj.groupInfo = struct('id',{},'name',{},'N',{},'neuronType',{},'isExcitatory',{},'isInhibitory',{},...
                'coordinateFrame',{},'start_ind',{},'end_ind',{},'targets',{},'connections',{},'connectionParams',{},...
                'std_noise',{},'mean_V0',{},'std_V0',{},'mean_Vreset',{},'std_Vreset',{},'mean_Vth0',{},'std_Vth0',{},...
                'mean_t_refrac',{},'std_t_refrac',{},'mean_tau_m',{},'std_tau_m',{},...
                'xcoords',{},'ycoords',{},'record',{});
            obj.spikeGeneratorInfo = struct('id',{},'name',{},'N',{},'neuronType',{},'isExcitatory',{},'isInhibitory',{},...
                'firingRate',{},'start_ind',{},'end_ind',{},'targets',{},'connections',{},'connectionParams',{});
            obj.coordinateFrames = struct('ID',{},'xmin',{},'xmax',{},'ymin',{},'ymax',{});
            obj.is_plastic = false;
        end
        
        function addGroup(obj,name,N,neuronType,coordinateFrame,varargin)
            % addGroup(name,N,neuronType,coordinateFrame,varargin)
            %   name            - string name of the group. Note, this string
            %                     is not used internally to identify the group
            %                     so it may be whatever the user likes.
            %
            %   N               - # of neurons in the group
            %
            %   neuronType      - either 'excitatory' or 'inhibitory'
            %
            %   coordinateFrame - ID # of the coordinate frame in which
            %                     this neuron group will reside
            %
            %   varargin        - there are many optional parameters that
            %                     the user may supply to customize the 
            %                     neuron group. They must be provided as 
            %                     key-value pairs e.g.
            %
            %     addGroup(name,N,neuronType,coordinateFrame,'mean_V0',-.06)
            %       where 'mean_V0' is the parameter for the mean initial 
            %       membrane potential and -.06 (-60mV) is the provided
            %       value. Each parameter governing a neuron group's
            %       behavior is given a mean and standard deviation from
            %       which the values for individual neurons will be chosen.
            %     
            %     Full list of optional parameters:
            %       1. 'std_noise'      - standard deviation of noise currents
            %       2. 'mean_V0'        - mean initial membrane voltage value
            %       3. 'std_V0'         - standard deviation of the initial
            %                             membrane voltage values
            %       4. 'mean_Vreset'    - mean membrane reset voltage
            %                             following a spike
            %       5. 'std_Vreset'     - standard deviation of the membrane
            %                             reset voltages following a spike
            %       6. 'mean_Vth0'      - mean baseline spike threshold
            %       7. 'std_Vth0'       - standard deviation of baseline spike
            %                             threshold
	    %       8. 'mean_t_refrac'  - mean time (across cells) membrane 
	    %				  voltage is held at 0 following a spike
	    %       9. 'std_t_refrac'   - standard deviation of refractory time
	    %  				  (across cells)
	    %	    10. 'mean_tau_m'    - mean (across cells) membrane time constant
	    %	    11. 'std_tau_m'     - standard deviation (across cells) of the
	    %				  membrane time constant
            %      54. 'record'         - can be true or false. Indicates
            %                             whether to record this groups spike
            %                             times to a file or not
            %      55. 'xmin'           - minimum x-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      56. 'xmax'           - maximum x-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      57. 'ymin'           - minimum y-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      58. 'ymax'           - maximum y-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            
            % ordered field names
            orderedFieldNames = {'id','name','N','neuronType','isExcitatory','isInhibitory',...
                'coordinateFrame','start_ind','end_ind','targets','connections','connectionParams',...
                'std_noise','mean_V0','std_V0','mean_Vreset','std_Vreset','mean_Vth0','std_Vth0',...
                'mean_t_refrac','std_t_refrac','mean_tau_m','std_tau_m','xcoords','ycoords','record'};
            % Create inputParser and assign default values and checks
            p = inputParser;
            validNeuronTypes = {'excitatory','inhibitory'};
            % default values
            defaultXmin = 0;
            defaultXmax = 1;
            defaultYmin = 0;
            defaultYmax = 1;
            default_std_noise = 50e-12; % 50pA*s
            default_mean_V0 = -.07; % -70mV
            default_std_V0 = 0;
            default_mean_Vreset = -.08; % -80mV
            default_std_Vreset = 0;
            default_mean_Vth0 = -.05; % -50mV
            default_std_Vth0 = 0;
	    default_mean_t_refrac = 2e-3; % 2ms
	    default_std_t_refrac = 0;
	    default_mean_tau_m = 20e-3; % 20ms
	    default_std_tau_m = 0;
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
	    addParameter(p,'mean_t_refrac',default_mean_t_refrac,validNumCheck);
	    addParameter(p,'std_t_refrac',default_std_t_refrac,validNumCheck);
	    addParameter(p,'mean_tau_m',default_mean_tau_m,validNumCheck);
	    addParameter(p,'std_tau_m',default_std_tau_m,validNumCheck);
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
        
        function addSpikeGenerator(obj,name,N,neuronType,firingRate,varargin)
            % addSpikeGenerator(name,N,neuronType,firingRate)
            %   name       - string name for this spikeGenerator group
            %
            %   N          - # of neurons in this spikeGenerator group
            %
            %   neuronType - 'excitatory' or 'inhibitory'
            %
            %   firingRate - firing rate of this group in Hz
            %
            %   varargin:
            p = inputParser;
            validNeuronTypes = {'excitatory','inhibitory'};
            checkNeuronType = @(x) any(validatestring(neuronType,validNeuronTypes));
            positiveNoInfCheck = @(x) x > 0 && ~isinf(x);
            nonNegativeNoInfCheck = @(x) x>= 0 && ~isinf(x);
            addRequired(p,'name',@ischar);
            addRequired(p,'N',positiveNoInfCheck);
            addRequired(p,'neuronType',checkNeuronType);
            addRequired(p,'firingRate',nonNegativeNoInfCheck);
            parse(p,name,N,neuronType,firingRate,varargin{:})
            
            obj.nSpikeGenerators = obj.nSpikeGenerators + 1;
            obj.nSpikeGeneratorNeurons = obj.nSpikeGeneratorNeurons + N;
            info.name = p.Results.name;
            info.id = -obj.nSpikeGenerators;
            info.isExcitatory = strcmp(neuronType,'excitatory');
            info.isInhibitory = strcmp(neuronType,'inhibitory');
            info.firingRate = p.Results.firingRate;
            info.start_ind = []; % to be assigned later
            info.end_ind = []; % to be assigned later
            
            names = fieldnames(p.Results);
            for i=1:length(names)
                if (any(strcmp(names{i},{'xmin','xmax','ymin','ymax'})))
                    continue;
                else
                    info.(names{i}) = p.Results.(names{i});
                end
            end
            
            % empty holders for connections
            info.targets = [];
            info.connections = [];
            info.connectionParams = {};
            obj.spikeGeneratorInfo(obj.nSpikeGenerators) = info;
        end
        
        function connect(obj,src_id,tgt_id,connType,connParams)
            % connect(src_id,tgt_id,connType,connParams)
            %   src_id  - ID # of presynaptic group
            %
            %   tgt_id  - ID # of postsynaptic group
            %
            %   connType - 'random', 'clustered', 'gaussian', or 'gradient'
            %
            %   connParams - structure containing the necessary inputs to
            %                build the connection object. Different
            %                connection types require different inputs.
            %                These are listed below.
            %       General connParams:
            %           is_plastic - true if this connection is subject to
            %                        synaptic plasticity (triplet STDP)
            %
            %       connParams by connType:
            %       'random':
            %           connParams.connProb - connection probability
            %           connParams.weightDistribution - weightDistribution
            %                                           object
            %       'clustered':
            %           connParams.nClusters - # of clusters to split the
            %                                  neuron group into
            %           connParams.intraConnProb - connection probability
            %                                      within a cluster
            %           connParams.interConnProb - connection probability
            %                                      between clusters
            %           connParams.intraWeightDist - weightDistribution
            %                                        object for
            %                                        intra-cluster
            %                                        connections
            %           connParams.interWeightDist - weightDistribution
            %                                        object for
            %                                        inter-cluster
            %                                        connections
            %       'gaussian':
            %           connParams.connProbFunction - function that takes
            %                                         as input a distance
            %                                         and returns a
            %                                         connection
            %                                         probability
            %           connParams.weightFunction - function that takes as
            %                                       input a distance and
            %                                       returns a synaptic
            %                                       weight
            %           connParams.useWrap - true or false. If true, when
            %                                calculating the distances
            %                                between neurons, the wrap
            %                                distance (distance as if the
            %                                coordinate frame wrapped
            %                                around on itself) is used.
            %       'gradient':
            %           connParams.connProbFunction - function that takes
            %                                         as input the
            %                                         post-synaptic
            %                                         x-coordinate
            %                                         and returns a
            %                                         connection
            %                                         probability
            %           connParams.weightFunction - function that takes as
            %                                       input the post-synaptic
            %                                       x-coordinate and
            %                                       returns a synaptic
            %                                       weight
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
                        conn=gradientConnector(src_id,tgt_id,connParams.connProbFunction,connParams.weightFunction,obj.groupInfo,obj.spikeGeneratorInfo);
                    end
                    obj.spikeGeneratorInfo(-src_id).targets = [obj.spikeGeneratorInfo(-src_id).targets tgt_id];
                    obj.spikeGeneratorInfo(-src_id).connections = [obj.spikeGeneratorInfo(-src_id).connections conn];
                    if (isfield(connParams,'is_plastic'))
                        if (connParams.is_plastic)
                            obj.is_plastic = true;
                        end
                    else
                        connParams.is_plastic = false;
                    end
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
                    if (isfield(connParams,'is_plastic'))
                        if (connParams.is_plastic)
                            obj.is_plastic = true;
                        end
                    else
                        connParams.is_plastic = false;
                    end
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
                fprintf('Total # of simulated neurons: %i \n',obj.nNeurons)
                fprintf('Excitatory        : %1$i across groups %2$s \n',sum([obj.groupInfo([obj.groupInfo.isExcitatory]).N]),sprintf('%d ',[obj.groupInfo([obj.groupInfo.isExcitatory]).id]))
                fprintf('Inhibitory        : %1$i across groups %2$s \n',sum([obj.groupInfo([obj.groupInfo.isInhibitory]).N]),sprintf('%d ',[obj.groupInfo([obj.groupInfo.isInhibitory]).id]))
                fprintf('\n')
                fprintf('Total # of SpikeGenerator neurons: %i \n',obj.nSpikeGeneratorNeurons)
                fprintf('Excitatory        : %1$i across generators %2$s \n',sum([obj.spikeGeneratorInfo([obj.spikeGeneratorInfo.isExcitatory]).N]),sprintf('%d ',[obj.spikeGeneratorInfo([obj.spikeGeneratorInfo.isExcitatory]).id]))
                fprintf('Inhibitory        : %1$i across generators %2$s \n',sum([obj.spikeGeneratorInfo([obj.spikeGeneratorInfo.isInhibitory]).N]),sprintf('%d ',[obj.spikeGeneratorInfo([obj.spikeGeneratorInfo.isInhibitory]).id]))
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
                        elseif (contains(fnames{j},'noise'))
                            fprintf('%1$15s = %2$10g \n',fnames{j},obj.groupInfo(i).(fnames{j}));
                        end
                    end
                    fprintf('________________________________________________________________\n')
                    fprintf('\n')
                end
                fprintf('__________________________________________________________________________________________________\n')
                fprintf('SPIKE GENERATOR INFO: \n')
                for i=1:obj.nSpikeGenerators
                    fprintf('GROUP ID: %i \n',obj.spikeGeneratorInfo(i).id)
                    fprintf('GROUP NAME: %s \n', obj.spikeGeneratorInfo(i).name)
                    fprintf('NEURON TYPE: %s \n', obj.spikeGeneratorInfo(i).neuronType)
                    fprintf('NEURON #: %i \n', obj.spikeGeneratorInfo(i).N)
                    fprintf('INDICES: %1$i ---> %2$i \n',obj.spikeGeneratorInfo(i).start_ind,obj.spikeGeneratorInfo(i).end_ind)
                    fprintf('FIRING RATE: %f \n', obj.spikeGeneratorInfo(i).firingRate)
                    fprintf('________________________________________________________________\n')
                    fprintf('\n')
                end
                fprintf('__________________________________________________________________________________________________\n')
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
                for i=1:obj.nSpikeGenerators
                    for j=1:length(obj.spikeGeneratorInfo(i).targets)
                        fprintf('%1$i ====> %2$i \n',obj.spikeGeneratorInfo(i).id,obj.spikeGeneratorInfo(i).targets(j))
                        fprintf('    %s \n',class(obj.spikeGeneratorInfo(i).connections(j)))
                        fnames = fieldnames(obj.spikeGeneratorInfo(i).connectionParams{j});
                        for k=1:length(fnames)
                            if (isa(obj.spikeGeneratorInfo(i).connectionParams{j}.(fnames{k}),'function_handle'))
                                fprintf('    %1$s : %2$s \n',fnames{k},func2str(obj.spikeGeneratorInfo(i).connectionParams{j}.(fnames{k})))
                            elseif (isa(obj.spikeGeneratorInfo(i).connectionParams{j}.(fnames{k}),'weightDistribution'))
                                w = obj.spikeGeneratorInfo(i).connectionParams{j}.(fnames{k});
                                fprintf('    %s : \n',fnames{k})
                                fprintf('    min : %1$10g       max : %2$10g \n',w.xrange(1),w.xrange(end))
                            end
                        end
                    end
                end
            end
        end
        
        function ecells = getEcells(self)
            for i=1:self.nGroups
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
