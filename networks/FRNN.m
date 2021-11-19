classdef FRNN < handle
    
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
        function obj = FRNN()
            % This class holds all of the necessary information about an
            % FRNN (firing rate neural network).
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
            obj.groupInfo = struct('id',{},'name',{},'N',{},'transferFunction',{},...
                'coordinateFrame',{},'depressed_synapses',{},'facilitating_synapses',{},'start_ind',{},'end_ind',{},...
                'targets',{},'connections',{},'connectionParams',{},...
                'std_noise',{},'mean_tau_D',{},'std_tau_D',{},'mean_tau_F',{},...
                'std_tau_F',{},'mean_p0',{},'std_p0',{},'mean_f_fac',{},'std_f_fac',{},...
                'mean_A2plus',{},'std_A2plus',{},'mean_A3plus',{},'std_A3plus',{},...
                'mean_A2minus',{},'std_A2minus',{}, 'mean_A3minus',{},'std_A3minus',{},...
                'mean_tau_plus',{},'std_tau_plus',{},'mean_tau_x',{},'std_tau_x',{},...
                'mean_tau_minus',{},'std_tau_minus',{},'mean_tau_y',{},'std_tau_y',{},...
                'xcoords',{},'ycoords',{},'record',{});
            obj.spikeGeneratorInfo = struct('id',{},'name',{},'N',{},'neuronType',{},...
                'isExcitatory',{},'isInhibitory',{},'firingRate',{},'facilitating_synapses',{},...
                'depressed_synapses',{},'mean_tau_D',{},'std_tau_D',{},'mean_tau_F',{},...
                'std_tau_F',{},'mean_p0',{},'std_p0',{},'mean_f_fac',{},'std_f_fac',{},...
                'mean_A2plus',{},'std_A2plus',{},'mean_A3plus',{},'std_A3plus',{},...
                'mean_tau_plus',{},'std_tau_plus',{},'mean_tau_x',{},'std_tau_x',{},...
                'start_ind',{},'end_ind',{},'targets',{},'connections',{},'connectionParams',{});
            obj.coordinateFrames = struct('ID',{},'xmin',{},'xmax',{},'ymin',{},'ymax',{});
            obj.is_plastic = false;
            obj.is_dynamic = false;
        end
        
        function addGroup(obj,name,N,transferFunction,coordinateFrame,varargin)
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
            %      31. 'mean_A2plus'    - mean doublet LTP factor for
            %                             triplet STDP
            %      32. 'std_A2plus'     - standard deviation of the doublet
            %                             LTP factor.
            %      33. 'mean_A3plus'    - mean triplet LTP factor
            %      34. 'std_A3plus'     - standard deviation of the triplet
            %                             LTP factor
            %      35. 'mean_A2minus'   - mean doublet LTD factor
            %      36. 'std_A2minus'    - standard deviation of the doublet
            %                             LTD factor
            %      37. 'mean_A3minus'   - mean triplet LTD factor
            %      38. 'std_A3minus'    - standard deviation of the triplet
            %                             LTD factor
            %      39. 'mean_tau_plus'  - mean decay time constant of the r1
            %                             (presynaptic) STDP variable
            %      40. 'std_tau_plus'   - standard deviation of the r1
            %                             decay time constant
            %      41. 'mean_tau_x'     - mean decay time constant of the
            %                             r2 (presynaptic) STDP variable
            %      42. 'std_tau_x'      - standard deviation of the r2
            %                             decay time constant
            %      43. 'mean_tau_minus' - mean decay time constant of the
            %                             o1 (postsynaptic) STDP variable
            %      44. 'std_tau_minus'  - standard deviation of the o1
            %                             decay time constant
            %      45. 'mean_tau_y'     - mean decay time constant of the
            %                             o2 (postsynaptic) STDP variable
            %      46. 'std_tau_y'      - standard deviation of the o2
            %                             decay time constant
            %      47. 'record'         - can be true or false. Indicates
            %                             whether to record this groups spike
            %                             times to a file or not
            %      48. 'xmin'           - minimum x-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      49. 'xmax'           - maximum x-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      50. 'ymin'           - minimum y-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      51. 'ymax'           - maximum y-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            
            % ordered field names
            orderedFieldNames = {'id','name','N','transferFunction',
                'coordinateFrame','depressed_synapses','facilitating_synapses','start_ind','end_ind',...
                'targets','connections','connectionParams',...
                'std_noise','mean_tau_D','std_tau_D','mean_tau_F','std_tau_F',...
                'mean_p0','std_p0','mean_f_fac','std_f_fac','mean_A2plus','std_A2plus',...
                'mean_A3plus','std_A3plus','mean_A2minus','std_A2minus',...
                'mean_A3minus','std_A3minus','mean_tau_plus','std_tau_plus',...
                'mean_tau_x','std_tau_x','mean_tau_minus','std_tau_minus',...
                'mean_tau_y','std_tau_y','xcoords','ycoords','record'};
            % Create inputParser and assign default values and checks
            p = inputParser;
            validTransferFunctions = {'linear','ReLU','tanh','sigmoid'};
            % default values
            defaultXmin = 0;
            defaultXmax = 1;
            defaultYmin = 0;
            defaultYmax = 1;
            default_std_noise = 50e-12;% 1pA*s
            default_mean_tau_D = .25; % 250ms (Paul Miller Textbook)
            default_std_tau_D = 0;
            default_mean_tau_F = .25; % 250ms (Paul Miller Textbook)
            default_std_tau_F = 0;
            default_mean_p0 = .4; % (Biro et al., 2005)
            default_std_p0 = 0;
            default_mean_f_fac = .25; % (Paul Miller Textbook)
            default_std_f_fac = 0;
            default_mean_A2plus = 8.8e-11; % Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_A2plus = 0;
            default_mean_A3plus = 5.3e-2; % Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_A3plus = 0;
            default_mean_A2minus = 6.6e-3; % Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_A2minus = 0;
            default_mean_A3minus = 3.1e-3; % Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_A3minus = 0;
            default_mean_tau_plus = .0168; % 16.8ms (Bi and Poo, 2001)
            default_std_tau_plus = 0;
            default_mean_tau_x = .714; % 714ms Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_tau_x = 0;
            default_mean_tau_minus = .0337; % 33.7ms (Bi and Poo, 2001)
            default_std_tau_minus = 0;
            default_mean_tau_y = .04; % 40ms Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_tau_y = 0;
            default_record = true;
            checkTransferFunction = @(x) any(validatestring(transferFunction,validTransferFunctions));
            validNumCheck = @(x) isnumeric(x) && ~isinf(x) && ~isnan(x);
            positiveNoInfCheck = @(x) x > 0 && ~isinf(x);
            nonNegativeNoInfCheck = @(x) x>= 0 && ~isinf(x);
            addRequired(p,'name',@ischar)
            addRequired(p,'N',positiveNoInfCheck);
            addRequired(p,'transferFunction',checkTransferFunction);
            addRequired(p,'coordinateFrame',positiveNoInfCheck);
            addParameter(p,'depressed_synapses',false,@islogical);
            addParameter(p,'facilitating_synapses',false,@islogical);
            addParameter(p,'xmin',defaultXmin,validNumCheck);
            addParameter(p,'xmax',defaultXmax,validNumCheck);
            addParameter(p,'ymin',defaultYmin,validNumCheck);
            addParameter(p,'ymax',defaultYmax,validNumCheck);
            addParameter(p,'std_noise',default_std_noise,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_D',default_mean_tau_D,positiveNoInfCheck);
            addParameter(p,'std_tau_D',default_std_tau_D,positiveNoInfCheck);
            addParameter(p,'mean_tau_F',default_mean_tau_F,positiveNoInfCheck);
            addParameter(p,'std_tau_F',default_std_tau_F,positiveNoInfCheck);
            addParameter(p,'mean_p0',default_mean_p0,nonNegativeNoInfCheck);
            addParameter(p,'std_p0',default_std_p0,positiveNoInfCheck);
            addParameter(p,'mean_f_fac',default_mean_f_fac,nonNegativeNoInfCheck);
            addParameter(p,'std_f_fac',default_std_f_fac,positiveNoInfCheck);
            addParameter(p,'mean_A2plus',default_mean_A2plus,nonNegativeNoInfCheck);
            addParameter(p,'std_A2plus',default_std_A2plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_A3plus',default_mean_A3plus,nonNegativeNoInfCheck);
            addParameter(p,'std_A3plus',default_std_A3plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_A2minus',default_mean_A2minus,nonNegativeNoInfCheck);
            addParameter(p,'std_A2minus',default_std_A2minus,nonNegativeNoInfCheck);
            addParameter(p,'mean_A3minus',default_mean_A3minus,nonNegativeNoInfCheck);
            addParameter(p,'std_A3minus',default_std_A3minus,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_plus',default_mean_tau_plus,positiveNoInfCheck);
            addParameter(p,'std_tau_plus',default_std_tau_plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_x',default_mean_tau_x,positiveNoInfCheck);
            addParameter(p,'std_tau_x',default_std_tau_x,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_minus',default_mean_tau_minus,positiveNoInfCheck);
            addParameter(p,'std_tau_minus',default_std_tau_minus,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_y',default_mean_tau_y,positiveNoInfCheck);
            addParameter(p,'std_tau_y',default_std_tau_y,nonNegativeNoInfCheck);
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
            
            % check if this group has dynamic synapses
            if (p.Results.depressed_synapses || p.Results.facilitating_synapses)
                obj.is_dynamic = true;
            end
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
            %       'mean_A2plus'    - mean doublet LTP factor for
            %                             triplet STDP
            %       'std_A2plus'     - standard deviation of the doublet
            %                             LTP factor.
            %       'mean_A3plus'    - mean triplet LTP factor
            %       'std_A3plus'     - standard deviation of the triplet
            %       'mean_tau_plus'  - mean decay time constant of the r1
            %                             (presynaptic) STDP variable
            %       'std_tau_plus'   - standard deviation of the r1
            %                             decay time constant
            %       'mean_tau_x'     - mean decay time constant of the
            %                             r2 (presynaptic) STDP variable
            %       'std_tau_x'      - standard deviation of the r2
            p = inputParser;
            validNeuronTypes = {'excitatory','inhibitory'};
            default_mean_tau_D = .25; % 250ms (Paul Miller Textbook)
            default_std_tau_D = 0;
            default_mean_tau_F = .25; % 250ms (Paul Miller Textbook)
            default_std_tau_F = 0;
            default_mean_p0 = .4; % (Biro et al., 2005)
            default_std_p0 = 0;
            default_mean_f_fac = .25; % (Paul Miller Textbook)
            default_std_f_fac = 0;
            default_mean_A2plus = 8.8e-11; % Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_A2plus = 0;
            default_mean_A3plus = 5.3e-2; % Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_A3plus = 0;
            default_mean_tau_plus = .0168; % 16.8ms (Bi and Poo, 2001)
            default_std_tau_plus = 0;
            default_mean_tau_x = .714; % 714ms Visual Cortex nearest spike full model (Pfister & Gerstner, 2006)
            default_std_tau_x = 0;
            checkNeuronType = @(x) any(validatestring(neuronType,validNeuronTypes));
            positiveNoInfCheck = @(x) x > 0 && ~isinf(x);
            nonNegativeNoInfCheck = @(x) x>= 0 && ~isinf(x);
            addRequired(p,'name',@ischar);
            addRequired(p,'N',positiveNoInfCheck);
            addRequired(p,'neuronType',checkNeuronType);
            addRequired(p,'firingRate',nonNegativeNoInfCheck);
            addParameter(p,'depressed_synapses',false,@islogical);
            addParameter(p,'facilitating_synapses',false,@islogical);
            addParameter(p,'mean_tau_D',default_mean_tau_D,positiveNoInfCheck);
            addParameter(p,'std_tau_D',default_std_tau_D,positiveNoInfCheck);
            addParameter(p,'mean_tau_F',default_mean_tau_F,positiveNoInfCheck);
            addParameter(p,'std_tau_F',default_std_tau_F,positiveNoInfCheck);
            addParameter(p,'mean_p0',default_mean_p0,nonNegativeNoInfCheck);
            addParameter(p,'std_p0',default_std_p0,positiveNoInfCheck);
            addParameter(p,'mean_f_fac',default_mean_f_fac,nonNegativeNoInfCheck);
            addParameter(p,'std_f_fac',default_std_f_fac,positiveNoInfCheck);
            addParameter(p,'mean_A2plus',default_mean_A2plus,nonNegativeNoInfCheck);
            addParameter(p,'std_A2plus',default_std_A2plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_A3plus',default_mean_A3plus,nonNegativeNoInfCheck);
            addParameter(p,'std_A3plus',default_std_A3plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_plus',default_mean_tau_plus,positiveNoInfCheck);
            addParameter(p,'std_tau_plus',default_std_tau_plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_x',default_mean_tau_x,positiveNoInfCheck);
            addParameter(p,'std_tau_x',default_std_tau_x,nonNegativeNoInfCheck)
            parse(p,name,N,neuronType,firingRate,varargin{:})
            
            obj.nSpikeGenerators = obj.nSpikeGenerators + 1;
            obj.nSpikeGeneratorNeurons = obj.nSpikeGeneratorNeurons + N;
            info.name = p.Results.name;
            info.id = -obj.nSpikeGenerators;
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
            
            % check if this group has dynamic synapses
            if (p.Results.depressed_synapses || p.Results.facilitating_synapses)
                obj.is_dynamic = true;
            end
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
                if (FRNN.checkConnInputs(connType,connParams))
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
                if (EVLIFnetwork.checkConnInputs(connType,connParams))
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
            % call net.print() for a short summary of net.print(true) for a
            % lengthy summary
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
                            else
                                fprintf('    %1$s : %2$f\n',fnames{k},obj.groupInfo(i).connectionParams{j}.(fnames{k}))
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