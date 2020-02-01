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
    end
    
    methods
        function obj = AEVLIFnetwork()
            % This class holds all of the necessary information about an
            % AEVLIFnetwork. The A in AELIF stands for 'adaptation', referring
            % to the fact that neurons in this network type have an adaptation
            % current. The E in AEVLIF stands for 'exponential' and
            % refers to the spike generation mechansim used. the V in AEVLIF
            % refers to the fact that the refractory period is modeled by
            % instantaneously increasing the spike threshold after a spike
            % (making it more difficult to spike again). This increased
            % threshold decays back to baseline with some time constant.
            % LIF stands for the standard leaky-integrate-and-fire model
            % neuron.
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
                'std_El',{},'mean_dth',{},'std_dth',{},'mean_A2plus',{},'std_A2plus',{},...
                'mean_A3plus',{},'std_A3plus',{},'mean_A2minus',{},'std_A2minus',{},...
                'mean_A3minus',{},'std_A3minus',{},'mean_tau_plus',{},'std_tau_plus',{},...
                'mean_tau_x',{},'std_tau_x',{},'mean_tau_minus',{},'std_tau_minus',{},...
                'mean_tau_y',{},'std_tau_y',{},'xcoords',{},'ycoords',{},'record',{});
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
            %       8. 'mean_Vth_max'   - mean maximum threshold value. A
            %                             neuron's spike threshold will be 
            %                             set to this value following a spike
            %       9. 'std_Vth_max'    - standard deviation of the maximum
            %                             threshold value
            %      10. 'mean_tau_ref'   - mean refractory period time
            %                             constant
            %      11. 'std_tau_ref'    - standard deviation of the
            %                             refractory period time constant
            %      12. 'mean_tau_sra'   - mean decay time constant for the
            %                             adaptation current
            %      13. 'std_tau_sra'    - standard deviation of the
            %                             adaptation current time constant
            %      14. 'mean_a'         - mean conductance of the
            %                             adaptation current
            %      15. 'std_a'          - standard deviation of the
            %                             conductance of the adaptation
            %                             current
            %      16. 'mean_b'         - mean increment to the adaptation
            %                             current following a spike
            %      17. 'std_b'          - standard deviation of the increment
            %                             to the adaptation current
            %                             following a spike
            %      18. 'mean_VsynE'     - mean excitatory synaptic reversal
            %                             potential
            %      19. 'std_VsynE'      - standard deviation of the
            %                             excitatory synaptic reversal potential
            %      20. 'mean_VsynI'     - mean inhibitory synaptic reversal
            %                             potential
            %      21. 'std_VsynI'      - standard deviation of the
            %                             inhibitory synaptic reversal
            %                             potential
            %      22. 'mean_max_GsynE' - mean maximum total excitatory
            %                             synaptic conductance
            %      23. 'std_max_GsynE'  - standard deviation of the maximum
            %                             total excitatory synaptic
            %                             conductance
            %      24. 'mean_max_GsynI' - mean maximum total inhibitory
            %                             synaptic conductance
            %      25. 'std_max_GsynI'  - standard deviation of the maximum
            %                             total inhibitory synaptic
            %                             conductance
            %      26. 'mean_tau_synE'  - mean excitatory synaptic time
            %                             constant
            %      27. 'std_tau_synE'   - standard deviation of the
            %                             excitatory synaptic time constants
            %      28. 'mean_tau_synI'  - mean inhibitory synaptic time
            %                             constant
            %      29. 'std_tau_synI'   - standard deviation of the
            %                             inhibitory synaptic time constants
            %      30. 'mean_Cm'        - mean membrane capacitance
            %      31. 'std_Cm'         - standard devation of membrane
            %                             capacitances
            %      32. 'mean_Gl'        - mean leak conductance
            %      33. 'std_Gl'         - standard deviation of leak
            %                             conductances
            %      34. 'mean_El'        - mean leak reversal potential
            %      35. 'std_El'         - standard deviation in leak
            %                             reversal potentials
            %      36. 'mean_dth'       - mean spike generation voltage range
            %      37. 'std_dth'        - standard deviation of the spike
            %                             generation voltage range
            %      38. 'record'         - can be true or false. Indicates
            %                             whether to record this groups spike
            %                             times to a file or not
            %      39. 'xmin'           - minimum x-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      40. 'xmax'           - maximum x-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      41. 'ymin'           - minimum y-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            %      42. 'ymax'           - maximum y-value of the coordinate
            %                             frame this group lies in. NOTE:
            %                             If using this optional parameter,
            %                             this neuron group must be the
            %                             first in this coordinate frame
            
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
                'std_El','mean_dth','std_dth','mean_A2plus','std_A2plus',...
                'mean_A3plus','std_A3plus','mean_A2minus','std_A2minus',...
                'mean_A3minus','std_A3minus','mean_tau_plus','std_tau_plus',...
                'mean_tau_x','std_tau_x','mean_tau_minus','std_tau_minus',...
                'mean_tau_y','std_tau_y','xcoords','ycoords','record'};
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
            default_mean_VsynI = -.07; % -80mV
            default_std_VsynI = 0;
            default_mean_max_GsynE = 10e-6; % 10uS
            default_std_max_GsynE = 0;
            default_mean_max_GsynI = 15e-6; % 15uS
            default_std_max_GsynI = 0;
            default_mean_tau_synE = 20e-3; % 10ms
            default_std_tau_synE = 0;
            default_mean_tau_synI = 10e-3; % 1ms
            default_std_tau_synI = 0;
            default_mean_Cm = 10e-9; % 10nF/mm^2 (Dayan & Abbott) .9e-6; %.uF/cm^2 (Gentet et al., 2000)
            default_std_Cm = 0;
            default_mean_Gl = 1e-6; % 1uS/mm^2 (Dayan & Abbott)
            default_std_Gl = 0;
            default_mean_El = -.07; % -70mV
            default_std_El = 0;
            default_mean_dth = .002; %2mV
            default_std_dth = 0;
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
            % addSpikeGenerator(name,N,neuronType,firingRate)
            %   name       - string name for this spikeGenerator group
            %
            %   N          - # of neurons in this spikeGenerator group
            %
            %   neuronType - 'excitatory' or 'inhibitory'
            %
            %   firingRate - firing rate of this group in Hz
            p = inputParser;
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
            addParameter(p,'mean_A2plus',default_mean_A2plus,nonNegativeNoInfCheck);
            addParameter(p,'std_A2plus',default_std_A2plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_A3plus',default_mean_A3plus,nonNegativeNoInfCheck);
            addParameter(p,'std_A3plus',default_std_A3plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_plus',default_mean_tau_plus,positiveNoInfCheck);
            addParameter(p,'std_tau_plus',default_std_tau_plus,nonNegativeNoInfCheck);
            addParameter(p,'mean_tau_x',default_mean_tau_x,positiveNoInfCheck);
            addParameter(p,'std_tau_x',default_std_tau_x,nonNegativeNoInfCheck)
            
            obj.nSpikeGenerators = obj.nSpikeGenerators + 1;
            obj.nSpikeGeneratorNeurons = obj.nSpikeGeneratorNeurons + N;
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