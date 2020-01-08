classdef ELIFnetwork < handle
    
    properties
       groupInfo
       nGroups
       nNeurons
       coordinateFrames
    end
    
    methods
        function obj = ELIFnetwork()
            obj.nGroups = 0;
            obj.groupInfo = struct('id',{},'isExcitatory',{},'isInhibitory',{},...
                                   'targets',{},'connections',{},...
                                   'connectionParams',{},'start_ind',{},...
                                   'end_ind',{},'coordinateFrame',{},...
                                   'xcoords',{},'ycoords',{});
            obj.nNeurons = 0;
            obj.coordinateFrames = struct('ID',{},'xmin',{},'xmax',{},'ymin',{},'ymax',{});
        end
        
        function addGroup(obj,name,N,neuronType,coordinateFrame,varargin)
            % Create inputParser and assign default values and checks
            p = inputParser;
            validNeuronTypes = {'excitatory','inhibitory'};
            % default values
            defaultXmin = 0;
            defaultXmax = 1;
            defaultYmin = 0;
            defaultYmax = 1;
            default_mean_V0 = -.07; % -70mV
            default_std_V0 = 0;
            default_mean_Vth = -.05; % -50mV
            default_std_Vth = 0;
            default_mean_dGref = ;
            checkNeuronType = @(x) any(validatestring(neuronType,validNeuronTypes));
            validNumCheck = @(x) isnumeric(x) && ~isinf(x) && ~isnan(x);
            positiveNoInfCheck = @(x) x > 0 && ~isinf(x);
            addRequired(p,'name',@ischar)
            addRequired(p,'N',positiveNoInfCheck);
            addRequired(p,'neuronType',checkNeuronType);
            addRequired(p,'coordinateFrame',positiveNoInfCheck);
            addParameter(p,'xmin',defaultXmin,validNumCheck);
            addParameter(p,'xmax',defaultXmax,validNumCheck);
            addParameter(p,'ymin',defaultYmin,validNumCheck);
            addParameter(p,'ymax',defaultYmax,validNumCheck);
            addParameter(p,'mean_V0',default_mean_V0,validNumCheck);
            addParameter(p,'std_V0',default_std_V0,validNumCheck);
            addParameter(p,'mean_Vth',default_mean_Vth,validNumCheck);
            addParameter(p,'std_Vth',default_std_Vth,validNumCheck);
            addParameter(p,'mean_dGref',default_mean_dGref,validNumCheck);
            addParameter(p,'std_dGref',default_std_dGref,validNumCheck);
            addParameter(p,'mean_tau_ref',default_mean_tau_ref,validNumCheck);
            addParameter(p,'std_tau_ref',default_std_tau_ref,validNumCheck);
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
            parse(p,name,N,neuronType,coordinateFrame,varargin{:})
            
            % Housekeeping to keep track of # of groups and neurons
            obj.nGroups = obj.nGroups + 1;
            obj.groupInfo(obj.nGroups).name = name;
            obj.groupInfo(obj.nGroups).id = obj.nGroups;
            obj.groupInfo(obj.nGroups).start_ind = obj.nNeurons+1;
            obj.groupInfo(obj.nGroups).end_ind = obj.nNeurons + p.Results.N;
            obj.nNeurons = obj.nNeurons + p.Results.N;
            
            % Determine whether the group is excitatory or inhibitory
            obj.groupInfo(obj.nGroups).isExcitatory = strcmp(p.Results.neuronType,'excitatory');
            obj.groupInfo(obj.nGroups).isInhibitory = strcmp(p.Results.neuronType,'inhibitory');
            
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
            obj.groupInfo(obj.nGroups).coordinateFrame = cf;
            
            % assign x and y coordinates to each new neuron
            obj.groupInfo(obj.nGroups).xcoords = rand(1,p.Results.N)*(p.Results.xmax - p.Results.xmin) + p.Results.xmin;
            obj.groupInfo(obj.nGroups).ycoords = rand(1,p.Results.N)*(p.Results.ymax - p.Results.ymin) + p.Results.ymin;
        end
        
        function connect(obj,src_id,tgt_id,connType,connParams)
            if (ELIFnetwork.checkConnInputs(connType,connParams))
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