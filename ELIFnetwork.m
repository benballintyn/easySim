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
            obj.coordinateFrames = [];
        end
        
        function addGroup(obj,N,neuronType,coordinateFrame)
            % Simple validity check of neuronType
            if (~strcmp(neuronType,'excitatory') && ~strcmp(neuronType,'inhibitory'))
                error('Invalid neuronType. Must be "excitatory" or "inhibitory"')
            end
                
            % Housekeeping to keep track of # of groups and neurons
            obj.nGroups = obj.nGroups + 1;
            obj.groupInfo(obj.nGroups).id = obj.nGroups;
            obj.start_ind = obj.nNeurons+1;
            obj.end_ind = obj.nNeurons + N;
            
            % Determine whether the group is excitatory or inhibitory
            obj.groupInfo(obj.nGroups).isExcitatory = strcmp(neuronType,'excitatory');
            obj.groupInfo(obj.nGroups).isInhibitory = strcmp(neuronType,'inhibitory');
            
            % Set coordinate frame of new group
            obj.groupInfo(obj.nGroups).coordinateFrame.ID = coordinateFrame;
        end
        
        function connect(obj,src_id,tgt_id,connType,connParams)
            if (checkConnInputs(connType,connParams))
                if (strcmp(connType,'random'))
                    conn=randomConnector(src_id,tgt_id,connParams.conn_prob,connParams.weight_distribution,obj.groupInfo);
                    obj.groupInfo(src_id).targets = [obj.groupInfo(src_id).targets tgt_id];
                    obj.groupInfo(src_id).connections = [obj.groupInfo(src_id).connections conn];
                    obj.groupInfo(src_id).connectionParams = [obj.groupInfo(src_id).connectionParams connParams];
                elseif (strcmp(connType,'clustered'))
                    if (src_id ~= tgt_id)
                        error('Source and target must be the same for clustered connections')
                    end
                    conn=clusterConnector(src_id,connParams.nClusters,connParams.intraConnProb,connParams.interConnProb,connParams.intraWeightDist,connParams.interWeightDist,obj.groupInfo);
                    obj.groupInfo(src_id).targets = [obj.groupInfo(src_id).targets tgt_id];
                    obj.groupInfo(src_id).connections = [obj.groupInfo(src_id).connections conn];
                    obj.groupInfo(src_id).connectionParams = [obj.groupInfo(src_id).connectionParams connParams];
                end
            end
        end
    end
    
    methods (Static)
        function all_ok=checkConnInputs(connType,connParams)
            if (strcmp(connType,'random'))
                reqFieldNames = {'conn_prob','weight_distribution'};
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
end