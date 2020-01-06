classdef ELIFnetwork < handle
    
    properties
       groupInfo
       nGroups
    end
    
    methods
        function obj = ELIFnetwork()
            obj.nGroups = 0;
            obj.groupInfo = struct('isExcitatory',{},'isInhibitory',{},...
                                   'targets',{},'connectionTypes',{},...
                                   'connectionParams',{},'start_ind',{},...
                                   'end_ind',{},'coordinateFrame',{},...
                                   'xcoords',{},'ycoords',{});
        end
        
        function addGroup(obj,neuronType,coordinateFrame)
            obj.nGroups = obj.nGroups + 1;
            if (strcmp(neuronType,'excitatory'))
                obj.groupInfo(obj.nGroups).isExcitatory = true;
                obj.groupInfo(obj.nGroups).isInhibitory = false;
            elseif (strcmp(neuronType,'inhibitory'))
                obj.groupInfo(obj.nGroups).isInhibitory = true;
                obj.groupInfo(obj.nGroups).isExcitatory = false;
            else
                error('neuronType is not valid. Must be "excitatory" or "inhibitory"')
            end
            obj.groupInfo(obj.nGroups).coordinateFrame.ID = coordinateFrame;
            
        end
        
    end
end