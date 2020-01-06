classdef clusterConnector < handle
    
    properties
        groupID
        intraConnProb
        interConnProb
        nClusters
        clusterIDs
        intraWeightDist
        interWeightDist
        nNeurons
    end
    
    methods
        function obj = clusterConnector(groupID,nClusters,intraConnProb,interConnProb,intraWeightDist,interWeightDist,groupInfo)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            if (preGroup ~= postGroup)
                error('clusterConnector only supports connecting a group to itself')
            end
            obj.groupID = groupID;
            obj.intraConnProb = intraConnProb;
            obj.interConnProb = interConnProb;
            obj.intraWeightDist = intraWeightDist;
            obj.interWeightDist = interWeightDist;
            obj.nNeurons = groupInfo(obj.groupID).end_ind - groupInfo(obj.groupID).start_ind + 1;
            neurons_per_cluster = floor(obj.nNeurons / nClusters);
            ids=1:nClusters-1;
            remainder = mod(obj.nNeurons,nClusters);
            obj.clusterIDs = repelem(ids,neurons_per_cluster);
            obj.clusterIDs = [obj.clusterIDs repelem(nClusters,remainder)];
        end
        
        function dGsynMat = genConn(obj)
            c = rand(obj.nNeurons,obj.nNeurons);
            dGsynMat = zeros(obj.nNeurons,obj.nNeurons);
            for i=1:obj.nNeurons
                for j=1:obj.nNeurons
                    if (i == j)
                        dGsynMat(j,i) = 0;
                        continue
                    end
                    if (obj.clusterIDs(i) == obj.clusterIDs(j))
                        if (c(j,i) < obj.intraConnProb)
                            dGsynMat(j,i) = obj.intraWeightDist.genSample([1 1]);
                        else
                            dGsynMat(j,i) = 0;
                        end
                    else
                        if (c(j,i) < obj.interConnProb)
                            dGsynMat(j,i) = obj.interWeightDist.genSample([1 1]);
                        else
                            dGsynMat(j,i) = 0;
                        end
                    end
                end
            end
        end
    end
end