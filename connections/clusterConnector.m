classdef clusterConnector < connectionType
    
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
            % This class derives from the class connectionType. It is used
            % to connect a group of neurons to itself. In a clustered
            % connection, each neuron in the group is assigned a cluster ID
            % (based on the number of clusters). Then, neurons in the same
            % cluster are more likely to connect to one another than to
            % neurons in another cluster.
            %
            % clusterConnector(groupID,nClusters,intraConnProb,interConnProb,intraWeightDist,interWeightDist,groupInfo)
            %   groupID         - ID # of group to connect with itself
            %   nClusters       - # of clusters to split the group into
            %   intraConnProb   - connection probability within a cluster
            %   interConnProb   - connection probability between clusters
            %   intraWeightDist - synaptic weight distribution to draw from
            %                     within a cluster
            %   interWeightDist - synaptic weight distribution to draw from
            %                     between clusters
            %   groupInfo       - struct array from network object
            %                     containing information about all groups
            
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
            % call this function to generate the synaptic connectivity
            % matrix defined by this connection. 
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