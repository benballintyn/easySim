classdef randomConnector < connectionType
    
    properties
        %preGroup %defined by superclass
        %postGroup %defined by superclass
        connProb
        weightDist
        nPre
        nPost
    end
    
    methods
        function obj = randomConnector(preGroup,postGroup,connProb,weightDist,groupInfo,spikeGeneratorInfo)
            % This class derives from the class connectionType. It is used
            % to randomly connect a group of neurons to another (or the same)
            % group. With a random connection, each presynaptic neuron is randomly
            % connected to a postsynaptic neuron with a certain fixed
            % probability.
            %
            % clusterConnector(groupID,nClusters,intraConnProb,interConnProb,intraWeightDist,interWeightDist,groupInfo)
            %   preGroup            - ID # of presynaptic group
            %   postGroup           - ID # of postsynaptic group
            %   connProb            - probability of connection between each
            %                         presynaptic-postsynaptic neuron pair
            %   weightDist          - synaptic weight distribution to draw from
            %   groupInfo           - struct array from network object
            %                         containing information about all groups
            %   spikeGeneratorInfo  - struct array from network object
            %                         containing information about all
            %                         spikeGenerator groups
            obj.preGroup = preGroup;
            obj.postGroup = postGroup;
            obj.connProb = connProb;
            obj.weightDist = weightDist;
            if (preGroup < 0)
                obj.nPre = spikeGeneratorInfo(-preGroup).N;
            else
                obj.nPre = groupInfo(preGroup).end_ind - groupInfo(preGroup).start_ind + 1;
            end
            obj.nPost = groupInfo(postGroup).end_ind - groupInfo(postGroup).start_ind + 1;
        end
        
        function dGsynMat = genConn(obj)
            % call this function to generate the synaptic connectivity
            % matrix defined by this connection
            c = rand(obj.nPost,obj.nPre);
            if (obj.preGroup == obj.postGroup)
                c(logical(eye(obj.nPre))) = 0;
            end
            c(c > obj.connProb) = 0; c(c ~= 0) = 1;
            dGsynMat = obj.weightDist.genSample([obj.nPost obj.nPre]);
            dGsynMat = max(0,c.*dGsynMat);
        end
    end
end

