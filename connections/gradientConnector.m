classdef gradientConnector < connectionType
    
    properties
        %preGroup % defined in superclass
        %postGroup % defined in superclass
        nPre
        nPost
        connProbFunction
        weightFunction
        xcoordsPost
    end
    
    methods
        function obj = gradientConnector(preGroup,postGroup,connProbFunction,weightFunction,groupInfo,spikeGeneratorInfo)
            % This class derives from the class connectionType. It is used
            % to connect a group of neurons to another (or the same) group.
            % With a gradient connection, each presynaptic neuron is
            % connected to a postsynaptic neuron with a probability that is
            % dependent on the postsynaptic neurons x-coordinate.
            %
            % clusterConnector(groupID,nClusters,intraConnProb,interConnProb,intraWeightDist,interWeightDist,groupInfo)
            %   preGroup            - ID # of presynaptic group
            %   postGroup           - ID # of postsynaptic group
            %   connProbFunction    - function that outputs the probability
            %                         of connection dependent on the
            %                         postsynaptic neuron's x-coordinate
            %   weightFunction      - function that outputs a synaptic
            %                         weight dependent on the postsynaptic
            %                         neuron's x-coordinate
            %   groupInfo           - struct array from network object
            %                         containing information about all groups
            %   spikeGeneratorInfo  - struct array from network object
            %                         containing information about all
            %                         spikeGenerator groups
            obj.preGroup = preGroup;
            obj.postGroup = postGroup;
            if (preGroup < 0)
                obj.nPre = spikeGeneratorInfo(-preGroup).N;
            else
                obj.nPre = groupInfo(preGroup).end_ind - groupInfo(preGroup).start_ind + 1;
            end
            obj.nPost = groupInfo(postGroup).end_ind - groupInfo(postGroup).start_ind + 1;
            obj.connProbFunction = connProbFunction;
            obj.weightFunction = weightFunction;
            obj.xcoordsPost = groupInfo(postGroup).xcoords;
        end
        
        function dGsynMat = genConn(obj)
            % call this function to generate the synaptic connectivity
            % matrix defined by this connection. 
            dGsynMat = zeros(obj.nPost,obj.nPre);
            for i=1:obj.nPre
                for j=1:obj.nPost
                    if (obj.preGroup == obj.postGroup && i == j)
                        continue;
                    else
                        connProb = obj.connProbFunction(obj.xcoordsPost(j));
                        if (rand < connProb)
                            dGsynMat(j,i) = max(0,obj.weightFunction(obj.xcoordsPost(j)));
                        end
                    end
                end
            end
        end
    end
end
