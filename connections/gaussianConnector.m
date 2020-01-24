classdef gaussianConnector < connectionType
    
    properties
        %preGroup % defined in superclass
        %postGroup % defined in superclass
        nPre
        nPost
        connProbFunction
        weightFunction
        xcoordsPre
        xcoordsPost
        ycoordsPre
        ycoordsPost
        useWrap
        xmax
        xmin
        ymax
        ymin
    end
    
    methods
        function obj = gaussianConnector(preGroup,postGroup,connProbFunction,weightFunction,useWrap,groupInfo)
            % This class derives from the class connectionType. It is used
            % to connect a group of neurons to another (or the same) group.
            % With a gaussian connection, each presynaptic neuron is
            % connected to a postsynaptic neuron with a probability that is
            % dependent on their distance in 2D space. Neuron groups
            % connected by a gaussianConnector MUST be in the same
            % coordinate frame.
            %
            % clusterConnector(groupID,nClusters,intraConnProb,interConnProb,intraWeightDist,interWeightDist,groupInfo)
            %   preGroup            - ID # of presynaptic group
            %   postGroup           - ID # of postsynaptic group
            %   connProbFunction    - function that outputs the probability
            %                         of connection dependent on distance
            %   weightFunction      - function that outputs a synaptic
            %                         weight dependent on distance
            %   useWrap             - true or false. If true, distances
            %                         between neurons will be computed as
            %                         if the 2D plane the neurons lie in
            %                         wraps around on itself.
            %   groupInfo           - struct array from network object
            %                         containing information about all groups
            if (preGroup < 0)
                error('gaussianConnector does not support spike generator groups')
            end
            if (groupInfo(preGroup).coordinateFrame.ID ~= groupInfo(postGroup).coordinateFrame.ID)
                error('gaussianConnector only supports connecting groups in the same coordinate frame')
            end
            obj.preGroup = preGroup;
            obj.postGroup = postGroup;
            obj.nPre = groupInfo(preGroup).end_ind - groupInfo(preGroup).start_ind + 1;
            obj.nPost = groupInfo(postGroup).end_ind - groupInfo(postGroup).start_ind + 1;
            obj.connProbFunction = connProbFunction;
            obj.weightFunction = weightFunction;
            obj.xcoordsPre = groupInfo(preGroup).xcoords;
            obj.xcoordsPost = groupInfo(postGroup).xcoords;
            obj.ycoordsPre = groupInfo(preGroup).ycoords;
            obj.ycoordsPost = groupInfo(postGroup).ycoords;
            obj.useWrap = useWrap;
            obj.xmax = groupInfo(preGroup).coordinateFrame.xmax;
            obj.xmin = groupInfo(preGroup).coordinateFrame.xmin;
            obj.ymax = groupInfo(preGroup).coordinateFrame.ymax;
            obj.ymin = groupInfo(preGroup).coordinateFrame.ymin;
        end
        
        function dGsynMat = genConn(obj)
            % call this function to generate the synaptic connectivity
            % matrix defined by this connection. 
            distances = zeros(obj.nPost,obj.nPre);
            xrng = obj.xmax - obj.xmin;
            yrng = obj.ymax - obj.ymin;
            for i=1:obj.nPre
                for j=1:obj.nPost
                    if (obj.preGroup == obj.postGroup && i == j)
                        continue;
                    else
                        dx = abs(obj.xcoordsPre(i) - obj.xcoordsPost(j));
                        dy = abs(obj.ycoordsPre(i) - obj.ycoordsPost(j));
                        if (obj.useWrap)
                            if (dx > .5*xrng)
                                dx = xrng - dx;
                            end
                            if (dy > .5*yrng)
                                dy = yrng - dy;
                            end
                        end
                        distances(j,i) = sqrt(dx.^2 + dy.^2);
                    end
                end
            end
            connProb = obj.connProbFunction(distances);
            dGsynMat = (rand(size(connProb)) < connProb);
            if (obj.preGroup == obj.postGroup)
                dGsynMat(logical(eye(obj.nPre))) = 0;
            end
            dGsynMat = dGsynMat.*obj.weightFunction(distances);
        end
    end
end