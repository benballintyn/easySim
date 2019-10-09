classdef gradientConnector < handle
    
    properties
        preGroup
        postGroup
        nPre
        nPost
        connProbFunction
        weightFunction
        xcoordsPre
        xcoordsPost
        ycoordsPre
        ycoordsPost
        xmax
        ymax
    end
    
    methods
        function obj = gaussianConnector(preGroup,postGroup,connProbFunction,weightFunction,netParams)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.preGroup = preGroup;
            obj.postGroup = postGroup;
            obj.nPre = netParams.groupInfo(preGroup).end_ind - netParams.groupInfo(preGroup).start_ind + 1;
            obj.nPost = netParams.groupInfo(postGroup).end_ind - netParams.groupInfo(postGroup).start_ind + 1;
            obj.connProbFunction = connProbFunction;
            obj.weightFunction = weightFunction;
            obj.xcoordsPre = netParams.groupInfo(preGroup).xcoords;
            obj.xcoordsPost = netParams.groupInfo(postGroup).xcoords;
            obj.ycoordsPre = netParams.groupInfo(preGroup).ycoords;
            obj.ycoordsPost = netPrams.groupInfo(postGroup).ycoords;
            obj.xmax = netParams.groupInfo(preGroup).coordinateFrame.xmax;
            obj.ymax = netParams.groupInfo(pregroup).coordinateFrame.ymax;
        end
        
        function dGsynMat = genConn(obj)
            dGsynMat = zeros(obj.nPost,obj.nPre);
            for i=1:obj.nPre
                for j=1:obj.nPost
                    if (obj.preGroup == obj.postGroup && i == j)
                        continue;
                    else
                        connProb = obj.connProbFunction(obj.xcoordsPost(j));
                        if (rand < connProb)
                            dGsynMat(j,i) = obj.weightFunction(obj.xcoordsPost(j));
                        end
                    end
                end
            end
        end
    end
end