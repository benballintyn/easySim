classdef gradientConnector < connectionType
    
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
        function obj = gradientConnector(preGroup,postGroup,connProbFunction,weightFunction,groupInfo)
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
            obj.xcoordsPre = groupInfo(preGroup).xcoords;
            obj.xcoordsPost = groupInfo(postGroup).xcoords;
            obj.ycoordsPre = groupInfo(preGroup).ycoords;
            obj.ycoordsPost = groupInfo(postGroup).ycoords;
            obj.xmax = groupInfo(preGroup).coordinateFrame.xmax;
            obj.ymax = groupInfo(pregroup).coordinateFrame.ymax;
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