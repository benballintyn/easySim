classdef gaussianConnector < handle
    
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
        useWrap
        xmax
        ymax
    end
    
    methods
        function obj = gaussianConnector(preGroup,postGroup,connProbFunction,weightFunction,useWrap,netParams)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            if (netParams.groupInfo(preGroup).coordinateFrame.ID ~= netParams.groupInfo(postGroup).coordinateFrame.ID)
                error('gaussianConnector only supports connecting groups in the same coordinate frame')
            end
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
            obj.useWrap = useWrap;
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
                        dx = abs(obj.xcoordsPre(i) - obj.xcoordsPost(j));
                        dy = abs(obj.ycoordsPre(i) - obj.ycoordsPost(j));
                        if (obj.useWrap)
                            if (dx > .5*obj.xmax)
                                dx = obj.xmax - dx;
                            end
                            if (dy > .5*obj.ymax)
                                dy = obj.ymax - dy;
                            end
                        end
                        distance = sqrt(dx.^2 + dy.^2);
                        connProb = obj.connProbFunction(distance);
                        if (rand < connProb)
                            dGsynMat(j,i) = obj.weightFunction(distance);
                        end
                    end
                end
            end
        end
    end
end