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
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
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
            dGsynMat = zeros(obj.nPost,obj.nPre);
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