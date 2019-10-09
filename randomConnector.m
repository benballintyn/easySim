classdef randomConnector < handle
    
    properties
        preGroup
        postGroup
        connProb
        weightDist
        nPre
        nPost
    end
    
    methods
        function obj = randomConnector(preGroup,postGroup,connProb,weightDist,netParams)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.preGroup = preGroup;
            obj.postGroup = postGroup;
            obj.connProb = connProb;
            obj.weightDist = weightDist;
            obj.nPre = netParams.groupInfo(preGroup).end_ind - netParams.groupInfo(preGroup).start_ind + 1;
            obj.nPost = netParams.groupInfo(postGroup).end_ind - netParams.groupInfo(postGroup).start_ind + 1;
        end
        
        function dGsynMat = genConn(obj)
            c = rand(obj.nPost,obj.nPre);
            if (obj.preGroup == obj.postGroup)
                c(logical(eye(obj.nPre))) = 0;
            end
            c(c > obj.connProb) = 0; c(c ~= 0) = 1;
            dGsynMat = obj.weightDist.genSample([obj.nPost obj.nPre]);
            dGsynMat = c.*dGsynMat;
        end
    end
end

