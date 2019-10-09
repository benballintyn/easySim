classdef clusterConnector < handle
    
    properties
        preGroup
        postGroup
        intraConnProb
        interConnProb
        clusterIDs
        intraWeightDist
        interWeightDist
        nPre
        nPost
    end
    
    methods
        function obj = clusterConnector(preGroup,postGroup,clusterIDs,intraConnProb,interConnProb,intraWeightDist,interWeightDist,netParams)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            if (preGroup ~= postGroup)
                error('clusterConnector only supports connecting a group to itself')
            end
            obj.preGroup = preGroup;
            obj.postGroup = postGroup;
            obj.clusterIDs = clusterIDs;
            obj.intraConnProb = intraConnProb;
            obj.interConnProb = interConnProb;
            obj.intraWeightDist = intraWeightDist;
            obj.interWeightDist = interWeightDist;
            obj.nPre = netParams.groupInfo(preGroup).end_ind - netParams.groupInfo(preGroup).start_ind + 1;
            obj.nPost = netParams.groupInfo(postGroup).end_ind - netParams.groupInfo(postGroup).start_ind + 1;
        end
        
        function dGsynMat = genConn(obj)
            c = rand(obj.nPost,obj.nPre);
            dGsynMat = zeros(obj.nPost,obj.nPre);
            for i=1:obj.nPre
                for j=1:obj.nPost
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