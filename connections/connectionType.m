classdef connectionType < matlab.mixin.Heterogeneous & handle
    
    properties (SetAccess = protected, GetAccess = protected)
        preGroup
        postGroup
    end
    
    methods (Abstract)
        genConn(obj)
    end
end

