classdef weightDistribution < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xrange
        px
    end
    
    methods
        function obj = weightDistribution(xrange,px)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            if (length(xrange) ~= length(px))
                error('xrange and px are of different lengths')
            else
                obj.xrange = xrange;
                obj.px = px;
            end
        end
        
        function randomSample = genSample(obj,sampleSize)
            randomSample = randpdf(obj.px,obj.xrange,sampleSize);
        end
    end
end

