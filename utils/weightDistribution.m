classdef weightDistribution < handle
    
    properties
        xrange
        px
    end
    
    methods
        function obj = weightDistribution(xrange,px)
            % This is a simple class to represent a probability
            % distribution over some range of values. This can be used to
            % represent an arbitrary distribution since the probabilities
            % of each value are inputs.
            %
            % weightDistribution(xrange,px)
            %   xrange - vector of values representing the domain of the
            %            probability distribution
            %   px     - must be the same length as xrange. Each value
            %            px(i) gives the associated probability density to
            %            the value xrange(i).
            %
            %   METHODS:
            %       genSample(sampleSize)
            %           - samples from the distribution and returns the
            %             samples in a matrix of size sampleSize.
            
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

