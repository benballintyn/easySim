function [netParams] = make_netParams(varnames,varvalues,structureParams)

for i=1:length(varnames)
    netParams.(varnames{i}) = varvalues(i);
end

netParams.nGroups  = length(structureParams.groupInfo);
start_ind = 1;

end

