function [netParams] = make_netParams(varnames,varvalues,structureParams)

for i=1:length(varnames)
    netParams.(varnames{i}) = varvalues(i);
end

netParams.nGroups  = length(structureParams.groupInfo);
start_ind = 1;
for i=1:length(structureParams.groupInfo)
    N = structureParams.groupInfo(i).nNeurons;
    netParams.groupInfo(i).start_ind = start_ind;
    netParams.groupInfo(i).end_ind = start_ind + N - 1;
    start_ind = start_ind+N;
    netParams.groupInfo(i).isExcitatory = structureParams.groupInfo(i).isExcitatory;
    netParams.groupInfo(i).isInhibitory = structureParams.groupInfo(i).isInhibitory;
    netParams.groupInfo(i).targets = structureParams.groupInfo(i).targets;
    netParams.groupInfo(i).target_conn_probs = structureParams.groupInfo(i).target_conn_probs;
    netParams.groupInfo(i).conn_types = structureParams.groupInfo(i).conn_types;
end
end

