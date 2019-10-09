function [V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt] = setupNet(netParams)
% Total # of neurons to be simulated
N = netParams.nOBcells + netParams.nGCcells + netParams.nPCEcells + netParams.nPCIcells;

% Variables that may change with time
V = ones(N,1)*netParams.El;
Gref = zeros(N,1);
GsynE = zeros(N,1);
GsynI = zeros(N,1);
Iapp = zeros(N,1);
Vth = ones(N,1)*netParams.Vth;
VsynE = ones(N,1)*netParams.VsynE;
VsynI = ones(N,1)*netParams.VsynI;

% Variables that will not change with time
dGref = netParams.dGref;
tau_ref = netParams.tau_ref;
dGsyn = zeros(N,N);
tau_synE = netParams.tau_synE;
tau_synI = netParams.tau_synI;
Cm  = netParams.Cm;
Gl = netParams.Gl;
El = netParams.El;
Ek = netParams.Ek;
dth = netParams.dth;
dt = netParams.dt;

for i=1:netParams.nGroups
    preStart = netParams.groupInfo(i).start_ind;
    preEnd = netParams.groupInfo(i).end_ind;
    targets = netParams.groupInfo(i).targets;
    conn_probs = netParams.groupInfo(i).target_conn_probs;
    mean_conn_strengths = netParams.groupInfo(i).mean_conn_strengths;
    std_conn_strengths = netParams.groupInfo(i).std_conn_strenghts;
    connectTypes = netParams.groupInfo(i).connection_types;
    for j=1:length(targets)
        postStart = netParams.groupInfo(targets(j)).start_ind;
        postEnd = netParams.groupInfo(targets(j)).end_ind;
        w = getWeightMatrix(i,j,netParams);
        dGsyn(postStart:postEnd,preStart:preEnd) = w;
    end
end


end