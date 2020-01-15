function [V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells] = setupNet(net,simobj,useGpu)
% Total # of neurons to be simulated
N = net.nNeurons;

if (useGpu)
    % Variables that may change with time
    V = gpuArray(zeros(N,1,'single'));
    Gref = gpuArray(zeros(N,1,'single'));
    GsynE = gpuArray(zeros(N,1,'single'));
    GsynI = gpuArray(zeros(N,1,'single'));
    maxGsynE = gpuArray(zeros(N,1,'single'));
    maxGsynI = gpuArray(zeros(N,1,'single'));
    Iapp = gpuArray(ones(N,1,'single'))*2e-10;
    Vth = gpuArray(zeros(N,1,'single'));
    VsynE = gpuArray(ones(N,1,'single'));
    VsynI = gpuArray(ones(N,1,'single'));

    % Variables that will not change with time
    dGref = gpuArray(zeros(N,1,'single'));
    tau_ref = gpuArray(zeros(N,1,'single'));
    dGsyn = gpuArray(zeros(N,N,'single'));
    tau_synE = gpuArray(zeros(N,1,'single'));
    tau_synI = gpuArray(zeros(N,1,'single'));
    Cm  = gpuArray(zeros(N,1,'single'));
    Gl = gpuArray(zeros(N,1,'single'));
    El = gpuArray(zeros(N,1,'single'));
    Ek = gpuArray(zeros(N,1,'single'));
    dth = gpuArray(zeros(N,1,'single'));
    dt = simobj.dt;
    ecells = gpuArray(zeros(net.nNeurons,1));
    icells = gpuArray(zeros(net.nNeurons,1));
else
    % Variables that may change with time
    V = zeros(N,1,'single');
    Gref = zeros(N,1,'single');
    GsynE = zeros(N,1,'single');
    GsynI = zeros(N,1,'single');
    maxGsynE = zeros(N,1,'single');
    maxGsynI = zeros(N,1,'single');
    Iapp = ones(N,1,'single')*2e-10;
    Vth = zeros(N,1,'single');
    VsynE = ones(N,1,'single');
    VsynI = ones(N,1,'single');

    % Variables that will not change with time
    dGref = zeros(N,1,'single');
    tau_ref = zeros(N,1,'single');
    dGsyn = zeros(N,N,'single');
    tau_synE = zeros(N,1,'single');
    tau_synI = zeros(N,1,'single');
    Cm  = zeros(N,1,'single');
    Gl = zeros(N,1,'single');
    El = zeros(N,1,'single');
    Ek = zeros(N,1,'single');
    dth = zeros(N,1,'single');
    dt = simobj.dt;
    ecells = zeros(N,1,'single');
    icells = zeros(N,1,'single');
end
for i=1:net.nGroups
    s = net.groupInfo(i).start_ind;
    e = net.groupInfo(i).end_ind;
    if (net.groupInfo(i).isExcitatory)
        ecells(s:e) = 1;
    elseif (net.groupInfo(i).isInhibitory)
        icells(s:e) = 1;
    else
        error('Group is neither excitatory or inhibitory?')
    end
end
ecells = logical(ecells);
icells = logical(icells);

for i=1:net.nGroups
    preStart = net.groupInfo(i).start_ind;
    preEnd = net.groupInfo(i).end_ind;
    groupN = preEnd - preStart + 1;
    % initialize all nonzero variables for current group (V(0), Vth, VsynE, VsynI, dGref,
    % tau_ref, tau_synE, tau_synI, Cm, Gl, El, Ek, dth)
    V(preStart:preEnd) = normrnd(net.groupInfo(i).mean_V0,net.groupInfo(i).std_V0,groupN,1);
    Vth(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Vth,net.groupInfo(i).std_Vth,groupN,1);
    VsynE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_VsynE,net.groupInfo(i).std_VsynE,groupN,1);
    VsynI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_VsynI,net.groupInfo(i).std_VsynI,groupN,1);
    maxGsynE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_max_GsynE,net.groupInfo(i).std_max_GsynE,groupN,1);
    maxGsynI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_max_GsynI,net.groupInfo(i).std_max_GsynI,groupN,1);
    dGref(preStart:preEnd) = normrnd(net.groupInfo(i).mean_dGref,net.groupInfo(i).std_dGref,groupN,1);
    tau_ref(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_ref,net.groupInfo(i).std_tau_ref,groupN,1);
    tau_synE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_synE,net.groupInfo(i).std_tau_synE,groupN,1);
    tau_synI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_synI,net.groupInfo(i).std_tau_synI,groupN,1);
    Cm(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Cm,net.groupInfo(i).std_Cm,groupN,1);
    Gl(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Gl,net.groupInfo(i).std_Gl,groupN,1);
    El(preStart:preEnd) = normrnd(net.groupInfo(i).mean_El,net.groupInfo(i).std_El,groupN,1);
    Ek(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Ek,net.groupInfo(i).std_Ek,groupN,1);
    dth(preStart:preEnd) = normrnd(net.groupInfo(i).mean_dth,net.groupInfo(i).std_dth,groupN,1);
    
    targets = net.groupInfo(i).targets;
    for j=1:length(targets)
        postStart = net.groupInfo(targets(j)).start_ind;
        postEnd = net.groupInfo(targets(j)).end_ind;
        curdGsyn = net.groupInfo(i).connections(j).genConn();
        dGsyn(postStart:postEnd,preStart:preEnd) = curdGsyn;
    end
end
end