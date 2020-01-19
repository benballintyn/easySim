function [V,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,dth,Iapp,dt,ecells,icells] = setupEVLIFNet(net,simobj,useGpu)
% Total # of neurons to be simulated
N = net.nNeurons;

if (useGpu)
    % Variables that may change with time
    V = gpuArray(zeros(N,1,'single'));
    GsynE = gpuArray(zeros(N,1,'single'));
    GsynI = gpuArray(zeros(N,1,'single'));
    maxGsynE = gpuArray(zeros(N,1,'single'));
    maxGsynI = gpuArray(zeros(N,1,'single'));
    Iapp = gpuArray(ones(N,1,'single'))*2e-10;
    Vth = gpuArray(zeros(N,1,'single'));
    VsynE = gpuArray(ones(N,1,'single'));
    VsynI = gpuArray(ones(N,1,'single'));

    % Variables that will not change with time
    Vth0 = gpuArray(zeros(N,1,'single'));
    Vth_max = gpuArray(zeros(N,1,'single'));
    tau_ref = gpuArray(zeros(N,1,'single'));
    dGsyn = gpuArray(zeros(N,N,'single'));
    tau_synE = gpuArray(zeros(N,1,'single'));
    tau_synI = gpuArray(zeros(N,1,'single'));
    Cm  = gpuArray(zeros(N,1,'single'));
    Gl = gpuArray(zeros(N,1,'single'));
    El = gpuArray(zeros(N,1,'single'));
    dth = gpuArray(zeros(N,1,'single'));
    dt = simobj.dt;
    ecells = gpuArray(zeros(net.nNeurons,1));
    icells = gpuArray(zeros(net.nNeurons,1));
else
    % Variables that may change with time
    V = zeros(N,1);
    GsynE = zeros(N,1);
    GsynI = zeros(N,1);
    maxGsynE = zeros(N,1);
    maxGsynI = zeros(N,1);
    Iapp = ones(N,1)*2e-10;
    Vth = zeros(N,1);
    VsynE = ones(N,1);
    VsynI = ones(N,1);

    % Variables that will not change with time
    Vth0 = zeros(N,1);
    Vth_max = zeros(N,1);
    tau_ref = zeros(N,1);
    dGsyn = zeros(N,N);
    tau_synE = zeros(N,1);
    tau_synI = zeros(N,1);
    Cm  = zeros(N,1);
    Gl = zeros(N,1);
    El = zeros(N,1);
    dth = zeros(N,1);
    dt = simobj.dt;
    ecells = zeros(N,1);
    icells = zeros(N,1);
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
    Vth(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Vth0,net.groupInfo(i).std_Vth0,groupN,1);
    Vth0(preStart:preEnd) = Vth(preStart:preEnd);
    Vth_max(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Vth_max,net.groupInfo(i).std_Vth_max,groupN,1);
    VsynE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_VsynE,net.groupInfo(i).std_VsynE,groupN,1);
    VsynI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_VsynI,net.groupInfo(i).std_VsynI,groupN,1);
    maxGsynE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_max_GsynE,net.groupInfo(i).std_max_GsynE,groupN,1);
    maxGsynI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_max_GsynI,net.groupInfo(i).std_max_GsynI,groupN,1);
    tau_ref(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_ref,net.groupInfo(i).std_tau_ref,groupN,1);
    tau_synE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_synE,net.groupInfo(i).std_tau_synE,groupN,1);
    tau_synI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_synI,net.groupInfo(i).std_tau_synI,groupN,1);
    Cm(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Cm,net.groupInfo(i).std_Cm,groupN,1);
    Gl(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Gl,net.groupInfo(i).std_Gl,groupN,1);
    El(preStart:preEnd) = normrnd(net.groupInfo(i).mean_El,net.groupInfo(i).std_El,groupN,1);
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