function [outputs] = easysim(net,nT,useGpu,varargin)
% easysim(net,nT,useGpu,varargin)
%   This function funnels the input network to the correct simulation
%   function, recompiling code as requested. It also opens and closes the
%   files in which spike times are saved and also saves the connection
%   matrix used in the simulation.
%
%   easysim(net,nT,useGpu,varargin)
%       Inputs:
%           net      - network object (e.g. EVLIFnetwork or AEVLIFnetwork)
%
%           nT       - # of timesteps to simulate
%
%           useGpu   - true or false. If true, easysim will attempt to simulate the
%                      network on a GPU using compiled CUDA code
%
%           varargin - There are several optional inputs that control if and where
%                      spike times are recorded as well as whether the simulation
%                      code should be recompiled. Each optional input must be
%                      specified as a key-value pair
%                      e.g.
%                      easysim(net,nT,useGpu,'recompile',true,'spikefile','myspikes.bin')
%
%               The optional parameters are:
%                   1. 'sim_dir'    - path to directory in which to save simulation
%                                     data. The default is ./tmp
%                   2. 'spikefile'  - filename in which spike times will be
%                                     stored. Default is no spikefile.
%                   3. 'recompile'  - true or false. If true, easysim will
%                                     call the relevant compilation
%                                     functions to recompile the simulation
%                                     code to be appropriate for the current
%                                     network. NOTE: you will need to do
%                                     this anytime you are using the GPU and
%                                     the # of neurons, spikeGenerating
%                                     neurons, or # of cells to be recorded
%                                     has changed from the last call to
%                                     easysim. Default value is false.
%                   4. 'recordVars' - true or false. If true, time series
%                                     of certain variables will be returned
%                                     (e.g. membrane voltage, spike
%                                     threshold, applied current). Default
%                                     value is false.
%                   5. 'timeStepSize' - a float value giving the desired
%                                       time step size (dt) to use instead
%                                       of the automatically determined one
%                                       (which is 10x smaller than the
%                                       smallest time constant). The
%                                       default value is the automatically
%                                       determined one in the setup
%                                       function.
%                   6. 'plasticity_type' - either 'all-to-all' or
%                                          'nearest'. 'all-to-all' will
%                                          include all spike pairs in
%                                          determining changes in synaptic
%                                          strength whereas 'nearest' only
%                                          includes the most recent spikes.
%                                          Default value is nearest.
p=inputParser;
isNetwork = @(x) isa(net,'EVLIFnetwork') || isa(net,'AEVLIFnetwork');
isPositiveNumber = @(x) isnumeric(x) && ~isinf(x) && ~isnan(x) && (x>0);
addRequired(p,'net',isNetwork)
addRequired(p,'nT',isPositiveNumber)
addRequired(p,'useGpu',@islogical)
addParameter(p,'sim_dir','tmp',@ischar)
addParameter(p,'spikefile','',@ischar)
addParameter(p,'recompile',false,@islogical)
addParameter(p,'recordVars',false,@islogical)
addParameter(p,'timeStepSize',.0001,isPositiveNumber)
addParameter(p,'plasticity_type','nearest',@ischar)
parse(p,net,nT,useGpu,varargin{:})
outputs.run_params = p.Results;
plasticity_type = p.Results.plasticity_type;

% make the directory to save simulation data if it does not exist
if (~isempty(p.Results.sim_dir))
    if (~exist(p.Results.sim_dir,'dir'))
        mkdir(p.Results.sim_dir)
    end
end
sim_dir = p.Results.sim_dir;

% Setup/initialize all relevant variables. Unused variables are set to []
[V,Vreset,Cm,Gl,El,Vth,Vth0,Vth_max,tau_ref,dth,p0,GsynE,GsynI,VsynE,VsynI,tau_synE,tau_synI,...
          Iapp,std_noise,GsynMax,Isra,tau_sra,a,b,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
          ecells,icells,spikeGenProbs,cells2record,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
          tau_plus,tau_x,tau_minus,tau_y,is_plastic,C,dt] = ...
          setupNet(net,useGpu);
outputs.GsynMax_pre = GsynMax;

if (p.Results.recompile)
    compileSimulator(net,useGpu,length(cells2record));
end

spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
if (useGpu)
    if (~any(strcmp(p.UsingDefaults,'timeStepSize')))
        dt = single(p.Results.timeStepSize);
    end
    switch class(net)
        case 'EVLIFnetwork'
            tic;
            [GsynMax] = runEVLIFNetGPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
              VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
            t=toc;
        case 'AEVLIFnetwork'
            tic;
            [GsynMax,V,Vth,Isra,GsynE,GsynI,D,F,r1,r2,o1,o2] = runAEVLIFNetGPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
              Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
            t=toc;
    end
    fclose(spkfid);
    disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
    save([p.Results.sim_dir '/net.mat'],'net','-mat')
    save([p.Results.sim_dir '/GsynMax.mat'],'GsynMax','-mat')
    cells2record = gather(cells2record);
else
    if (~any(strcmp(p.UsingDefaults,'timeStepSize')))
        dt = p.Results.timeStepSize;
    end
    switch class(net)
        case 'EVLIFnetwork'
            recordVars = p.Results.recordVars;
            tic;
            [GsynMax,recordV,recordVth,recordIsyn,recordIapp,recordD,recordF] = runEVLIFNetCPU_mex(V,Vreset,...
              tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid,recordVars);
            t=toc;         
        case 'AEVLIFnetwork'
            recordVars = p.Results.recordVars;
            tic;
            size(D)
            size(F)
            [GsynMax,recordV,recordVth,recordIsyn,recordIapp,recordD,recordF] = runAEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
              Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,GsynMax,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
              p0,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid,recordVars);
            t=toc;
    end
    outputs.GsynMax = GsynMax;
    outputs.V_recording = recordV;
    outputs.Vth_recording = recordVth;
    outputs.Isyn_recording = recordIsyn;
    outputs.Iapp_recording = recordIapp;
    outputs.D_recording = recordD;
    outputs.F_recording = recordF;
    disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
    fclose(spkfid);
    save([p.Results.sim_dir '/net.mat'],'net','-mat')
    save([p.Results.sim_dir '/GsynMax.mat'],'GsynMax','-mat')
end
outputs.dt = dt;
outputs.cells2record = cells2record;
end

