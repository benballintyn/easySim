function [dt,cells2record,sim_dir] = easysim(net,nT,useGpu,varargin)
% easysim(net,nT,useGpu,varargin)
%   This function funnels that input network to the correct simulation
%   function, recompiling code as requested. It also opens and closes the
%   files in which spike times are saved and also saves the connection
%   matrix used in the simulation.
%
%   net      - network object (e.g. EVLIFnetwork or AEVLIFnetwork)
%
%   nT       - # of timesteps to simulate
%
%   useGpu   - true or false. If true, easysim will attempt to simulate the
%              network on a GPU using compiled CUDA code
%
%   varargin - There are a few optional inputs that control if and where
%              spike times are recorded as well as whether the simulation
%              code should be recompiled. Each optional input must be
%              specified as a key-value pair
%              e.g.
%              easysim(net,nT,useGpu,'recompile',true,'spikefile','myspikes.bin')
%
%               The optional parameters are:
%                   1. 'sim_dir'   - path to directory in which to save simulation
%                                  data. The default is ./tmp
%                   2. 'spikefile' - filename in which spike times will be
%                                    stored.
%                   3. 'recompile' - true or false. If true, easysim will
%                                    call the relevant compilation
%                                    functions to recompile the simulation
%                                    code to be appropriate for the current
%                                    network. NOTE: you will need to do
%                                    this anytime you are using the GPU and
%                                    the # of neurons, spikeGenerating
%                                    neurons, or # of cells to be recorded
%                                    has changed from the last call to
%                                    easysim.
p=inputParser;
isNetwork = @(x) isa(net,'EVLIFnetwork') || isa(net,'AEVLIFnetwork');
isPositiveNumber = @(x) isnumeric(x) && ~isinf(x) && ~isnan(x) && (x>0);
addRequired(p,'net',isNetwork)
addRequired(p,'nT',isPositiveNumber)
addRequired(p,'useGpu',@islogical)
addParameter(p,'sim_dir','tmp',@ischar)
addParameter(p,'spikefile','',@ischar)
addParameter(p,'recompile',false,@islogical)
parse(p,net,nT,useGpu,varargin{:})

if (~isempty(p.Results.sim_dir))
    if (~exist(p.Results.sim_dir,'dir'))
        mkdir(p.Results.sim_dir)
    end
end
sim_dir = p.Results.sim_dir;

if (useGpu)
    switch class(net)
        case 'EVLIFnetwork'
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
              tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
              setupEVLIFNet(net,useGpu);

            if (p.Results.recompile)
                compile_loopUpdateEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
            end

            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            fprintf('Calling loopUpdateEVLIFNetGPU_fast_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
            tic;
            loopUpdateEVLIFNetGPU_fast_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
            t=toc;
            disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
            fclose(spkfid);
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
            cells2record = gather(cells2record);
        case 'AEVLIFnetwork'
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
              tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
              setupAEVLIFNet(net,useGpu);

            if (p.Results.recompile)
                compile_loopUpdateAEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
            end

            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            fprintf('Calling loopUpdateAEVLIFNetGPU_fast_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
            tic;
            loopUpdateAEVLIFNetGPU_fast_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
            t=toc;
            disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
            fclose(spkfid);
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
            cells2record = gather(cells2record);
    end
else
    switch class(net)
        case 'EVLIFnetwork'
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
              tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
              setupEVLIFNet(net,useGpu);

            if (p.Results.recompile)
                compile_loopUpdateEVLIFNetCPU();
            end

            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            fprintf('Calling loopUpdateEVLIFNetCPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
            tic;
            loopUpdateEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
            t=toc;
            disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
            fclose(spkfid);
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
        case 'AEVLIFnetwork'
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
              tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
              setupAEVLIFNet(net,useGpu);

            if (p.Results.recompile)
                compile_loopUpdateAEVLIFNetCPU();
            end

            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            fprintf('Calling loopUpdateAEVLIFNetCPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
            tic;
            loopUpdateAEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
            t=toc;
            disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
            fclose(spkfid);
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
    end
end
end

