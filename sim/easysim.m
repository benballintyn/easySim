function [dt,cells2record,sim_dir] = easysim(net,nT,useGpu,varargin)

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
    end
end
end

