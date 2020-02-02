function [outputs] = easysim(net,nT,useGpu,varargin)
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

if (~isempty(p.Results.sim_dir))
    if (~exist(p.Results.sim_dir,'dir'))
        mkdir(p.Results.sim_dir)
    end
end
sim_dir = p.Results.sim_dir;

switch class(net)
    case 'EVLIFnetwork'
        if (net.is_plastic)
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
            tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
            is_plastic,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,tau_plus,tau_x,tau_minus,tau_y] = ...
            setup_plasticEVLIFNet(net,useGpu);
            outputs.dGsyn_pre = dGsyn;
        else
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
              tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
              setupEVLIFNet(net,useGpu);
        end
    case 'AEVLIFnetwork'
        if (net.is_plastic)
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
            tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
            is_plastic,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,tau_plus,tau_x,tau_minus,tau_y] = ...
            setup_plasticAEVLIFNet(net,useGpu);
            outputs.dGsyn_pre = dGsyn;
        else
            [V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
              tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
              setupAEVLIFNet(net,useGpu);
        end
end
if (useGpu)
    if (~any(strcmp(p.UsingDefaults,'timeStepSize')))
        dt = single(p.Results.timeStepSize);
    end
    switch class(net)
        case 'EVLIFnetwork'
            if (p.Results.recompile)
                if (net.is_plastic)
                    compile_loopUpdate_plasticEVLIFNetGPU(size(V,1),length(spikeGenProbs),length(cells2record));
                else
                    compile_loopUpdateEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
                end
            end
            
            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            if (net.is_plastic)
                fprintf('Calling loopUpdate_plasticEVLIFNetGPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                plasticity_type = p.Results.plasticity_type;
                tic;
                loopUpdate_plasticEVLIFNetGPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
                                              VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,tau_synI,...
                                              Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
                                              is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
                                              tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
                t=toc;
                disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
                fclose(spkfid);
            else
                fprintf('Calling loopUpdateEVLIFNetGPU_fast_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                tic;
                dGsyn=loopUpdateEVLIFNetGPU_fast_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
                t=toc;
                disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
                fclose(spkfid);
            end
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
            cells2record = gather(cells2record);
        case 'AEVLIFnetwork'
            if (p.Results.recompile)
                if (net.is_plastic)
                    compile_loopUpdate_plasticAEVLIFNetGPU(size(V,1),length(spikeGenProbs),length(cells2record));
                else
                    compile_loopUpdateAEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
                end
            end

            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            if (net.is_plastic)
                fprintf('Calling loopUpdated_plasticAEVLIFNetGPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                plasticity_type = p.Results.plasticity_type;
                tic;
                dGsyn=loopUpdate_plasticAEVLIFNetGPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
                                               Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,tau_synI,...
                                               Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
                                               is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
                                               tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
                t=toc;
                disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
                fclose(spkfid);
            else
                fprintf('Calling loopUpdateAEVLIFNetGPU_fast_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                tic;
                loopUpdateAEVLIFNetGPU_fast_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                                dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                                dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
                t=toc;
                disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
                fclose(spkfid);
            end
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
            cells2record = gather(cells2record);
    end
else
    if (~any(strcmp(p.UsingDefaults,'timeStepSize')))
        dt = p.Results.timeStepSize;
    end
    switch class(net)
        case 'EVLIFnetwork'
            if (p.Results.recompile)
                if (p.Results.recordVars)
                    if (net.is_plastic)
                        compile_loopUpdate_plasticEVLIFNetCPU_recordVars();
                    else
                        compile_loopUpdateEVLIFNetCPU_recordVars();
                    end
                else
                    if (net.is_plastic)
                        compile_loopUpdate_plasticEVLIFNetCPU();
                    else
                        compile_loopUpdateEVLIFNetCPU();
                    end
                end
            end

            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            
            if (p.Results.recordVars)
                if (net.is_plastic)
                else
                    fprintf('Calling loopUpdateEVLIFNetCPU_recordV_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                    tic;
                    [recordV,recordVth,iappRecord]=loopUpdateEVLIFNetCPU_recordVars_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                        dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                        dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
                    t=toc;
                    outputs.voltage_recording = recordV;
                    outputs.threshold_recording = recordVth;
                    outputs.iapp_recording = iappRecord;
                end
            else
                if (net.is_plastic)
                    plasticity_type = p.Results.plasticity_type;
                    fprintf('Calling loopUpdate_plasticEVLIFNetCPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                    tic;
                    dGsyn = loopUpdate_plasticEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
                                                VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,tau_synI,...
                                                Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
                                                is_plastic,plasticity_type,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
                                                tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
                    t=toc;
                else
                    fprintf('Calling loopUpdateEVLIFNetCPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                    tic;
                    loopUpdateEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                        dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                        dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
                    t=toc;
                end
            end
            disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
            fclose(spkfid);
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
        case 'AEVLIFnetwork'
            if (p.Results.recompile)
                if (p.Results.recordVars)
                    if (net.is_plastic)
                        compile_loopUpdate_plasticAEVLIFNetCPU_recordVars();
                    else
                        compile_loopUpdateAEVLIFNetCPU_recordVars();
                    end
                else
                    if (net.is_plastic)
                        compile_loopUpdate_plasticAEVLIFNetCPU();
                    else
                        compile_loopUpdateAEVLIFNetCPU();
                    end
                end
            end

            spkfid = fopen([p.Results.sim_dir '/' p.Results.spikefile],'W');
            
            if (p.Results.recordVars)
                if (net.is_plastic)
                else
                    fprintf('Calling loopUpdateAEVLIFNetCPU_recordV_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                    tic;
                    [recordV,recordVth,iappRecord]=loopUpdateAEVLIFNetCPU_recordVars_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                        dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                        dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
                    t=toc;
                    outputs.voltage_recording = recordV;
                    outputs.threshold_recording = recordVth;
                    outputs.iapp_recording = iappRecord;
                end
            else
                if (net.is_plastic)
                    plasticity_type = p.Results.plasticity_type;
                    fprintf('Calling loopUpdate_plasticAEVLIFNetCPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                    tic;
                    dGsyn=loopUpdate_plasticAEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,...
                                                   Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,...
                                                   maxGsynE,maxGsynI,dGsyn,tau_synE,tau_synI,...
                                                   Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,...
                                                   spikeGenProbs,cells2record,is_plastic,plasticity_type,...
                                                   C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
                                                   tau_plus,tau_x,tau_minus,tau_y,nT,spkfid);
                    t=toc;
                else
                    fprintf('Calling loopUpdateAEVLIFNetCPU_mex for %1$i timesteps with dt = %2$e. \nSpikes will be saved in %3$s\n',nT,dt,[p.Results.sim_dir '/' p.Results.spikefile])
                    tic;
                    loopUpdateAEVLIFNetCPU_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                                        dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                                        dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
                    t=toc;
                end
            end
            disp(['Total sim time: ' num2str(t) '. Time per timestep = ' num2str(t/nT) ' --> ' num2str((t/nT)/dt) 'x real time'])
            fclose(spkfid);
            save([p.Results.sim_dir '/net.mat'],'net','-mat')
            save([p.Results.sim_dir '/dGsyn.mat'],'dGsyn','-mat')
    end
end
outputs.dt = dt;
outputs.cells2record = cells2record;
end

