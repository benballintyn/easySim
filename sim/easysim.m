function [] = easysim(net,nT,useGpu,varargin)

p=inputParser;
isNetwork = @(x) isa(net,'EVLIFnetwork') || isa(net,'AEVLIFnetwork');
isPositiveNumber = @(x) isnumeric(x) && ~isinf(x) && ~isnan(x) && (x>0);
addRequired(p,'net',isNetwork)
addRequired(p,'nT',isPositiveNumber)
addRequired(p,'useGpu',@islogical)
addParameter(p,'spikefile','',@ischar)
addParameter(p,'recompile',false,@islogical)
parse(p,net,nT,useGpu,varargin{:})


switch class(net)
    case 'EVLIFnetwork'
        [V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
          setupEVLIFNet(net,useGpu);
      
        size(V,1)
        length(spikeGenProbs)
        length(cells2record)
        if (p.Results.recompile)
            compile_loopUpdateEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
        end
        
        spkfid = fopen(p.Results.spkfile,'W');
        fprintf('Calling loopUpdateEVLIFNetGPU_fast_mex for %i timesteps. Spikes will be saved in %s',nT,p.Results.spkfile)
        loopUpdateEVLIFNetGPU_fast_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
        fclose(spkfid);
         
    case 'AEVLIFnetwork'
        [V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
          setupAEVLIFNet(net,useGpu);
      
        if (p.Results.recompile)
            compile_loopUpdateAEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
        end
        
        spkfid = fopen(p.Results.spkfile,'W');
        fprintf('Calling loopUpdateEVLIFNetGPU_fast_mex for %i timesteps. Spikes will be saved in %s',nT,p.Results.spkfile)
        loopUpdateAEVLIFNetGPU_fast_mex(V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
        fclose(spkfid);
end
end

