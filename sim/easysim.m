function [] = easysim(net,nT,useGpu,varargin)
p=inputParser;
isNetwork = @(x) isa(net,'EVLIFnetwork') || isa(net,'AEVLIFnetwork');
isPositiveInteger = @(x) isnumeric(x) && ~isinf(x) && ~isnan(x) && (x>0) && isinteger(x);
addRequired(p,'net',isNetwork)
addRequired(p,'nT',isPositiveInteger)
addRequired(p,'useGpu',@islogical)
addParameter(p,'spikefile','',@ischar)
addParamter(p,'recompile',false,@islogical)
parse(p,net,nT,useGpu,varargin{:})


switch class(net)
    case 'EVLIFnetwork'
        [V,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
          setupEVLIFNet(net,useGpu);
        
        if (p.Results.recompile)
            compile_loopUpdateEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
        end
        
        spkfid = fopen(p.Results.spkfile,'W');
        loopUpdateEVLIFNetGPU_fast(V,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
         fclose(spkfid);
         
    case 'AEVLIFnetwork'
        [V,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
          setupAEVLIFNet(net,useGpu);
      
        if (p.Results.recompile)
            compile_loopUpdateAEVLIFNetGPU_fast(size(V,1),length(spikeGenProbs),length(cells2record));
        end
        
        spkfid = fopen(p.Results.spkfile,'W');
        loopUpdateAEVLIFNetGPU_fast(V,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,...
                            dt,ecells,icells,spikeGenProbs,cells2record,nT,spkfid);
        fclose(spkfid);
end
end

