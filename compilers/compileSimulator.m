function [] = compileSimulator(net,useGpu,nCells2record)

if (isa(net,'EVLIFnetwork'))
    if (useGpu)
        compile_runEVLIFNetGPU(net,nCells2record);
    else
        compile_runEVLIFNetCPU();
    end 
elseif (isa(net,'AEVLIFnetwork'))
    if (useGpu)
        compile_runAEVLIFNetGPU(net,nCells2record);
    else
        compile_runEVLIFNetCPU();
    end
    error('Unable to compile. Class of network not recognized')
end
end