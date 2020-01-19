function [] = loopUpdateNetGPU_fast(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells,nT) %#codegen

coder.gpu.kernelfun;
% V=gpuArray(V);
% Gref = gpuArray(Gref);
% dGref = gpuArray(dGref);
% tau_ref = gpuArray(tau_ref);
% Vth = gpuArray(Vth);
% VsynE = gpuArray(VsynE);
% VsynI = gpuArray(VsynI);
% GsynE = gpuArray(GsynE);
% GsynI = gpuArray(GsynI);
% maxGsynE = gpuArray(maxGsynE);
% maxGsynI = gpuArray(maxGsynI);
% dGsyn = gpuArray(dGsyn);
% tau_synE = gpuArray(tau_synE);
% tau_synI = gpuArray(tau_synI);
% Cm = gpuArray(Cm);
% Gl = gpuArray(Gl);
% El = gpuArray(El);
% Ek = gpuArray(Ek);
% dth = gpuArray(dth);
% Iapp = gpuArray(Iapp);
% dt = gpuArray(dt);
% ecells = gpuArray(ecells);
% icells = gpuArray(icells);
fmt1 = '%i ';
for i=1:nT
    spiked = (V > Vth);
    %spiked = find(V > Vth);
    areSpikes = any(spiked);
    if (areSpikes)
        V(spiked) = -.08;
    end

    e_spiked = logical(spiked.*ecells);
    %e_spiked = intersect(spiked,ecells);
    i_spiked = logical(spiked.*icells);
    %i_spiked = intersect(spiked,icells);
    dGsynEdt = arrayfun(@rdivide,-GsynE,tau_synE);
    dGsynIdt = arrayfun(@rdivide,-GsynI,tau_synI);
    %GsynE = GsynE + dGsynEdt*dt;
    %GsynI = GsynI + dGsynIdt*dt;
    GsynE = arrayfun(@plus,GsynE,dGsynEdt*dt);
    GsynI = arrayfun(@plus,GsynI,dGsynIdt*dt);

    %dGsynE_sum = sum(dGsyn.*e_spiked,2);
    %dGsynI_sum = sum(dGsyn.*i_spiked,2);
    if (areSpikes)
        dGsynE_sum = sum(dGsyn(:,e_spiked),2);
        dGsynI_sum = sum(dGsyn(:,i_spiked),2);
        GsynE = arrayfun(@plus,GsynE,dGsynE_sum);
        GsynE = arrayfun(@min,GsynE,maxGsynE);
        GsynI = arrayfun(@plus,GsynI,dGsynI_sum);
        GsynI = arrayfun(@min,GsynI,maxGsynI);
    end
    %GsynE = GsynE + dGsynE_sum
    %GsynI = GsynI + dGsynI_sum;
    
    Isyn = arrayfun(@times,GsynE,arrayfun(@minus,VsynE,V)) + arrayfun(@times,GsynI,arrayfun(@minus,VsynI,V));
    %dGrefdt = arrayfun(@rdivide,-Gref,tau_ref);
    %dGrefdt = -Gref./tau_ref;
    %Gref = Gref + dGrefdt*dt;
    %Gref(spiked) = Gref(spiked) + dGref(spiked);
    %Gref = Gref + dGref.*spiked;
    %Gref = Gref + arrayfun(@times,dGref,spiked);
    %dVdt = (1./Cm).*(Gl.*(El - V + dth.*exp((V - Vth)./dth)) + Gref.*(Ek - V) + Isyn + Iapp);
    %dVdt = (1./Cm).*(Gl.*(El - V + dth.*exp((V - Vth)./dth)) + Isyn + Iapp);
    f1 = (1./Cm);
    f2 = arrayfun(@minus,El,V);
    f3 = arrayfun(@minus,V,Vth);
    f4 = arrayfun(@rdivide,f3,dth);
    f5 = arrayfun(@exp,f4);
    f6 = arrayfun(@times,dth,f5);
    f7 = arrayfun(@plus,f2,f6);
    f8 = arrayfun(@times,Gl,f7);
    f9 = arrayfun(@plus,f8,Isyn);
    f10 = arrayfun(@plus,f9,Iapp);
    dVdt = arrayfun(@times,f1,f10);
    %{
    f1 = (1./Cm);
    f2 = arrayfun(@minus,El,V);
    f3 = arrayfun(@minus,V,Vth);
    f4 = arrayfun(@rdivide,f3,dth);
    f5 = arrayfun(@exp,f4);
    f6 = arrayfun(@times,dth,f5);
    f7 = arrayfun(@plus,f2,f6);
    f8 = arrayfun(@times,Gl,f7);
    f9 = arrayfun(@minus,Ek,V);
    f10 = arrayfun(@times,Gref,f9);
    f11 = arrayfun(@plus,f8,f10);
    f12 = arrayfun(@plus,f11,Isyn);
    f13 = arrayfun(@plus,f12,Iapp);
    dVdt = arrayfun(@times,f1,f13);
    %}
    V = V + dVdt*dt;
    %{
    if (areSpikes)
        fprintf(spkfid,fmt1,int16(i));
        spkInds = find(spiked);
        for j=1:length(spkInds)
            fprintf(spkfid,fmt1,int16(spkInds(j)));
        end
        fprintf(spkfid,'\n');
    end
    %}
    fprintf('%i',int32(i))
end
end