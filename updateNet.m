function [V,Gref,GsynE,GsynI] = updateNet(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells)
spiked = (V > Vth);
e_spiked = spiked.*ecells;
i_spiked = spiked.*icells;
downstreamE = dGsyn(:,e_spiked);
downstreamI = dGsyn(:,i_spiked);
dGsynE_sum = sum(downstreamE,2);
dGsynI_sum = sum(downstreamI,2);
dGsynEdt = -GsynE/tau_synE;
dGsynIdt = -GsynI/tau_synI;
GsynE = GsynE + dGsynEdt*dt;
GsynI = GsynI + dGsynIdt*dt;
GsynE(e_spiked) = GsynE(e_spiked) + dGsynE_sum;
GsynI(i_spiked) = GsynI(i_spiked) + dGsynI_sum;
Isyn = GsynE.*(VsynE - V) + GsynI.*(VsynI - V);
dGrefdt = -Gref./tau_ref;
Gref = Gref + dGrefdt*dt;
Gref(spiked) = Gref(spiked) + dGref;
dVdt = (1./Cm)*(Gl.*(El - V + dth.*exp((V - Vth)./dth)) + Gref.*(Ek - V) + Isyn + Iapp);
V = V + dVdt*dt;
end

