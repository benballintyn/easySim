function [V,Gref,GsynE,GsynI,spiked] = updateNet(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells) %#codegen

%coder.gpu.kernelfun;

spiked = (V > Vth);
disp(['num spiked = ' num2str(sum(spiked))])
V(spiked) = -.08;
disp(['still spiking = ' num2str(sum(V > Vth))])
e_spiked = logical(spiked.*ecells);
i_spiked = logical(spiked.*icells);
dGsynEdt = -GsynE./tau_synE;
dGsynIdt = -GsynI./tau_synI;
GsynE = GsynE + dGsynEdt*dt;
GsynI = GsynI + dGsynIdt*dt;


dGsynE_sum = sum(dGsyn(:,e_spiked),2);
dGsynI_sum = sum(dGsyn(:,i_spiked),2);

GsynE = GsynE + dGsynE_sum;
GsynI = GsynI + dGsynI_sum;
Isyn = GsynE.*(VsynE - V) + GsynI.*(VsynI - V);
dGrefdt = -Gref./tau_ref;
Gref = Gref + dGrefdt*dt;
%Gref(spiked) = Gref(spiked) + dGref(spiked);
Gref = Gref + dGref.*spiked;
dVdt = (1./Cm).*(Gl.*(El - V + dth.*exp((V - Vth)./dth)) + Gref.*(Ek - V) + Isyn + Iapp);
V = V + dVdt*dt;
end