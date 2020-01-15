function [V,Gref,GsynE,GsynI,spiked] = updateNet(V,Gref,dGref,tau_ref,Vth,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,...
                            dGsyn,tau_synE,tau_synI,Cm,Gl,El,Ek,dth,Iapp,dt,ecells,icells) %#codegen

%coder.gpu.kernelfun;

spiked = (V > Vth);
V(spiked) = -.08;

e_spiked = logical(spiked.*ecells);
i_spiked = logical(spiked.*icells);
dGsynEdt = arrayfun(@rdivide,-GsynE,tau_synE);
dGsynIdt = arrayfun(@rdivide,-GsynI,tau_synI);
GsynE = GsynE + dGsynEdt*dt;
GsynI = GsynI + dGsynIdt*dt;


%dGsynE_sum = sum(dGsyn.*e_spiked,2);
%dGsynI_sum = sum(dGsyn.*i_spiked,2);
dGsynE_sum = sum(dGsyn(:,e_spiked),2);
dGsynI_sum = sum(dGsyn(:,i_spiked),2);

GsynE = GsynE + dGsynE_sum; GsynE = arrayfun(@min,GsynE,maxGsynE);
GsynI = GsynI + dGsynI_sum; GsynI = arrayfun(@min,GsynI,maxGsynI);
Isyn = arrayfun(@times,GsynE,arrayfun(@minus,VsynE,V)) + arrayfun(@times,GsynI,arrayfun(@minus,VsynI,V));
%dGrefdt = arrayfun(@rdivide,-Gref,tau_ref);
%dGrefdt = -Gref./tau_ref;
%Gref = Gref + dGrefdt*dt;
%Gref(spiked) = Gref(spiked) + dGref(spiked);
%Gref = Gref + dGref.*spiked;
%Gref = Gref + arrayfun(@times,dGref,spiked);
%dVdt = (1./Cm).*(Gl.*(El - V + dth.*exp((V - Vth)./dth)) + Gref.*(Ek - V) + Isyn + Iapp);
dVdt = (1./Cm).*(Gl.*(El - V + dth.*exp((V - Vth)./dth)) + Isyn + Iapp);
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
end