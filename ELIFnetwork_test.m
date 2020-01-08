% ELIFnetwork_test
clear all;
net = ELIFnetwork();
net.addGroup('group1',1000,'excitatory',1);
net.addGroup('group2',1000,'excitatory',2);
net.addGroup('group3',1000,'excitatory',3);

weightRange = 0:.001:1;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
randConnProb = .01;
weightFunction = @() rand();
gaussConnProbFunction = @(x)(sqrt(2)/(.2*sqrt(pi)))*exp(-x.^2/(2*.2^2));

gaussConnParams.connProbFunction = gaussConnProbFunction;
gaussConnParams.weightFunction = weightFunction;
gaussConnParams.useWrap = true;

randConnParams.connProb = randConnProb;
randConnParams.weightDistribution = uniformWeightDist;

net.connect(1,1,'gaussian',gaussConnParams);
net.connect(2,2,'gaussian',gaussConnParams);
net.connect(3,3,'gaussian',gaussConnParams);
net.connect(1,2,'random',randConnParams);
net.connect(1,3,'random',randConnParams);
net.connect(2,3,'random',randConnParams);


