% readSpikes debugger
nTimeSteps = 2000;
cells = 1:2000;
maxSpikeFrac = .5;
makeFakeSpikeData(nTimeSteps,cells,maxSpikeFrac)
tic;
data = readSpikes('results/test_spikes.bin',cells);
time1 = toc;
disp(time1)

tic;
data = readSpikes2('results/test_spikes.bin',cells);
time2 = toc;
disp(time2)
disp(['Time ratio = ' num2str(time1/time2)])