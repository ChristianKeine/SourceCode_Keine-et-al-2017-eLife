function activityDependentSubtractiveInhibition_Figure2S1(Raster,NTrials)

PSTH = sum(Raster);
Tau = 0.003;
I0 = 0.05;
O = 85;
S = 5;
T = [0:0.0005:0.025];

Kexp = [0,exp(-T/Tau)];
IntPSTH  = conv(PSTH,Kexp);
IntPSTH = IntPSTH(1:length(PSTH));


SubRate = I0./(1+exp(-S*(IntPSTH-prctile(IntPSTH,O))));        
NRemove = NTrials * SubRate;
NRemove(NRemove>NTrials) = NTrials;
RelInd = (NRemove-NTrials*PSTH)>0;
NRemove(RelInd) = PSTH(RelInd);
NRemove = round(NRemove);
RemInd = find(NRemove>0 & PSTH > 0);
for iB = 1:length(RemInd)
  SpikePos = find(Raster);
  Raster(RemInd(iB),SpikePos(1:NRemove(RemInd(iB)))) = 0;
end






end