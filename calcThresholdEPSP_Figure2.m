function ThresholdEPSP_Figure2 = calcThresholdEPSP(EPSPsizeBins,probability,nEventsPerBin)

w = nEventsPerBin/sum(nEventsPerBin);

ft = fittype('1/(1+exp((thalf-x)/b))');
opt = fitoptions('method','NonlinearLeastSquares','Robust','off','Alg','trust',...
    'Display','off','Startpoint',[],'Lower',[0 0],'upper',[],'weight',w);


[F,GoF,Out]=fit(EPSPsizeBins(:),probability(:),ft,opt);

ThresholdEPSP = f.thalf;


end