function [wlpumpCut, interceptPhase, SlopePhase,phasecut]=linearphaseFit(wlpump,amplFFT, phase,HighWlFit,LowWlFit)
cut=find(wlpump>LowWlFit & wlpump<HighWlFit);
wlpumpCut=wlpump(cut);
amplFFTcut=amplFFT(cut);
phasecut=phase(cut);

yy=fliplr(phasecut);
xx=fliplr(300000./wlpumpCut); %%freq cut in THz

P = polyfit(xx,yy',1);
interceptPhase=P(2);
SlopePhase=P(1)/6.28*1000;
% figure(7)
% plot(xx,yy)
% hold on
% plot(xx,P(1).*xx+P(2),'r')
end