%%%%Calc_Reso%%% function
    function [AverageFreq,GroupDelay,Resolution]=Calc_Reso(OptFreq,InterceptCal, SlopeCal,sizeRotInterf) %%%size of rotated Interferogram oaming from FFT phasing unwrap
    PseudoFreq=(OptFreq-InterceptCal)/SlopeCal; %%?pseudofrequencies
    Tp=PseudoFreq.*sizeRotInterf./OptFreq*1000;
    kk=Tp(1);
    SS=OptFreq(1);
    
for i=2:length(Tp)
TT(i-1)=(Tp(i)-kk)/(OptFreq(i)-SS);
SS=OptFreq(i);
kk=OptFreq(i);
end
TT(length(Tp))=TT(length(Tp)-1); %%%--> fare un check

Tg=TT.*OptFreq+Tp;  %? we are doing tg=dtp/dw*w+tp
GroupDelay=mean(Tg);
AverageFreq=mean(OptFreq);
Resolution=3e8/(AverageFreq^2*GroupDelay); %?Delta_nu=Delta_t*lambda*2/c
    end