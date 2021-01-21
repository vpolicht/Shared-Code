function [fit_result] = PiecewiseLinearFit2D(data,data_starts,data_ends)
%assume data is of form: x,y, and we're fitting along dim2.
dimX = size(data,1);
Ls = length(data_starts);
Le = length(data_ends);
assert(Ls == Le);
L = Ls;
fit_result = 0*data;

for k=1:dimX;
    for q=1:L;
        l = data_ends(q)-data_starts(q)+1;
        fit_result(k,data_starts(q):data_ends(q)) = polyval(polyfit(1:l,squeeze(data(k,data_starts(q):data_ends(q))),1),1:l);            
    end
end
