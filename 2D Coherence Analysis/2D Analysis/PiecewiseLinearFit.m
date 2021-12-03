function [fit_result] = PiecewiseLinearFit(data,data_starts,data_ends)
%assume data is of form: x,y,z, and we're fitting along z.
dimX = size(data,1);
dimY = size(data,2);
Ls = length(data_starts);
Le = length(data_ends);
assert(Ls == Le);
L = Ls;
fit_result = 0*data;

for k=1:dimX;
    for j=1:dimY;
        for q=1:L;
%             l = data_ends(q)-data_starts(q)+1;
%             fit_result(k,j,data_starts(q):data_ends(q)) = polyval(polyfit(1:l,squeeze(data(k,j,data_starts(q):data_ends(q)))',2),1:l);            
              fit_result(k,j,data_starts(q):data_ends(q)) = smooth(data(k,j,data_starts(q):data_ends(q)),20);
        end
    end
end
