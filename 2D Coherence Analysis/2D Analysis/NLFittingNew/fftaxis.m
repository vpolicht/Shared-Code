function [r] = fftaxis(N,TSpac)
    if abs(mod(N,2))>0
        a = (0:1:floor((N-1)/2))';
        b = (-(N-1)/2:1:-1)';
        r = cat(1,a,b)/(N*TSpac);
    else
        a = (0:1:floor(N/2-1))';
        b = (-floor(N/2):1:-1)';
        r = cat(1,a,b)/(N*TSpac);
    end
end