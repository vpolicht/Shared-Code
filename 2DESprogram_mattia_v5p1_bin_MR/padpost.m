function [ B ] = padpost( PhIntCal_filt,padsize)
A=zeros(padsize,1);
B=[PhIntCal_filt;A];

end

