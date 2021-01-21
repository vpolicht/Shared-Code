function [ y ] = smoothmr( vett,span )
y(1)=vett(1);
y(length(vett))=vett(length(vett));
for i=2:1:length(vett)-1
    y(i)=(vett(i-1)+vett(i)+vett(i+1))/span;
end

end

