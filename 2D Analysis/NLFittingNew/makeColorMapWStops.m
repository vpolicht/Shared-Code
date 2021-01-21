function cMap = makeColorMapWStops(stops,colors)
    M = length(stops);
    cellcmap = cell(1);
    for k=1:(M-1)
        cellcmap{k} = interpMap(colors(k,:), colors(k+1,:),floor(stops(k+1)-stops(k)));
    end
    cMap = cat(1,cellcmap{:});
end
    

function cMap = interpMap(colorStart, colorEnd, n)
    cMap = zeros(n,3);
    for i = 1:3
        cMap(1:n,i) = linspace(colorStart(i), colorEnd(i), n);
    end
end