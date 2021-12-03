function [ Z ] = UncompressData( S, Bx, By)
    Z = By * S * Bx';
end

