function Z = sigage5(Xdist, Yage, Vref, B, b, S, d, Yref, A)
    dref = 100;
    Z = zeros(size(Xdist));
    for i = 1:numel(Xdist)
        Z(i) = (Vref - S*log10(Xdist(i)/dref) + B*((1+(dref/d)^b)/(1+(Xdist(i)/d)^b)))*(1-A*exp(Yage(i)/Yref));
    end
end