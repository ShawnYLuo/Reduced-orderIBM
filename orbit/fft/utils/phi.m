function val = phi(x, y)
    r = sqrt(x.*x + y.*y);
    
    val = (r < 1d0).*(1 + (2*r - 3).*r.*r);
    
    const = 0.3*pi;
    val = val/const;
end