function [hhat_miny] = height_miny(yc, uc, z)
l = length(uc);
j=1; j1=1;
hhat_miny = zeros(l,1);
for i=1:l-1
    if uc(i) < 0 && uc(i+1) >= 0
        k(j)   = i;
        k1(j1) = i;
        j      = j+1;
        j1     = j1+1;
    end
    if uc(i) > 0 && uc(i+1) <= 0
        k(j)  = i;
        k2(j) = i;
        j     = j+1;
    end
end
period = 2*k(end)/length(k);  % Estimate cantilever oscillation period
for i = 1:floor(l/period)-1
    hhat_miny(round(k2(1)+(i-0.75)*period)) = ...
        min(yc(round((i-0.75)*period)+1:round((i+0.25)*period)))...
        - z(round(k2(1)+(i-0.75)*period));
end

for i=2:l
    if hhat_miny(i) == 0
        hhat_miny(i) = hhat_miny(i-1);
    end
end
end