function [h_sample,h_radius] = construct_sample(base,height,f_drive,scanspeed,scansize,steps,lines)

%% Sample height
h12      = height*[zeros(4,f_drive/scanspeed*scansize/steps), ones(4,f_drive/scanspeed*scansize/steps)]+base;     % Sample height in nm
h13      = height*[zeros(4,f_drive/scanspeed*scansize/steps*3), ones(4,f_drive/scanspeed*scansize/steps)]+base;   % Sample height in nm
h122     = repmat(h12,1,steps/2);
h132     = repmat(h13,1,steps/4);
h_sample = [zeros(4,f_drive/scanspeed*scansize)+base;h122;zeros(4,f_drive/scanspeed*scansize)+base;h132;zeros(4,f_drive/scanspeed*scansize)+base];

%% Sample height considering tip radius
ps_i = scansize/lines;
ps_j = scansize/length(h_sample(1,:));

tr = 1;

il = floor(2*tr/ps_i);
jl = floor(2*tr/ps_j);
clear O
for i = 0:il
    for j = 0:jl
        O(i+1,j+1) = tr - sqrt(max(0, tr^2 - ((max(0,(i-0.5))*ps_i)^2 + (max(0,(j-0.5))*ps_j)^2) ));
    end
end
O2 = [flip(flip(O(2:end,:),2),1), flip(O(2:end,2:end),1);
    flip(O,2), O(:,2:end)];

for i = 1:length(h_sample(:,1))
    for j = 1:length(h_sample(1,:))
        h_radius(i,j) = max(h_sample(max(1,i-il):min(lines,i+il),max(1,j-jl):min(length(h_sample(1,:)),j+jl)) - O2(max(1,il-i+2):min(2*il+1,lines - i+il+1),max(1,jl-j+2):min(2*jl+1,length(h_sample(1,:)) - j+jl+1)),[], 'all');
    end
end


end

