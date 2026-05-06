function [mean_deg, std_deg] = circular_stats_deg(angles_deg)
% CIRCULAR_STATS_DEG  Circular mean and std for orientation data (degrees)
%   [MEAN_DEG, STD_DEG] = CIRCULAR_STATS_DEG(ANGLES_DEG)
%   Computes circular statistics for orientation data where the period is 180°.
%   ANGLES_DEG may be any vector of angles (degrees). MEAN_DEG is returned in [0,180).
%   STD_DEG is the circular standard deviation in degrees (orientation-aware).

if isempty(angles_deg)
    mean_deg = NaN; std_deg = NaN; return;
end

% convert to radians and double angles for orientation (period pi -> use 2*theta)
theta = deg2rad(angles_deg(:));
z = exp(1i * 2 * theta);
R = mean(z);
ang2 = angle(R); % doubled mean angle in radians
mean_rad = 0.5 * ang2; % map back
mean_deg = mod(rad2deg(mean_rad), 180);
% resultant length
Rlen = abs(R);
% circular std (radians) for doubled angles: S = sqrt(-2*log(Rlen))/2
if Rlen <= 0
    std_deg = NaN;
else
    std_rad = 0.5 * sqrt(max(0, -2 * log(Rlen)));
    std_deg = rad2deg(std_rad);
end
end