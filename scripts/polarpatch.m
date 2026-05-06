function p = polarpatch(theta, rho, color, alpha_val)
% polarpatch  Simple polar-area-like plot (returns polarplot handle)
%   p = polarpatch(theta_rad, rho, color, alpha_val)
% Note: keeps behavior compatible with previous script (draws a thick line).
theta = [theta(:); theta(1)];
rho = [rho(:); rho(1)];
% draw on current polaraxes
p = polarplot(theta, rho, 'Color', color, 'LineWidth', 2);
end