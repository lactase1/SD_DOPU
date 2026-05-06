function rgb = generate_hsv_oa(oa_deg)
% generate_hsv_oa  Map OA angles (0-180 deg) to RGB using HSV hue
%   rgb = generate_hsv_oa(oa_deg)
%   oa_deg can be scalar or array (degrees in [0,180])
H = double(oa_deg) / 180;    % normalize to [0,1]
S = ones(size(H));
V = ones(size(H));
hsv = cat(3, H, S, V);
rgb = hsv2rgb(hsv);
end