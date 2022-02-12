function p = make_ellipse_gait(r)
r_len = length(r);
p.phi_def = cell(1,2);
p.dphi_def = cell(1,2);
p.ddphi_def = cell(1,2);
if r_len == 2
    p.phi_def{1} = @(t) r(1)*cos(t);
    p.phi_def{2} = @(t) r(2)*cos(t-pi/4);
    p.dphi_def{1} = @(t) -r(1)*sin(t);
    p.dphi_def{2} = @(t) -r(2)*sin(t-pi/4);
    p.ddphi_def{1} = @(t) -r(1)*cos(t);
    p.ddphi_def{2} = @(t) -r(2)*cos(t-pi/4);
elseif r_len == 3
    p.phi_def{1} = @(t) r(1)*cos(t);
    p.phi_def{2} = @(t) r(2)*cos(t-r(3));
    p.dphi_def{1} = @(t) -r(1)*sin(t);
    p.dphi_def{2} = @(t) -r(2)*sin(t-r(3));
    p.ddphi_def{1} = @(t) -r(1)*cos(t);
    p.ddphi_def{2} = @(t) -r(2)*cos(t-r(3));
end

end