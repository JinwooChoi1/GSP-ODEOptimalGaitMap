function p = make_ellipse_gait(r)
p.phi_def = cell(1,2);
p.dphi_def = cell(1,2);
p.ddphi_def = cell(1,2);

p.phi_def{1} = @(t) r(1)*cos(t)*cos(pi/4) - r(2)*sin(t)*sin(pi/4);
p.phi_def{2} = @(t) r(1)*cos(t)*sin(pi/4) + r(2)*sin(t)*cos(pi/4);
p.dphi_def{1} = @(t) -r(1)*sin(t)*cos(pi/4) - r(2)*cos(t)*sin(pi/4);
p.dphi_def{2} = @(t) -r(1)*sin(t)*sin(pi/4) + r(2)*cos(t)*cos(pi/4);
p.ddphi_def{1} = @(t) -r(1)*cos(t)*cos(pi/4) + r(2)*sin(t)*sin(pi/4);
p.ddphi_def{2} = @(t) -r(1)*cos(t)*sin(pi/4) - r(2)*sin(t)*cos(pi/4);
end