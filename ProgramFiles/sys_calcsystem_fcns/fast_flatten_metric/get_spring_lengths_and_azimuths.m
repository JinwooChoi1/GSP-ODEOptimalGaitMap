function [L, az, delta] = get_spring_lengths_and_azimuths(springs,x,y)

	% Get the x and y displacements of each spring's endpoints
	delta_x = x(springs(:,1))-x(springs(:,2));
	delta_y = y(springs(:,1))-y(springs(:,2));
	
	% Get the length of each spring
	L = sqrt(delta_x.^2 + delta_y.^2);
	
	% Get the azimuth of each spring
	az = atan2(delta_y,delta_x);
    
    % Turn the deltas into vectors
    delta = [delta_x delta_y];

end