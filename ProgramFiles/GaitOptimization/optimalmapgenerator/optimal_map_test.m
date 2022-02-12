function optimal_map_test(y,s,direction,ode_fun,nfparam)
figure(2);
plot(y(:,2),y(:,3));

[y1,y2]=meshgrid(y(:,2),y(:,3));
size_y = length(y(:,2));

disp = y1;
cost = y1;
for i = 1:size_y
    for j = 1:size_y
        temp_y =  y(i,:);
        p = path_from_fourier(y(1,:));
        [~, temp_disp, temp_cost] = evaluate_displacement_and_cost1(s,p,[0, 2*pi],'interpolated','fixed_step');
        disp(i,j) = temp_disp(direction);
        cost(i,j) = temp_cost;
    end
end

hold on;
index_y=floor(length(y)/4);
disp_lvl_list = zeros(1,5);
cost_lvl_list = zeros(1,5);
for i = 1:5
    p = make_ellipse_gait(y(index_y*(i-1)+1,:));
    [~, temp_disp, temp_cost] = evaluate_displacement_and_cost1(s,p,[0, 2*pi],'interpolated','fixed_step');
    disp_lvl_list(i) = temp_disp(direction);
    cost_lvl_list(i) = temp_cost;
end
contour(r1,r2,disp,'LevelList',disp_lvl_list);
contour(r1,r2,cost,'--','LevelList',cost_lvl_list);

end