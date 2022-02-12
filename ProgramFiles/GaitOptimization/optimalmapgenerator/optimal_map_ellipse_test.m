function optimal_map_ellipse_test(y,s,direction,ode_fun)
figure(2);
plot(y(:,1),y(:,2));

r1lim=xlim;
r1lim(1)=r1lim(1)*0.8;
r1lim(2)=r1lim(2)*1.2;
r2lim=ylim;
r2lim(1)=r2lim(1)*0.8;
r2lim(2)=r2lim(2)*1.2;
size_r=41;
[r1,r2]=meshgrid(linspace(r1lim(1),r1lim(2),size_r),linspace(r2lim(1),r2lim(2),size_r));

disp = r1;
cost = r1;
for i = 1:size_r
    for j = 1:size_r
        temp_r(1)=r1(i,j);
        temp_r(2)=r2(i,j);
        p = make_ellipse_gait(temp_r);
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

figure(3);
plot(y(:,1),y(:,2));
hold on;
[disp_r1,disp_r2]=gradient(disp);
[cost_r1,cost_r2]=gradient(cost);
lambda = cost_r1;
for i = 1:size_r
    for j = 1:size_r
        temp_cost_r = [cost_r1(i,j);cost_r2(i,j)];
        temp_disp_r = [disp_r1(i,j);disp_r2(i,j)];
        lambda(i,j) = pinv(temp_disp_r)*temp_cost_r;
    end
end

pdot_r1 = cost_r1-lambda.*disp_r1;
pdot_r2 = cost_r2-lambda.*disp_r2;
xdot_r1 = zeros(4,4);
xdot_r2 = zeros(4,4);
for i = 1:4
    temp_xdot = ode_fun(0,y(index_y*i,:));
    xdot_r1(i,i) = temp_xdot(1);
    xdot_r2(i,i) = temp_xdot(2);
end
streamslice(r1,r2,-pdot_r1.^3,-pdot_r2.^3);
quiver(y(index_y*(1:4),1),y(index_y*(1:4),2),xdot_r1,xdot_r2);
end