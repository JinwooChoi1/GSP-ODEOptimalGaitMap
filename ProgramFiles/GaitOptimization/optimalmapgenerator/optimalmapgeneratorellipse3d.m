function y=optimalmapgeneratorellipse3d(s,dimension,npoints,a,lb,ub,stretch,direction,costfunction,handles)
%%%%%%%%%%%%%%
% This function takes an input gait and runs fmincon to find the neareast locally
% optimal gait

%Inputs:
%
%s: System file which contains the connection vector field, CCF's and
%   metric data
%dimension: Indicates the number of shape variables of the system
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
% a: n x m array, for n waypoints in m dimensions
% lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% direction: Direction to optimize travel: 1=x,2=y,3=theta
% costfunction: costfunction type to optimize for
%           Options: Path length (metric), torque squared (metric),
%           covariant acceleration (metric), path length (coordinate),
%           acceleration (coordinate)
%
%
% Outputs:
%
% y: Matrix whose values indicate coordinates of points which form the optimal gait
%%%%%%%%%%%%

A=[];
b=[];
Aeq=[];
beq=[];

lb1=[];
ub1=[];

writerObj = [];
% % Uncomment this section if you'd like to record a video of the
% % optimizer's steps in the shape space
% writerObj = VideoWriter('cost_as_time_period_circle_start.mp4','MPEG-4');
% writerObj.FrameRate = 5;
% % set the seconds per image
% % open the video writer
% open(writerObj);
% figure(5);
% subplot(1,2,1)
% contour(s.grid.eval{1},s.grid.eval{2},s.DA_optimized{1},'LineWidth',1.5)
% axis square
% hold on

s.costfunction = costfunction;
global bestCost bestDisp bestEff stepOptimalGaits;
bestEff = 0;
stepOptimalGaits = cell(4,1);
%Suppress warning for annoying thing in jacobianeqicalculator
warning('off','MATLAB:sqrtm:SingularMatrix');

try
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','iter','Algorithm','sqp','CheckGradients',false,'FiniteDifferenceType','central','MaxIter',4000,'MaxFunEvals',20000,'TolCon',10^-3,'StepTolerance',1e-6);
catch
    error('This code requires the global optimization toolbox to run')
end

objective_function_gradient = @(x) solvedifffmincon_ellipse_gait(x,s,npoints,dimension,direction,lb,ub,writerObj);
constraint_function = @(x) nonlcon_ellipse_gait(x,s,npoints,dimension,lb,ub,direction);

[x0, ~,~,~]=fmincon(objective_function_gradient,[2;2],A,b,Aeq,beq,lb1,ub1,constraint_function,options);
odeoption=odeset('Events',@(t,x) event_optimal_map_ellipse_gait(x,s,npoints,dimension,direction));
x0=[x0;pi/2];
ode_fun = @(t,x) ode_optimal_map_ellipse_gait(x,s,npoints,dimension,direction);

[~,x,~,xe,ie]=ode45(ode_fun,[0 inf],x0,odeoption);
stepOptimalGaits{1} = x0;
t=linspace(0,2*pi,100).';
p=make_ellipse_gait(stepOptimalGaits{1});

y = [p.phi_def{1}(t) p.phi_def{2}(t)];
end


