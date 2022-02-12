# Geometric System Plotter

1. Open /ProgramFiles/GaitOptimization/optimalmapgenerator/optimalmapgenerator.m
2. Put the breakpoint after running ode45.
3. StepOptimalGaits may contain the fourier coefficient, displacement and cost at the event point.

Ex)

| Fourier Coeff| Displacement | Cost          |
| :---         |    :----:    |          ---: |
| 10x3 Double  | 0.8888       | 0.7777        |
| 10x3 Double  | 0.6666       | 0.5555        |
| 10x3 Double  | 0.5555       | 0.3333        |
| 10x3 Double  | 0.4444       | 0.1111        |

## Try this code to plot

Open the CCF Plot by sysplotter, and then

    y1=cell(4,1);
    for i = 1:4
      y1{i}=path_from_fourier(stepOptimalGaits{i,1},npoints,dimension);
    end
    hold on;
    for i = 1:4
      plot3(y1{i}(:,1),y1{i}(:,2),y1{i}(:,3),'Linewidth',4);
    end
