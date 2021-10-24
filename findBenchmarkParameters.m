load('measurePoints.mat');
load('measureResults.mat');
p0=[0;0;0];
eps=1.0e-6;
alpha0=1.0e-3;
beta=100;
myP=levMarqDescent(@(p)leastSquaresObjective(@benchmarkModel,p,measurePoints,measureResults), p0, eps, alpha0, beta, true);
disp(myP);
  