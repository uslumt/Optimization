function [] = solveBenchmark(verbose)
  %SOLVEBENCHMARK is a script to solve the benchmark problem of Optimization for Engineers
  
  %% Purpose: 
  % final script to test all previous programming homework
  
  %% Input Definition:
  % verbose: bool, if set to true, verbose information is displayed
  
  %% Output Definition:
  % GMP: column vector in R^3, GMP of the benchmark problem
  % GMPValue: real value, objective evaluation at the GMP
  
  %% Required files:
  % all programming homework files
  
  %% Test cases:
  % <not available>
  
  %% Input verification:
  if nargin < 1
    verbose = false;
  end
  
  %% Implementation:
  if verbose
    disp('Start solveBenchmark...');
  end  

  load('measurePoints.mat');
  load('measureResults.mat');
  p0=[0;0;0];
  eps=1.0e-3;
  delta=1.0e-6;
  alpha0=1.0e-3;
  beta=100;
  myP=levMarqDescent(@(p)leastSquaresObjective(@benchmarkModel,p,measurePoints,measureResults), p0, eps, alpha0, beta, verbose);
  
  myObjective=@(x)benchmarkObjective(x,myP);
    
  myBoxConstraints=@(x)projectionInBox(x,[0;-4;-1],[8;4;1],delta);

  RandomPointInBox=[8*rand(1);8*rand(1)-4;2*rand(1)-1]; 

  disp('For the starting point');
  disp(RandomPointInBox);
     
  A = [2, 0, 0; 0, 2, 0; 0, 0, 2];
  b = [-8;0;0];
  c = 7;
  myQuadraticConstraint = @(x)quadraticConstraint(x,A,b,c);
  
  [LMPonOmega, lambda] = augmentedLagrangianDescent(@(x)myObjective(x),@(x)myBoxConstraints(x),@(x)myQuadraticConstraint(x), RandomPointInBox, 0, eps, delta, verbose);
  LMPValue = myObjective(LMPonOmega);
  
  disp('the benchmark problem is solved by');
  disp(LMPonOmega);
  disp('and the objective value is');
  disp(LMPValue);
end




