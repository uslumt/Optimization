function [errorVector, jacobian] = leastSquaresObjective(model, p, xData, fData)
%LEASTSQUARESOBJECTIVE builds the least squares objective of comparing a model function with meassure data

%% Purpose:
% provides errorVector and jacobian of the least squares mapping p -> 0.5*sum_k (model(xData_k,p)-fData_k)^2

%% Input Definition:
% model: function handle of type [value, gradient_x, gradient_p] = f(x,p), objective model.
% p: column vector in R^m (parameter space)
% xData: matrix in R^nxN (measure points). xData(:,k) returns the k-th measure point as column vector.
% fData: row vector in R^1xN (measure results). fData(:,k) returns the k-th measure result as a scalar.

%% Output Definition:
% errorVector: column vector in R^N. errorVector(k,:) returns: model(xData_k,p)-fData_k
% jacobian: matrix in R^Nxm. jacobian(k,j) returns: partial derivative with respect to p_j of (model(xData_k,p)-fData_k)

%% Required files:
% <none>

%% Test cases:
% [myError,myJacobian]=leastSquaresObjective(@simpleValleyModel,[2;3],[0, 0, 1, 2; 1, 2, 3, 4],[2, 3, 2.54, 4.76]);
% should return
% myError close to [2; 3; 10; 20], myJacobian = [0, 1; 1, 1; 4, 1; 9, 1];

% [myError,myJacobian]=leastSquaresObjective(@exponentialModel,[1;-1;0],[0, 1, 2, 3],[0,0,0,0]);
% should return
% myError = [0;0;0;0] and myJacobian =[1,1,0;1,1,-1; 1,1,-2;1,1,-3]

% [myError,myJacobian]=leastSquaresObjective(@benchmarkModel,[2;1;3],[3, 0, 0;-1, -1, 0;-1, -1, -1],[3, 2, 1]);
% should return
% myError = [4;2;2] and myJacobian =[0,-1,2;0,-1,1; 0, 0,1]

%% Input verification:
try
    [value, gradient_x,gradient_p] = model(xData(:,1),p);
catch
    error('evaluation of model handle failed!');
end

[~,N]=size(xData);
if (~isequal(size(fData), [1,N]))
    error('dimensions of xData and fData do not match!');
end

%% Implementation:
% Hints: 
% 1. Do a loop over N
% 2. Each row of the Jacobian is a (transposed) gradient with respect to p of the model


for k = 1:N

  errorVector(k,:) = model(xData(:,k),p) - fData(:,k);
  [value, gradient_x ,gradient_p] = model(xData(:,k),p);
  jacobian(k,:) = gradient_p';
  
endfor
end
