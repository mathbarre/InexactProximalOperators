clear all; clc;
% In this example, we use a relatively inexact proximal point algorithm
% for solving the non-smooth mu-strongly convex minimization problem
%   min_x F(x); for notational convenience we denote xs = argmin_x F(x);
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying F(x0)-F(xs)<=1.
%
% The inexact proximal point algorithm produces a x_{k+1} as
%   x_{k+1} = x_k - gamma * ( F'(x_{k+1}) + e )
% with ||e||^2 <= sigma^2/gamma^2 * (x_{k+1}-x_k)^2.
% 
%
% The example is used in [1, Section 5].

% [1] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators


% (0) Initialize an empty PEP
P = pep();
params.L = Inf;
params.mu = 0.1;
% (1) Set up the objective function
F = P.DeclareFunction('SmoothStronglyConvex',params);

% (2) Set up the starting point and initial condition
x0 = P.StartingPoint();		 % x0 is some starting point
F0 = F.value(x0);
[xs,Fs] = F.OptimalPoint(); 	 % xs is an optimal point, and fs=F(xs)
P.InitialCondition(F0 - Fs <= 1)

% (3) Set up the method
N               = 5;
sigma           = sqrt(.5);
opt.criterion   = 'PD_gapII';

x           = cell(N+1,1);
Fx          = cell(N+1,1);
x{1}        = x0;
lambda      = ones(1,N);
% lambda      = rand(1,N);
rho         = zeros(1,N);

for i = 1:N
    [x{i+1},~,Fx{i+1},~,~,~,epsVar] = inexact_proximal_step(x{i},F,lambda(i),opt);
    P.AddConstraint( epsVar <= (sigma/lambda(i))^2 * (x{i}-x{i+1})^2);    
    rho(i) = (1 + sigma) / (1 + sigma + lambda(i)*params.mu);
end


% (4) Set up the objective
P.PerformanceMetric(Fx{N+1}-Fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve(2)

% (6) Evaluate the output
[double(Fx{N+1}-Fs), prod(rho.^2)]
