clear all; clc;
% In this example, we use the optimized relatively inexact proximal point
% algorithm (ORI-PPA) introduced in [1] for solving the non-smooth convex
% minimization problem:
%
% min_x F(x); for notational convenience we denote xs in argmin_x F(x);
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying ||x0-xs||<=1.
%
% For a given sequence of step sizes lambda_k,
% the (ORI-PPA) algorithm defines iterates as follows
%   A_{k+1} = A_k + lambda_k/2 + sqrt(4*lambda_k*A_k+lambda_k^2)/2;
%   y_k     = x_k + lambda_k / (A_{k+1} - A_k) (z_k - x_k);
%   x_{k+1} = y_k - lambda_{k} * (g_{k+1} + e_{k+1})
% with PD_{lambda_k F}(x_{k+1},g_{k+1};y_k) <= sigma^2/2 * (x_{k+1}-y_k)^2;
%   z_{k+1} = z_k + 2/(1+sigma)(A_{k+1} - A_k) g_{k+1};
%
% This method is presented in [1, Section 4.3].

% [1] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
F = P.DeclareFunction('Convex');

% (2) Set up the starting point and initial condition
x0=P.StartingPoint();		 % x0 is some starting point
[xs,Fs]=F.OptimalPoint(); 	 % xs is an optimal point, and Fs=F(xs)
P.InitialCondition((xs-x0)^2<=1)

% (3) Set up the method
N               = 5;
sigma           = sqrt(.5);
opt.criterion   = 'PD_gapI';

x           = cell(N+1,1);
y           = cell(N,1);
z           = cell(N+1,1);
Fx          = cell(N+1,1);
lambda      = ones(1,N);
% lambda      = rand(1,N);
Ak          = zeros(1,N+1);
Ak(1)       = 0;
x{1}        = x0;
z{1}        = x0;

for i = 1:N
    Ak(i+1) = Ak(i) + lambda(i)/2 + sqrt(4*lambda(i)*Ak(i)+lambda(i)^2)/2;
    y{i} = x{i} + lambda(i)/(Ak(i+1)-Ak(i))*(z{i}-x{i});
    [x{i+1},~,Fx{i+1},~,g,~,epsVar] = inexact_proximal_step(y{i},F,lambda(i),opt);
    P.AddConstraint( epsVar <= sigma^2/2 * (y{i}-x{i+1})^2);
    z{i+1} = z{i} - 2*(Ak(i+1)-Ak(i))/(1+sigma)*g;
end


% (4) Set up the objective
P.PerformanceMetric(Fx{N+1}-Fs); % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve(2)

% (6) Evaluate the output
[double(Fx{N+1}-Fs), (1+sigma)/4/Ak(N+1)]
