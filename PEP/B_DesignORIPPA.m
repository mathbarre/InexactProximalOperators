clear all; clc;
% Design of optimized relatively inexact proximal point algorithm
% (ORI-PPA) introduced in [1].
% Detail on the design process presented in [1, Section 4.2].
% Solving an approximate version of problem (14), where the inner maximization is
% dualized and only a subset convex interpolation inqualities are used.
% x_k is replaced by its expression in (12) and bilinear terms are removed 
% by introducing eta and gamma variables. 
%
% Require Yalmip and Mosek (or any SDP solver but then change the solver 
% options )
%
%
% ORI-PPPA: optimized relatively inexact proximal point algorithm (see [1])
% [1] M. Barre, A. Taylor, F. Bach. Principled analyses and design of
%     first-order methods with inexact proximal operators


clc; clear all;

verbose   = 1;
tolerance = 1e-8;
N = 5;
lambda = rand(N,1);
sigma = 0.4;



dimG = 3*N+1;%x0,e1..en,v1..vn,u1..un
dimF = 2*N;%hx1..hxn,hu1..hun
x0   = zeros(1,dimG); x0(:,1) = 1; 
xs   = zeros(1,dimG);
u    = zeros(N,dimG); u(:,2:(N+1)) = eye(N);
v    = zeros(N,dimG); v(:,(N+2):(2*N+1)) = eye(N);
e    = zeros(N,dimG); e(:,(2*N+2):(3*N+1)) = eye(N);


hx  = zeros(N, dimF); hx(:,1:N) = eye(N);
hu  = zeros(N, dimF); hu(:,(N+1):(2*N)) = eye(N);
hs  = zeros(1, dimF);



tau = sdpvar(1); % multiplier associated with initial condition 
nu_s = sdpvar(N,1); % multiplier associated with interpolation between xs and u_k
nu = sdpvar(N-1,1); % multiplier associated with interpolation between x_k and u_{k+1}
nu_ = sdpvar(N,1); % multiplier associated with interpolation between x_k and u_k
mu = sdpvar(N,1); % multiplier associated with inexactness criterion

hxN = hx(N,:);


cons_SDP = -1/4*(x0 - xs).'*(x0 - xs);

cons_LIN = tau*(hxN - hs) - nu_s(N)*(hu(N,:) - hs) - mu(N)*(hx(N,:) - hu(N,:)) - nu_(N)*(hu(N,:) - hx(N,:));


for i = 1:(N-1)
    cons_LIN = cons_LIN - nu_s(i)*(hu(i,:) - hs)...
               - nu(i)*(hu(i+1,:) - hx(i,:))...
               - mu(i)*(hx(i,:) - hu(i,:))...
               - nu_(i)*(hu(i,:) - hx(i,:));
    cons_SDP = cons_SDP -nu_s(i)*v(i,:).'*(xs - u(i,:))...
               - nu(i)*(v(i+1,:).'*(x0 - lambda(i)*(v(i,:) + e(i,:)) - u(i+1,:)))...
               - mu(i)*(lambda(i)/2*e(i,:).'*e(i,:) - lambda(i)*sigma^2/2*(v(i,:) + e(i,:)).'*(v(i,:) + e(i,:))...
                        - v(i,:).'*(x0 - lambda(i)*(v(i,:) + e(i,:)) - u(i,:)))...
               - nu_(i)*(v(i,:).'*(x0 - lambda(i)*(v(i,:) + e(i,:)) - u(i,:))); 
end

cons_SDP = cons_SDP - nu_s(N)*v(N,:).'*(xs - u(N,:))...
           - mu(N)*(lambda(N)/2*e(N,:).'*e(N,:) - lambda(N)*sigma^2/2*(v(N,:) + e(N,:)).'*(v(N,:) + e(N,:))...
                    - v(N,:).'*(x0-lambda(N)*(v(N,:)+e(N,:))-u(N,:)))...
           - nu_(N)*(v(N,:).'*(x0 - lambda(N)*(v(N,:) + e(N,:)) - u(N,:))); 

    
       
% Bilinear terms in the multiplicator and method parameters are replaced by
% eta dans gamma variables.

eta = sdpvar(N-1,N-1);
gamma = sdpvar(N-1,N-1);

for i = 1:(N-1)
    for j = 1:i
        cons_SDP = cons_SDP + eta(i,j)*(v(i+1,:).'*v(j,:)) + gamma(i,j)*(v(i+1,:).'*e(j,:));    
    end
end


cons_SDP = (cons_SDP.' + cons_SDP)/2;
cons = cons_SDP <= 0;
cons = cons+ (cons_LIN == 0);
cons = cons + (nu_s >=0);
cons = cons + (nu >= 0);
cons = cons + (nu_ >= 0);
cons = cons + (mu >= 0);





obj = -tau;


solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tolerance);
solverDetails = optimize(cons,obj,solver_opt);

Ak = zeros(N+1,1);
for i = 1:N
    Ak(i+1) = Ak(i) + (lambda(i) + sqrt(4*lambda(i)*Ak(i) + lambda(i)^2))/2;
end


% The rate obtained is equal to that of Theorem 4.1
[double(tau), Ak(end)/(1+sigma)]


% Recovers parameters alpha and beta in the expression of x_k (11).

alpha = zeros(N,N);
beta = zeros(N,N);
for i = 1:(N-1)
    alpha(i+1,i) = -double(eta(i,i)/(mu(i+1) - nu_(i+1)));
    beta(i+1,i) = -double(gamma(i,i)/(mu(i+1) - nu_(i+1)));
    for j = (i-1):(-1):1
       alpha(i+1,j) = double((nu(i)*alpha(i,j) - eta(i,j))/((mu(i+1) - nu_(i+1)))); 
       beta(i+1,j) = double((nu(i)*beta(i,j) - gamma(i,j))/(mu(i+1) - nu_(i+1)));      
    end
end
 
