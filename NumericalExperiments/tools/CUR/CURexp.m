function [info] = CURexp(dataset,maxIter,lambda_row,lambda_col,initMaxProxIter,sigma,zeta,xi_,initialL)
%####################
%
% Mainly based on Mark Schmidt's code
% https://www.cs.ubc.ca/~schmidtm/Software/inexactProx.zip
%
%####################


switch xi_
    case "zero"
        xi = @(k) 0;
    case "poly1"    
        xi = @(k) 1/k;
    case "poly2"    
        xi = @(k) 1/k^2;
    case "poly3"    
        xi = @(k) 1/k^3;
    case "poly4"    
        xi = @(k) 1/k^4;
end        
%xi = @(k) 1/(k);

optTol = 1e-8;
backtrack = 1;



%% Load data
load(dataset);
data = data(:,1:end);

W    = data - mean( data(:) );
W    = W / norm(W,'fro');

[n,p] = size(W);

%% Initialize
X = zeros(p,n);
Z = zeros(p,n);
L = initialL;
lambda_col = lambda_col*sqrt(p/n); % Roughly put rows and columns on same scale

%% Define Objective Function and Prox
g = @(X)(1/2)*norm(W-W*X*W,'fro')^2;
gp = @(X)-W'*(W-W*X*W)*W';
h = @(X)lambda_row*sum(sqrt(sum(X'.^2))) + lambda_col*sum(sqrt(sum(X.^2)));
prox = @(X,g,L,proxIter,xik)CUR_prox(X,g,L,lambda_row,lambda_col,proxIter,sigma,zeta,xik);

%% Run proximal-gradient method
gX = g(X);
f = gX+h(X);
gY = gX;
fY = f;  
info = [f 0 0];
Ak = 0;
tk = 1;

k = 1;

maxProxIter = initMaxProxIter;
while 1
	X_old = X;
    
        
    lambda_k = (1-sigma^2)/L;
    eta_k = (1-zeta^2)*lambda_k;
    ak = (eta_k +sqrt(4*eta_k*Ak+eta_k^2))/2;

    Y = X_old + ak/(ak +Ak)*(Z-X_old); 

    % Compute proximal gradient step
    gpY = gp(Y);
    gY = g(Y);
    [X,V,proxIters,proxGap] = prox(Y - lambda_k*gpY,gpY,L,maxProxIter,xi(k));
    gX = g(X);
    f = gX+h(X);
    fprintf('k = %d, proxIters = %d, proxGap = %.10f, f = %.10f\n',k,proxIters,proxGap,f);
    info(end+1,:) = [f proxIters proxGap];
    gpX = gp(X);
    if backtrack == 1
        % Backtrack if Lipschitz inequality not satisfied
        while gY  < gX + gpX(:)'*(Y(:)-X(:)) + 1/2/L*sum((gpX(:)-gpY(:)).^2) 
            fprintf('Lipschitz estimate too small, increasing...\n');
            L = L*2;
            lambda_k = (1-sigma^2)/L;
            eta_k = (1-zeta^2)*lambda_k;
            ak = (eta_k +sqrt(4*eta_k*Ak+eta_k^2))/2;

            Y = X_old + ak/(ak +Ak)*(Z-X_old); 

            % Compute proximal gradient step
            gpY = gp(Y);
            [X,V,proxIters,proxGap] = prox(Y - lambda_k*gpY,gpY,L,maxProxIter,xi(k));
            gX = g(X);
            gY = g(Y);
            f = gX+h(X);
            fprintf('k = %d, proxIters = %d, proxGap = %.10f, f = %.10f, L = %f\n',k,proxIters,proxGap,f,L);
            info(end+1,:) = [f proxIters proxGap];
            gpX = gp(X);
            if sum(info(:,2)) >= maxIter
                return;
            end
        end
    end




    % Compute Z

 
    Z = Z - ak*(V+gpY);


    Ak = ak+Ak;
        
           
    
	if sum(info(:,2)) >= maxIter
		break;
	end
	k = k + 1;
    L = L/1.1;

    
end


