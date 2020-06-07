function [info,X] = TVexp(dataset,b,noise,maxIter,lambda_reg,initMaxProxIter,sigma,zeta,xi_,initialL)
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
W = data;

W    = (imfilter(W, ones(b)/b^2) + noise) ;
figure();
subplot(1,2,1)
imshow(W,[])
[n,p] = size(W);


%% Initialize
X = zeros(n,p);
Y = zeros(n,p);
Z = zeros(n,p);
L = initialL;


%% Define Objective Function and Prox
g = @(X)(1/2)*norm(W-imfilter(X, ones(b)/b^2),'fro')^2;
gp = @(X)imfilter((imfilter(X, ones(b)/b^2)-W), ones(b)/b^2);
h = @(X)lambda_reg*sum(sum( sqrt(sum( grad(X).^2,3 )) ));
prox = @(X,g,L,proxIter,xik)TV_prox(X,g,L,lambda_reg,proxIter,sigma,zeta,xik);

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
    [X,V,proxIters,proxGap] = prox(Y - lambda_k*gpY,gpY,L,maxProxIter,xi(k));
    gX = g(X);
    f = gX+h(X);
    fprintf('k = %d, proxIters = %d, proxGap = %.10f, f = %.10f\n',k,proxIters,proxGap,f);
    info(end+1,:) = [f proxIters proxGap];
    gpX = gp(X);
    gY = g(Y);
    if backtrack == 1
        % Backtrack if Lipschitz inequality not satisfied
        while gY < gX + gpX(:)'*(Y(:)-X(:)) + 1/2/L*sum((gpX(:)-gpY(:)).^2)
            fprintf('Lipschitz estimate too small, increasing...\n');
            L = L*2;
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


