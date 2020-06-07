function [W,V,i,gap] = TV_prox(W0,gpY,L,lambda_reg,maxIter,sigma,zeta,epsilon)
lambda = (1-sigma^2)/L;
p_old = zeros([size(W0) 2]);
u = zeros([size(W0) 2]); 
t = 1;
for i = 1:maxIter
	% Prox with respect to columns
	p = projection(u - 1/8*grad(-div(u(:,:,1),u(:,:,2))-W0),lambda,lambda_reg);
    t_ = (1+sqrt(1+4*t^2))/2;
    u = p + (t-1)/t_*(p-p_old);
    t = t_;
    p_old = p;
    
    V =  -div(p(:,:,1),p(:,:,2))/lambda;
    W = W0 -  lambda*V;
	f_primal = (1/2)*norm(W-W0,'fro')^2 + lambda*lambda_reg*sum(sum( sqrt(sum( grad(W).^2,3 )) ));
	f_dual = (1/2)*norm(W0,'fro')^2 - (1/2)*norm(W0 - V*lambda ,'fro')^2;
	gap = f_primal - f_dual;
    
	if  gap < sigma^2/2*norm(W-(W0+lambda*gpY),'fro')^2+zeta^2/2*lambda^2*norm(V+gpY,'fro')^2+epsilon/2
		break;
	end
end
end

%% Prox solvers with respect to blocks
function [p] = projection(p,lambda,lambda_reg)
p = p./max(1,sqrt(sum(p.^2,3))/lambda/lambda_reg);
end

