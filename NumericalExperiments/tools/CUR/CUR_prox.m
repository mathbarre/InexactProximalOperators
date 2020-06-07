function [W,V,i,gap] = CUR_prox(W0,gpY,L,lambda_row,lambda_col,maxIter,sigma,zeta,epsilon)
%####################
%
% Mainly based on Mark Schmidt's code
% https://www.cs.ubc.ca/~schmidtm/Software/inexactProx.zip
%
%####################

lambda = (1-sigma^2)/L;
W = W0;
Xi_row = zeros(size(W));
Xi_col = zeros(size(W));
for i = 1:maxIter
	% Prox with respect to rows
    W      = W0 - Xi_col;
	Xi_row = projection_row(W,1/lambda,lambda_row);
    
    
	% Prox with respect to columns
    W      = W0 - Xi_row;
	Xi_col = projection_col(W,1/lambda,lambda_col);
	
    % Dual estimate
    V =  (Xi_row+Xi_col)/lambda;
    % Primal estimate
    W     = W0 - lambda * V;
    
	% Duality gap
	f_primal = (1/2)*norm(W-W0,'fro')^2 + (lambda_row*lambda)*sum(sqrt(sum(W'.^2))) + (lambda_col*lambda)*sum(sqrt(sum(W.^2)));
    f_dual = (1/2)*norm(W0,'fro')^2 - (1/2)*norm(W0 - Xi_row-Xi_col,'fro')^2;
	
    gap = f_primal - f_dual;
    epssss=sigma^2/2*norm(W-(W0+lambda*gpY),'fro')^2+zeta^2/2*lambda^2*norm(V+gpY,'fro')^2+epsilon/2;
	if  gap < epssss
		break;
	end
end
end

%% Prox solvers with respect to blocks
function [W] = projection_row(W,L,lambda)
nRows = size(W,1);
for r = 1:nRows
	nrm = norm(W(r,:));
	if nrm > 0
		W(r,:) = (W(r,:)/nrm)*min(nrm,lambda/L);
	end
end
end

function [W] = projection_col(W,L,lambda)
nCols = size(W,2);
for c = 1:nCols
	nrm = norm(W(:,c));
	if nrm > 0
		W(:,c) = (W(:,c)/nrm)*min(nrm,lambda/L);
	end
end
end








