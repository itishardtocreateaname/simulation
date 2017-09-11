function weights = PP_likelihood_bychen(y,lambda,sym)

C = size(y,2);
N = length(lambda);
tmp = (lambda>=1);
lambda(tmp) = 1-eps(1);

y = repmat(y,N,1);
if C==size(lambda,2)
    r = lambda;
elseif C>size(lambda,2)
    r = repmat(lambda,1,C);
end

switch sym
    case 'poisson'
        global period
%         weights = poisspdf(y,r);
         weights = poisspdf(y,r*period);

        %tmp = sum(weights,1); idx = tmp>0;
        %tmp = repmat(tmp(idx),N,1);
        %weights = weights(:,idx)./tmp;
        %weights = prod(weights*N,2);
        weights = prod(weights,2);
        %weights = log(weights');
        weights = weights';
    case 'bino'
        tmp1 = (y==1);
        tmp0 = (y==0);
        weights = zeros(size(lambda));
        weights(tmp1) = r(tmp1);
        weights(tmp0) = 1-r(tmp0);
        weights = prod(weights,2);
        weights = weights';
end

