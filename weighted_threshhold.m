function [Z1] = weighted_threshhold(X1, lambda)
    Power = abs(X1);    
    A = Power - lambda;
    A(A < 0) = 0;
    Z1 = A .* sign(X1);
end