function y_factors = factors(y)

[a,b]=size(y);
y_factors = zeros(a,b);
for i = 1:a
    for j = 1:b
        y_factors(i,j) = factor(y(i,j));
    end
end