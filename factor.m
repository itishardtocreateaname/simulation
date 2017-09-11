function y_factor = factor(y)

global FACTORS
if y==0
    y_factor = 1;
elseif FACTORS(y)
    y_factor = FACTORS(y);
elseif y==1
    y_factor = 1;
    FACTORS(y) = 1;
else
    y_factor = y*factor(y-1);
    FACTORS(y) = y_factor;
end
    