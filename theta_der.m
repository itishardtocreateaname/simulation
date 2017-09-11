function der = theta_der(sr, angle)
%% derivative of sd w.r.t sr

t = exp(sr);
der = (angle(2)-angle(1))*t./((1+t).*(1+t));