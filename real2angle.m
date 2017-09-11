function sd = real2angle(sr,angle)

t = exp(sr);
p = t./(1+t);
sd = p*(angle(2)-angle(1))+angle(1);