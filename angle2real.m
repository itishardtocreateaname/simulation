function sr = angle2real(sd,angle)

p = (sd-angle(1))/(angle(2)-angle(1));
sr = log(p./(1-p));