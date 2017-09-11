function gradients = gradient(sr, angle, y, x, pl, alpha_t, mag, beta1_t, beta2_t, neuron) 
%% sr: 1*1, y: 1*C, x: 1*2, pl: 1*C, alpha_t,beta1_t:1*C, mag: 1*1
%% y= 0 or 1, sr: real value.

C = size(y,2);
theta = real2angle(sr, angle);

lambdas = lambda(x, sr,angle, pl, alpha_t, mag, beta1_t, beta2_t, neuron);
lambda_ders = lambda_der(x, sr,angle, pl, alpha_t, mag, beta1_t, beta2_t, neuron);
theta_ders = theta_der(sr, angle);

tmp = (lambdas >= 1);
lambdas(tmp) = 1-eps(1);

%% derivative of negative log-likehood loss function (binary distribution)
global period
gradients = (-y./lambdas+period).*lambda_ders.*theta_ders;

%  gradients = -y.*lambda_ders.*theta_ders./lambdas + (1-y).*lambda_ders.*theta_ders./(1-lambdas);

gradients = gradients(:,neuron);