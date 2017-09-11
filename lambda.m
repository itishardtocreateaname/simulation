function lambdas = lambda(x, sr,angle, pl, alpha_t, mag, beta1_t, beta2_t, neuron)
%% lambdas: 1*C, x: 1*2, sr:1*1, pl: 1*C, alpha_t,beta_t:1*C,mag: 1*1
%% sr: real value.

C = size(alpha_t,2);
theta = real2angle(sr, angle);

beta1_t(:,neuron) = mag*cos(theta);
beta2_t(:,neuron) = mag*sin(theta);

k = alpha_t + repmat(x(:, 1), 1, C) .* beta1_t + ...
	repmat(x(:, 2), 1, C) .* beta2_t;
	
lambdas = exp(k) .* pl;