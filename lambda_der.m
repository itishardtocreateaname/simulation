function derivatives = lambda_der(x, sr,angle, pl, alpha_t, mag, beta1_t, beta2_t, neuron)
%% lambdas: 1*C, x: 1*2, theta:1*1,pl: 1*C, alpha_t,beta1_t:1*C,mag: 1*1
%% theta: angel.
%% derivative of lambda w.r.t theta
N = size(alpha_t,1);
C = size(alpha_t,2);
theta = real2angle(sr, angle);

beta1_t(:,neuron) = mag*cos(theta);
beta2_t(:,neuron) = mag*sin(theta);

k = alpha_t + repmat(x(:, 1), 1, C) .* beta1_t + ...
	repmat(x(:, 2), 1, C) .* beta2_t;
	
beta1_der = zeros(N,C);
beta2_der = zeros(N,C);
	
beta1_der(:,neuron) = mag*(-sin(theta));
beta2_der(:,neuron) = mag*(cos(theta));	
	
derivatives = pl .* exp(k) .* (repmat(x(:, 1), 1, C) .* beta1_der + ...
								repmat(x(:, 2), 1, C) .* beta2_der);