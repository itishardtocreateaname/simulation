function [X,PD] = PP_ParticleFilter_bychen2_modified(freq_Y, testY,mu0,cov0,F,R,alpha,w0,N,p,neuron,top,bot,pl,mag)
%% freq_Y: T*C, mu0: 1*2, cov0:2*2, alpha:1*1, w0:1*C, neuron:1*1, pl:1*C, mag:1*1

T = size(freq_Y,1);
C = size(freq_Y,2);
X = zeros(T,2); %% estimated x
PD = zeros(T,1); %% estimated theta

sw = 1/N*ones(1,N); %% sw: w_t^i (1*N)
angle = [0 pi];
sx = mvnrnd(mu0,cov0,N);%% sx: x_t^i (N*2)
sd = w0(neuron)+1e-4*randn(1,1); %% sd: theta_t (1*1)
sr = angle2real(sd,angle); %% sr: theta_t in real number.
%%sd = repmat(w0(neuron),N,1)+1e-4*randn(N,1); %% sd: theta_t^i (N*1)
% sr = sd;


PD(1) = sd; %% theta_0
%%PD(1) = sw*sd;
X(1,:) = sw*sx; %% X_0

%######################################
%{
T0 = 200; %% proper value?????????????
phi1 = 0.95;
phi2 = 0.9975;
eta = 0.008;
%}
global T0
global PHI1
global PHI2
global ETA
global GLOBAL

t0 = T0;
phi1 = PHI1;
phi2 = PHI2;
eta = ETA;
global_sampling = GLOBAL;
%######################################

m_t = zeros(1, 1);
v_t = zeros(1, 1);
gt = zeros(1,T);
mt = zeros(1,T);

alpha_t = alpha*ones(1,C);
beta1_t = mag*cos(w0); %%1*C
beta2_t = mag*sin(w0);

adam_total_time = 0;
likelihood_total_time = 0;

for t = 2:T
%     if global_sampling
%         if mod(t, t0) == 0
%             sx2 = sx*F;
% 
%             theta = unifrnd(bot, top, N, 1);
% 
%             alpha1 = repmat(alpha_t, N, 1);
%             beta1 = repmat(beta1_t, N, 1);
%             beta2 = repmat(beta2_t, N, 1);
%             beta1(:, neuron) = mag*cos(theta);
%             beta2(:, neuron) = mag*sin(theta);
% 
%             k = alpha1+repmat(sx2(:,1),1,C).*beta1+...
%                 repmat(sx2(:,2),1,C).*beta2;
%             lambda1 = exp(k).*repmat(pl,N,1);
%             tmp = (lambda1>=1);
%             lambda1(tmp) = 1-eps(1);
%             weights1 = PP_likelihood_bychen(freq_Y(t,:),lambda1,'bino');
%             weights1 = weights1/sum(weights1);
% 
%             theta_hat = weights1*theta;
%             real_theta_hat = angle2real(theta_hat, angle);
% 
%             m_t = sr-real_theta_hat;
%             v_t = (sr-real_theta_hat).*(sr-real_theta_hat);
% 
%         end
%     end
% 	
    
    %% update sd&sr using ADAM
    %
    pred_x = X(t-1,:)*F;
    
    tic;
	[sr, m_t, v_t, g_t] = Adam(sr, angle, m_t, v_t,t, freq_Y(t,:), pred_x, phi1, phi2, eta, pl, alpha_t, mag, beta1_t, beta2_t, neuron);
    adam_time = toc;
    adam_total_time = adam_total_time+adam_time;
    
    gt(t)=g_t;
    mt(t)=m_t;
    
    sd = real2angle(sr, angle);    
	
    if sum(isnan(sd))>0
        disp('pause');
    end
    %}
    %% another update way using ADAM
    %{
    sr_hat=repmat(sr,10,1);
    m_t_hat=repmat(m_t,10,1);
    v_t_hat=repmat(v_t,10,1);
    for i=1:10
        [sr_hat(i,:),m_t_hat(i,:),v_t_hat(i,:)]=Adam(sr,angle,m_t,v_t,t,freq_Y(t,:),mean(sx((i-1)*(N/10)+1:i*(N/10),:))*F,phi1,phi2,eta,pl,alpha_t,mag,beta1_t,beta2_t,neuron);
    end
    
    sr=mean(sr_hat);
    sd=real2angle(sr,angle);
    m_t=mean(m_t_hat);
    v_t=mean(v_t_hat);
    
    if sum(isnan(sd))>0
        disp('pause');
    end
    %} 

    %% reweight
    sx2 = sx*F; 
    lambda2 = lambda(sx2, repmat(sr,N,1),angle, repmat(pl,N,1), repmat(alpha_t,N,1), mag, repmat(beta1_t,N,1), repmat(beta2_t,N,1), neuron);
    tmp = (lambda2>=1);
    lambda2(tmp) = 1-eps(1);
%     weights2 = PP_likelihood_bychen(Y(t,:),lambda2,'bino');
    
    tic;
    weights2 = PP_likelihood_bychen(testY(t,:),lambda2,'poisson');
%     weights2 = PP_likelihood_bychen(freq_Y(t,:),lambda2,'poisson');
    likelihood_time = toc;
    likelihood_total_time = likelihood_total_time+likelihood_time*2;
    
    %% propagate x
    k = resample(weights2./sum(weights2));
    sx = sx(k,:); 
    sx2 = sx*F + mvnrnd([0 0],R,N);
    lambda3 = lambda(sx2, repmat(sr,N,1),angle, repmat(pl,N,1), repmat(alpha_t,N,1), mag, repmat(beta1_t,N,1), repmat(beta2_t,N,1), neuron);
    tmp = (lambda3>=1);
    lambda3(tmp) = 1-eps(1);
%     weights3 = PP_likelihood_bychen(Y(t,:),lambda3,'bino');
    weights3 = PP_likelihood_bychen(testY(t,:),lambda3,'poisson');
%     weights3 = PP_likelihood_bychen(freq_Y(t,:),lambda3,'poisson');

    sw = weights3./weights2;
    sw = sw/sum(sw);

    PD(t) = sd;
    X(t,:) = sw*sx2;

    kk = resample(sw);
    sx = sx2(kk,:);
    sw = 1/N*ones(1,N);

end
adam_total_time
likelihood_total_time
% plot(gt);hold on;
% plot(mt);