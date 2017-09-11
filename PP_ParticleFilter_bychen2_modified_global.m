function [X,PD] = PP_ParticleFilter_bychen2_modified_global(Y,mu0,cov0,F,R,alpha,w0,N,p,neuron,top,bot,pl,mag)
%% Y: T*C, mu0: 1*2, cov0:2*2, alpha:1*1, w0:1*C, neuron:1*1, pl:1*C, mag:1*1

T = size(Y,1);
C = size(Y,2);
X = zeros(T,2); %% estimated x
PD = zeros(T,1); %% estimated sd

sw = 1/N*ones(1,N); %% sw: w_t^i (1*N)
angle = [0 pi];
sx = mvnrnd(mu0,cov0,N);%% sx: x_t^i (N*2)
sd = w0(neuron)+1e-4*randn(1,1); %% sd: sd_t (1*1)
sr = angle2real(sd,angle); %% sr: sd_t in real number.
%%sd = repmat(w0(neuron),N,1)+1e-4*randn(N,1); %% sd: sd_t^i (N*1)
% sr = sd;


PD(1) = sd; %% sd_0
%%PD(1) = sw*sd;
X(1,:) = sw*sx; %% X_0

%######################################
%{
T0 = 200; %% proper value?????????????
phi1 = 0.95;
phi2 = 0.9975;
eta = 0.008;
%}
global TRIAL
global T0
global PHI1
global PHI2
global ETA
global GLOBAL
global MOD

t0 = T0;
phi1 = PHI1;
phi2 = PHI2;
eta = ETA;
global_sampling = GLOBAL;
%######################################

m_t = zeros(1, 1);
v_t = zeros(1, 1);

alpha_t = alpha*ones(1,C);
beta1_t = mag*cos(w0); %%1*C
beta2_t = mag*sin(w0);

for t = 2:T
        
    %% update sd&sr using ADAM %%%%%%%%%%%%%%%%%%%%%%%%%
    pred_x = X(t-1,:)*F;
	[sr, m_t, v_t] = Adam(sr, angle, m_t, v_t,t, Y(t,:), pred_x, phi1, phi2, eta, pl, alpha_t, mag, beta1_t, beta2_t, neuron);
    sd = real2angle(sr, angle);    
	
    if sum(isnan(sd))>0
        disp('pause');
    end
    
	%% global sampling%%%%%%%%%%%%%%%%%%%%%%%
	if global_sampling && mod(t,t0)==mod(MOD,t0) && t>=T0
	
		sx2 = sx*F; % pred_x; N*2

		rnd_theta = unifrnd(bot, top, N, 1); %N*1

		alpha1 = repmat(alpha_t, N, 1); %N*C
		beta1 = repmat(beta1_t, N, 1);
		beta2 = repmat(beta2_t, N, 1);
		
		beta1(:, neuron) = mag*cos(rnd_theta);
		beta2(:, neuron) = mag*sin(rnd_theta);

		k = alpha1+repmat(sx2(:,1),1,C).*beta1+...  
			repmat(sx2(:,2),1,C).*beta2;			%N*C
			
		lambda2 = exp(k).*repmat(pl,N,1); %% lambda's samples; N*C
		tmp = (lambda2>=1);
		lambda2(tmp) = 1-eps(1);
		
		weights2 = PP_likelihood_bychen(Y(t,:),lambda2,'bino');
		%weights2 = weights2/sum(weights2);
		
		
        sd_adam = sd; 
		sd = repmat(sd, N, 1); %duplicate sd to be N samples.
		beta1(:, neuron) = mag*cos(sd);
		beta2(:, neuron) = mag*sin(sd);

		k = alpha1+repmat(sx2(:,1),1,C).*beta1+...  
			repmat(sx2(:,2),1,C).*beta2;			%N*C
			
		lambda1 = exp(k).*repmat(pl,N,1); %% lambda's samples; N*C
		tmp = (lambda1>=1);
		lambda1(tmp) = 1-eps(1);
		
		weights1 = PP_likelihood_bychen(Y(t,:),lambda1,'bino');
		%weights1 = weights1/sum(weights1);
		
		w0 = weights1.*sw; %1*N
		wA = weights2.*sw;
		
		weights3 = [(1-p)*w0 p*wA];
		weights3 = weights3/sum(weights3);
		
		k = resample2(weights3);
		k1 = k(k<=N);
		k2 = k(k>N);
		k2 = k2-N;
		k = [k1 k2];
		wX = [weights1(k1) weights2(k2)];
		
		%%update sd
		sd = [sd(k1); rnd_theta(k2)];
		%sr = angle2real(sd, angle);
		if sum(isnan(sd))>0
			disp('pause');
		end
		
		%%propagate x
		sx = sx(k, :);
		sx2 = sx*F + mvnrnd([0 0],R,N);
		beta1(:,neuron) = mag*cos(sd);
		beta2(:,neuron) = mag*sin(sd);
		k = alpha1 + repmat(sx2(:,1),1,C).*beta1 + ...
			repmat(sx2(:,2),1,C).*beta2;
			
		lambda3 = repmat(pl,N,1).*exp(k);
		tmp = (lambda3>=1);
		lambda3(tmp) = 1-eps(1);
		
	%   sw = PP_likelihood_bychen(Y(t,:),lambda3,'poisson');
		sw = PP_likelihood_bychen(Y(t,:),lambda3,'bino');
		
		%sw = exp(sw-max(sw))./wX;
		sw = sw./wX;
		sw = sw/sum(sw);
		
		sd = sw*sd;
        g_t = sd_adam - sd;
        if g_t
            O = round(log(abs(m_t/g_t)));
        else
            O = -1;
        end
        g_t = g_t * exp(O+2);  
        m_t = phi1*m_t + (1-phi1)*g_t;
        v_t = phi2*v_t + (1-phi2)*g_t.*g_t;
        
		sr = angle2real(sd,angle);
		PD(t) = sd;
		X(t,:) = sw*sx2;
		
		kk = resample(sw);
		%SD(:,t) = sd(kk);
		sx = sx2(kk,:);
		%sd = sd(kk,:);
		sw = 1/N*ones(1,N);
		
		%sr = angle2real(sd,angle);
	%     sr = sd;
		%if sum(isinf(sr))>0
		%	disp('pause');
		%end
		        
    else
	     
		%% reweight
		sx2 = sx*F; 
		lambda2 = lambda(sx2, repmat(sr,N,1),angle, repmat(pl,N,1), repmat(alpha_t,N,1), mag, repmat(beta1_t,N,1), repmat(beta2_t,N,1), neuron);
		tmp = (lambda2>=1);
		lambda2(tmp) = 1-eps(1);
		weights2 = PP_likelihood_bychen(Y(t,:),lambda2,'bino');
		
		%% propagate x
		k = resample(weights2./sum(weights2));
		sx = sx(k,:); 
		sx2 = sx*F + mvnrnd([0 0],R,N);
		lambda3 = lambda(sx2, repmat(sr,N,1),angle, repmat(pl,N,1), repmat(alpha_t,N,1), mag, repmat(beta1_t,N,1), repmat(beta2_t,N,1), neuron);
		tmp = (lambda3>=1);
		lambda3(tmp) = 1-eps(1);
		weights3 = PP_likelihood_bychen(Y(t,:),lambda3,'bino');

		%% calculate weights for new samples
		sw = weights3./weights2;
		sw = sw/sum(sw);

		%% store estimated value
		PD(t) = sd;
		X(t,:) = sw*sx2;
	 
		%% resample to get equal weights
		kk = resample(sw);
		sx = sx2(kk,:);
		sw = 1/N*ones(1,N);	
		
	end
end
