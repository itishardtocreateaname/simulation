function [X,PD] = PP_ParticleFilter_bychen2(Y,mu0,cov0,F,R,alpha,w0,N,p,neuron,top,bot,pl,mag)

T = size(Y,1);
C = size(Y,2);
X = zeros(T,2);
PD = zeros(T,1);

delta = 0.99;
h2 = 1-((3*delta-1)/(2*delta))^2;
a = sqrt(1-h2);

sw = 1/N*ones(1,N);
angle = [0 pi];
sx = mvnrnd(mu0,cov0,N);%#########cov0应该是什么？
sd = repmat(w0(neuron),N,1)+1e-4*randn(N,1);
sr = angle2real(sd,angle);
% sr = sd;

SD = zeros(N,T);
SD(:,1) = sd;

PD(1) = sw*sd;
X(1,:) = sw*sx;

alpha_t = alpha*ones(N,C);
beta1_t = mag*repmat(cos(w0),N,1);
beta2_t = mag*repmat(sin(w0),N,1);

for t = 2:T
%     if mod(t,2000)==0
%        disp(t); 
%     end
%     
%     K1_t = repmat(K1t(t,:),N,1);
%     K2_t = repmat(K2t(t,:),N,1);
%     K3_t = repmat(K3t(t,:),N,1); 
    
    pr_mea = sw*sr;
    pr_cov = sw*(sr-pr_mea).^2;
    pr_mus = a*sr + (1-a)*pr_mea;
    
    sx2 = sx*F;
    sd = real2angle(pr_mus,angle); 
%     sd = pr_mus;
    if sum(isnan(sd))>0
        disp('pause');
    end
    
%     alpha_t(:,neuron) = repmat(alpha,N,1);
    beta1_t(:,neuron) = mag*cos(sd);
    beta2_t(:,neuron) = mag*sin(sd);

    k = alpha_t + repmat(sx2(:,1),1,C).*beta1_t + ...
        repmat(sx2(:,2),1,C).*beta2_t;
    lambda1 = exp(k).*repmat(pl,N,1);%%%%为什么会有大于1的情况？
    tmp = (lambda1>=1);
    lambda1(tmp) = 1-eps(1);
    %weights1 = PP_likelihood_bychen(Y(t,:),lambda1,'poisson');
     weights1 = PP_likelihood_bychen(Y(t,:),lambda1,'bino');
    
    sd = unifrnd(bot,top,N,1);
    beta1_t(:,neuron) = mag*cos(sd);
    beta2_t(:,neuron) = mag*sin(sd);

    k = alpha_t + repmat(sx2(:,1),1,C).*beta1_t + ...
        repmat(sx2(:,2),1,C).*beta2_t;
    lambda2 = exp(k).*repmat(pl,N,1);
    tmp = (lambda2>=1);
    lambda2(tmp) = 1-eps(1);
   % weights2 = PP_likelihood_bychen(Y(t,:),lambda2,'poisson');
     weights2 = PP_likelihood_bychen(Y(t,:),lambda2,'bino');
    
%     w3 = exp(weights1-max(weights1));
%     w0 = w3.*sw;
%     w0 = w0/sum(w0);
%     
%     w4 =exp(weights2-max(weights2));
%     wA = w4.*sw;
%     wA = wA/sum(wA);

% debug
%     if t > 3000 && Y(t)==1
%         c = 1;
%     end
%     close all
%     figure(1)
%     subplot(2,1,1);plot(weights1);
%     subplot(2,1,2);plot(lambda1);
%     figure(2)
%     subplot(2,1,1);plot(weights2);
%     subplot(2,1,2);plot(lambda2);

    w0 = weights1.*sw;
    wA = weights2.*sw;
    
    weights3 = [(1-p)*w0 p*wA];
    weights3 = weights3/sum(weights3);
    k = resample2(weights3);
    k1 = k(k<=N);
    k2 = k(k>N);
    k2 = k2-N;
    if length(k2)>10
        c = 1;
    end
    k = [k1 k2];
    wX = [weights1(k1) weights2(k2)];
    
    sx = sx(k,:);%#######这里不应该是sx = sx2(k,:);？？？？
    if pr_cov<=0||isnan(pr_cov)
        pr_cov = 1e-8;
    end
    sr = pr_mus(k1)+randn(length(k1),1)*sqrt(h2*pr_cov);
    %sr = pr_mus(k1)+randn(length(k1),1)*chol(h2*pr_cov);
    sd = [real2angle(sr,angle);sd(k2)];
%     sd = [sr;sd(k2)];
    if sum(isnan(sd))>0
        disp('pause');
    end
    
    sx2 = sx*F + mvnrnd([0 0],R,N);
    beta1_t(:,neuron) = mag*cos(sd);
    beta2_t(:,neuron) = mag*sin(sd);

    k = alpha_t + repmat(sx2(:,1),1,C).*beta1_t + ...
        repmat(sx2(:,2),1,C).*beta2_t;
    lambda3 = repmat(pl,N,1).*exp(k);
    tmp = (lambda3>=1);
    lambda3(tmp) = 1-eps(1);
    
%    sw = PP_likelihood_bychen(Y(t,:),lambda3,'poisson');
     sw = PP_likelihood_bychen(Y(t,:),lambda3,'bino');
    
    %sw = exp(sw-max(sw))./wX;
    sw = sw./wX;
    sw = sw/sum(sw);
    
    PD(t) = sw*sd;
    X(t,:) = sw*sx2;
    
    kk = resample(sw);
    SD(:,t) = sd(kk);
    sx = sx2(kk,:);
    sd = sd(kk,:);
    sw = 1/N*ones(1,N);
    
    sr = angle2real(sd,angle);
%     sr = sd;
    if sum(isinf(sr))>0
        disp('pause');
    end
end
