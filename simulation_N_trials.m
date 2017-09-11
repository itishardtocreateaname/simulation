clear all;

%%
global TRIAL
global PHI1
global PHI2
global ETA
global T0
global GLOBAL
global MOD
global period

TRIAL = 9;
PHI1 = 0.90;
PHI2 = 0.9975;
ETA = 0.0025;
T0 = 200;
GLOBAL = false;
MOD = 0;
period = 70;

%% calculate PDs && PXs
tic;
% simulation_bychen;
simulation_bychen_graduate_change;

N = size(PD,1);
PDs = zeros(N, TRIAL);

[xsize1, xsize2] = size(PX);
PXs = zeros(xsize1, xsize2, TRIAL);

PDs(:, 1) = PD;
PXs(:,:,1) = PX;

if TRIAL > 1
    for tr = 2:TRIAL
%         clearvars -except TRIAL PHI1 PHI2 ETA T0 GLOBAL PDs tr MOD period PXs
        MOD = (tr-1)*round(T0/TRIAL);
%         simulation_bychen;
        simulation_bychen_graduate_change;
        PDs(:, tr) = PD;
        PXs(:,:,tr) = PX;
    end
end
%% calculate PD_mean && PX_mean && CC && MSE
PD_mean = mean(PDs, 2);
PX_mean = mean(PXs, 3);
ground_truth = ws(1,:)';

MSE_theta = immse(PD_mean, ground_truth);
MSE_x = immse(PX_mean, testX);

CC_theta = corrcoef(PD_mean, ground_truth);
CC_x = corrcoef(PX_mean, testX);
%% MSE vs TRIAL_num %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSE = zeros(1, TRIAL); %%MSE(tr): mse over whole time period during trial 1~tr.
% for tr = 1:TRIAL
%     PD_mean_tr = mean(PDs(:, 1:tr), 2); %% average value of PD from trial 1 to tr.
%     MSE(tr) = (PD_mean_tr - ground_truth)' * (PD_mean_tr - ground_truth) / size(PD_mean, 1);
% end
% plot(MSE, '.-');
% if GLOBAL
%     title(sprintf('T0=%d;phi1=%.4f;phi2=%.4f;eta=%.4f',T0,PHI1,PHI2,ETA));
% else
%     title(sprintf('phi1=%.4f;phi2=%.4f;eta=%.4f',PHI1,PHI2,ETA));
% end
% 
% xlabel('Trial Num');
% ylabel('MSE');
%%
% if GLOBAL
%     savefig(sprintf('MSE_Global_T0=%d_Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.fig',T0, TRIAL,PHI1,PHI2,ETA));
% else
%     savefig(sprintf('MSE_Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.fig', TRIAL,PHI1,PHI2,ETA));
% end

%% AVERAGE-PD with error bar over trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ground_truth);
hold on;
plot(PD_mean);

% S = std(PDs,0,2);
% kk = [1:500:6000,6000:250:8000,8000:500:T];
% errorbar(kk,PD_mean(kk), S(kk));

% PD_beta1=mag*cos(PDs(:,end));
% PD_beta2=mag*sin(PDs(:,end));
% PD_lambda=pl(:,1).*exp(alpha+PX(:,1).*PD_beta1+PX(:,2).*PD_beta2);
% plot(PD_lambda);
% plot(lambda(:,1));
% plot(Y(:,1));

if GLOBAL
    title(sprintf('T0=%d;Trial=%d;phi1=%.4f;phi2=%.4f;eta=%.4f',T0, TRIAL, PHI1,PHI2, ETA));
else
    title(sprintf('Trial=%d;phi1=%.4f;phi2=%.4f;eta=%.4f', TRIAL, PHI1,PHI2, ETA));
end
xlabel('Time');
ylabel('Theta');

hold off;
%%
if GLOBAL
    savefig(sprintf('ADAM_theta_Global_T0=%d_Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.fig', T0, TRIAL, PHI1,PHI2, ETA));
else
    savefig(sprintf('ADAM_theta_Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.fig', TRIAL, PHI1,PHI2, ETA));
end

%% plot PX_mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(testX);
hold on;
plot(PX_mean);
if GLOBAL
    title(sprintf('T0=%d;Trial=%d;phi1=%.4f;phi2=%.4f;eta=%.4f',T0, TRIAL, PHI1,PHI2, ETA));
else
    title(sprintf('Trial=%d;phi1=%.4f;phi2=%.4f;eta=%.4f', TRIAL, PHI1,PHI2, ETA));
end
xlabel('Time');
ylabel('Theta');

hold off;
%%
if GLOBAL
    savefig(sprintf('ADAM_X_Global_T0=%d_Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.fig', T0, TRIAL, PHI1,PHI2, ETA));
else
    savefig(sprintf('ADAM_X_Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.fig', TRIAL, PHI1,PHI2, ETA));
end
%%
time_period = toc;
%% save workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GLOBAL
    save(sprintf('T0=%d_Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.mat',T0,TRIAL,PHI1,PHI2,ETA));
else
    save(sprintf('Trial=%d_phi1=%.4f_phi2=%.4f_eta=%.4f.mat',TRIAL,PHI1,PHI2,ETA));
end
%% find min mse
% mse = zeros(1,TRIAL);
% for i = 1:TRIAL
%     mse(i) = immse(PDs(:,i), ground_truth);
% end
% [m,i] = min(mse);
% 
% plot(ground_truth);
% hold on;
% plot(PDs(:,i));
% 
