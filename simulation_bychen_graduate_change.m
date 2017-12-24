%clearvars -except PDs -global;
%clc
%close all;
%% simulation
% %????????????????????????
% load('data\data1\data1_1\B04-PMD-am001-03-binned-10ms','pos_bin');
% X = pos_bin(:,20001:50000)';
% X = InterpPos(X);
% %??????
% fn = 100;
% ap=0.1;
% as=100;
% wp=0.1;
% ws=10; %?????????????????????
% wpp=wp/(fn/2);wss=ws/(fn/2); %?????????;
% [n,wn]=buttord(wpp,wss,ap,as); %????????????????????????
% [b,a]=butter(n,wn); %??????N?????????????????????????????????????????????????????????????????????????????????b,a
% % figure(2)
% % freqz(b,a,512,fn);%??????H(z)??????????????????
% X2=filter(b,a,X);
% % figure(3)
% % plot(X2(:,1),'r','LineWidth',2);
% % hold on
% % plot(X(:,1))
% 
% X = X2(5001:25000,:);
% X = guiyi(X,'2')';
% T = length(X);

% ??????????????????1
% N = 1000;
% T = 20000;
% C = 3;
% delta = 2*pi/N;
% t = 0:delta:2*pi*T/N-delta;
% R = 2.5e-5;
% T = length(t);
% X = 1*[cos(t)' sin(0.3*t)'];
% X = X + sqrt(R)*randn(T,2);

% ??????????????????2
%??????????????????
N = 5000/2;
T = 5000;
delta = 2*pi/N;
t = 0:delta:2*pi*T/N-delta;
R = 1e-3;
T = length(t);
X = repmat((1-2*abs(cos(t).*sin(t)).^2)',1,2).*[cos(t)' sin(t)'];
X = X + sqrt(R)*randn(T,2);

% freq = 1./[500 900 1000 1200];
% X2 = [0:freq(1):1 1-freq(1):-freq(1):-1];
% for i = 1:length(freq)
%    X2 = [X2 -1+freq(i):freq(i):1 1-freq(i):-freq(i):-1]; 
% end
% X2 = X2(1:10000);
% freq = 1./[1500 600 800 1000];
% X3 = [0:freq(1):1 1-freq(1):-freq(1):-1];
% for i = 1:length(freq)
%    X3 = [X3 -1+freq(i):freq(i):1 1-freq(i):-freq(i):-1]; 
% end
% X3 = X3(5001:15000);
% X2 = [X2;X3];
% 
% R = 2.5e-5;
% T = length(X2);
% X2 = X2 + sqrt(R)*randn(2,T);
% X2 = X2(:,1:10000)';
% T = 20000;
% 
% X = [X; X2];
% X(4000:8000) = X(4000:8000)*0.5;
% X(10000:14000) = X(10000:14000)*0.8;

%%%%??????????????????
% C = 8;

% %??????????????????-v1
% % w1 = rand(1,10)*2-1;
% % w1 = repmat(w1,4000,1);
% % k1 = w1(:);
% % k2 = repmat(rand(1,1)*2-1,40000,1);
% % k3 = repmat(rand(1,1)*2-1,40000,1);
% % K1 = [k1,k2,k3]';
% load('E:\algorithms\mcpp3_byxu\bindata\data1\data1_14\1-16_binned_pm001_10ms-tuning1','K2');
% for i = 1:C
%     K2(2,i,:) = 0;
%     K2(3,i,:) = 0;
% end
% % point = [1 1354 5534 8634 11523 16527 18123 23667 26217 29632 31020 34663 37666 40000
% %      1 1286 4355 7636 9739 10945 12965 17532 20052 22713 26495 32767 37515 40000
% %      1 2832 5692 7811 9768 13096 15483 22631 22915 26706 28575 37516 38237 40000
% %      1 1980 4765 8181 11811 15189 17641 20137 24019 27183 32045 36609 38363 40000
% %      1 1013 4864 9321 13027 15467 19211 22859 25791 28406 31662 34207 38788 40000
% %      1 3659 5757 9997 12312 16365 19319 22403 26805 28467 31565 35700 38479 40000
% %      1 2980 5765 8781 12811 16189 19641 21137 23019 26183 30045 33609 37363 40000
% %      1 1486 4355 8636 12739 16945 18965 21532 23052 27713 31495 36767 38515 40000
% %      1 2659 5757 8997 11312 16365 18319 21403 25805 27467 32565 34700 37479 40000
% %      1 957 5113 10154 15605 19307 22053 27692 28430 32381 34862 36951 38502 40000];
% point = [1 4354 7534 10000;
%         1 3286  6636 10000;
%         1 2800 7500 10000;
%         1 3000 8000 10000;
%         1 5000 7900 10000];
% va = [2 3 2.5;
%     2.8 2 2.4;
%     1.7 3 2.4;
%     2.2 3 2;
%     2.4 3 1.8];
% for i = 1:C
%     for j = 1:size(point,2)-1
%         K2(1,i,point(i,j):point(i,j+1)) = va(i,j);
%     end
% end
% K2 = K2(:,:,1:T);
% X1 = [X,ones(length(X),1)];
% lambda = zeros(T,C);
% pl = 0.01;
% for i = 1:C
%     tmpK = squeeze(K2(:,i,:))';
%     lambda(:,i) = pl*exp(sum(X1.*tmpK,2));
% end
% rnd = rand(size(lambda));
% Y = lambda>rnd;

% %??????????????????-V2
% C = 1;
% mag = 6.5;
% pl = 1e-3;
% alpha = 0;
% beta = zeros(2,T);
% lambda = zeros(T,C);
% 
% chpt1 = 4354;
% chpt2 = 7534;
% % plot(chpt1,1:500,'linewidth',2);
% % plot(chpt2,1:500,'linewidth',2);
% 
% w0 = 10/180*pi;
% w1 = 170/180*pi;
% w2 = 60/180*pi;
% ws = [w0*ones(1,chpt1) w1*ones(1,chpt2-chpt1)  w2*ones(1,T-chpt2)];
% K2 = zeros(3,C,T);
% for t = 1:T
%     beta(:,t) = mag*[cos(ws(t)); sin(ws(t))];
%     lambda(t,:) = exp(alpha+X(t,:)*beta(:,t));
%     K2(1,:,t) = repmat(beta(1,t),C,1);
%     K2(2,:,t) = repmat(beta(2,t),C,1);
%     K2(3,:,t) = alpha;
% end
% 
% rndnum = rand(size(lambda));
% Y = lambda*pl>rndnum;

%??????????????????-V3
mag = 1;
pl = repmat([4e-2 4e-2 4.1e-2 3.1e-2 9e-2 9e-2 8.1e-2 6.1e-2 4e-2 4e-2 4.1e-2 3.1e-2 6e-2 6e-2 6.1e-2 6.1e-2],T,1);
pl = repmat(exp(-1.2),T,16);
alpha = 0.2;
% chapt = [1 3000 6000 10000;
%     1 4000 7000 10000;
%     1 2500 6500 10000;
%     1 3500 7500 10000;
%     1 2800 7300 10000];
% w = [50/180*pi 150/180*pi 60/180*pi;
%     150/180*pi 20/180*pi 90/180*pi;
%     10/180*pi 120/180*pi 175/180*pi;
%     45/180*pi 160/180*pi 20/180*pi;
%     90/180*pi 170/180*pi 30/180*pi];
chapt = [1 4000 12000 20000;
    1 8000 14000 20000;
    1 5000 13000 20000;
    1 10000 14000 20000;
    1 6000 12000 20000;
    1 8000 14000 20000;
    1 5000 13000 20000;
    1 10000 14000 20000;
    1 6000 12000 20000;
    1 8000 14000 20000;
    1 5000 13000 20000;
    1 10000 14000 20000;
    1 6000 12000 20000;
    1 8000 14000 20000;
    1 5000 13000 20000;
    1 10000 14000 20000];
chapt = chapt./repmat([1 4 4 4],16,1);
% w = [100/180*pi 40/180*pi 80/180*pi;
%     10/180*pi 80/180*pi 70/180*pi;
%     20/180*pi 80/180*pi 30/180*pi;];
w = [90/180*pi 150/180*pi 30/180*pi;
    135/180*pi 135/180*pi 135/180*pi;
    225/180*pi 225/180*pi 225/180*pi;
    315/180*pi 315/180*pi 315/180*pi;
    0/180*pi 0/180*pi 0/180*pi;
    90/180*pi 90/180*pi 90/180*pi;
    180/180*pi 180/180*pi 180/180*pi;
    270/180*pi 270/180*pi 270/180*pi;
    20/180*pi 20/180*pi 20/180*pi;
    70/180*pi 70/180*pi 70/180*pi;
    110/180*pi 110/180*pi 110/180*pi;
    160/180*pi 160/180*pi 160/180*pi;
    200/180*pi 200/180*pi 200/180*pi;
    250/180*pi 250/180*pi 250/180*pi;
    290/180*pi 290/180*pi 290/180*pi;
    340/180*pi 340/180*pi 340/180*pi;];
C = size(w,1);
% w = [45/180*pi 90/180*pi 135/180*pi;
%     135/180*pi 45/180*pi 90/180*pi;
%     30/180*pi 70/180*pi 45/180*pi;
%     125/180*pi 60/180*pi 90/180*pi;
%     0/180*pi 0/180*pi 0/180*pi;
%     90/180*pi 90/180*pi 90/180*pi;
%     180/180*pi 180/180*pi 180/180*pi;
%     360/180*pi 360/180*pi 360/180*pi;];
%% graduate change 

for i = 1:C
    if i > 1
        for j =1:size(chapt,2)-1
            ws(i,chapt(i,j):chapt(i,j+1)) = w(i,j);
        end
    else %% graduate change for first neuron
        for j = chapt(i,1):chapt(i,2)
            ws(i,j) = w(i,1)-(j-chapt(i,1))*(8/180*pi)/(chapt(i,2)-chapt(i,1));
        end
        for j = chapt(i,2):chapt(i,3)
            ws(i,j) = w(i,2)-(j-chapt(i,2))*(4/180*pi)/(chapt(i,3)-chapt(i,2));
        end
        for j = chapt(i,3):chapt(i,4)
            ws(i,j) = w(i,3)+(j-chapt(i,3))*(16/180*pi)/(chapt(i,4)-chapt(i,3));
        end      
    end
end
%%
% alpha = repmat(alpha,10000,1);
beta1 = mag*cos(ws)';
beta2 = mag*sin(ws)';   
lambda = pl.*exp(alpha+repmat(X(:,1),1,C).*beta1+repmat(X(:,2),1,C).*beta2);
rndnum = rand(size(lambda));
Y = lambda>rndnum;
%% decode
%train
trainT = 1:T;
tmpy = X(trainT(2:end),:);
tmpx = X(trainT(1:(end-1)),:);
F = (tmpx'*tmpx)\(tmpx'*tmpy);%%%% train F
tmpz = tmpx*F; %% prediction
R = var(tmpy-tmpz); %%%% train R

testT = 1:T;
testX = X(testT,:);
testY = Y(testT,:);
% for ii = 12:4:16
%     subNeu = 1:ii;
%     % % static decode
%     % % K = mean(squeeze(K2(:,:,testT)),3);
%     % % K = K2(:,:,testT);
%     % subNeu = 1:1;
%     % % XS = static_ParticleFilter_bychen(testY(:,subNeu),testX(1,:),cov(testX),F,R,K(:,subNeu,:),1000,mag);
%     % 
%     w0 = ws(subNeu,testT(1))';
%     XS = static_ParticleFilter_bychen2(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,w0,10000,mag,pl(1,subNeu));
%     [CCstatic, RMSEstatic] = evaluate(XS,testX);
% 
%     %dynamic decode
%     % testT = 1:20000;
%     % testX = X(testT,:);
%     % testY = Y(testT,:);
%     % subNeu = 2;
%     XD = dynamic_ParticleFilter_bychen2(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,ws(subNeu,:),1000,mag,pl(1,subNeu));
%     [CCdynamic, RMSEdynamic] = evaluate(XD,testX);
%     save(['data\data2\',num2str(ii),'Neu']);
% end

%decode
neuron = 1;

 global T0
 global GLOBAL
if GLOBAL
    p = 0.0001*T0; 
else
%     p = 0.0001*20;
    p=2/T;
%     p=0.003;
end


% K1t = repmat((K2(1,subNeu,testT(1))),length(testT),1);
% K2t = repmat((K2(2,subNeu,testT(1))),length(testT),1);
% K3t = repmat((K2(3,subNeu,testT(1))),length(testT),1);
% top = max(squeeze(K2(1,neuron,testT)));
% bot = min(squeeze(K2(1,neuron,testT)));
% [PX,PD] = PP_ParticleFilter_bychen(testY(:,subNeu),testX(1,:),cov(testX),F,R,K1t,K2t,K3t,1000,p,neuron,top,bot,pl);
% realK1 = squeeze(K2(1,neuron,testT));
% Ktmp = K2;
% Ktmp(1,neuron,:) = PD;
% for i = 1:C
%     tmpK = squeeze(Ktmp(:,i,:))';
%     lambda2(:,i) = pl*exp(sum(X1.*tmpK,2));
% end
top = max(ws(neuron,:));
bot = min(ws(neuron,:));
%for ii = 4:4:16
    ii = 16;
    subNeu = 1:ii;
    subNeu = 1:4;
    
    
    % frequency of Y(:,1)
    Y_bar = zeros(size(Y,1),C);
    global period
    
    for i = 1:T
        if i > T-period
            Y_bar(i,:) = sum(Y(T-2*(T-i+1)+1:T,:),1)/(T-i+1)*period/2;
        elseif i < period
            Y_bar(i,:) = sum(Y(1:2*i,:),1)/i*period/2;
        else 
            Y_bar(i,:) = sum(Y(i-period+1:i+period,:),1)/2;
        end
    end
    
 
% for i = 1:T
%     if i > T-2*period
%         Y_bar(i,:) = sum(Y(T-2*period+1:T,:),1)/2;
% 
%     else 
%         Y_bar(i,:) = sum(Y(i:i+2*period-1,:),1)/2;
%     end
% end
      

    [PX,PD] = PP_ParticleFilter_bychen2_modified(Y_bar(:,subNeu),testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,ws(subNeu,testT(1))',1000,p,neuron,top,bot,pl(1,subNeu),mag);
    %[PX,PD] = PP_ParticleFilter_bychen2_modified_global2(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,ws(subNeu,testT(1))',1000,p,neuron,top,bot,pl(1,subNeu),mag);
%     [PX,PD] = PP_ParticleFilter_bychen2_modified_global(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,ws(subNeu,testT(1))',N,p,neuron,top,bot,pl(1,subNeu),mag);
    %[PX,PD] = PP_ParticleFilter_bychen2_modified(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,ws(subNeu,testT(1))',N,p,neuron,top,bot,pl(1,subNeu),mag);
%     [PX,PD] = PP_ParticleFilter_bychen2(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,ws(subNeu,testT(1))',1000,p,neuron,top,bot,pl(1,subNeu),mag);
    %[PX,PD] = PP_ParticleFilter_bychen2_only_LW(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,ws(subNeu,testT(1))',10000,p,neuron,top,bot,pl(1,subNeu),mag);
    %[CCadaptive, RMSEadaptive] = evaluate(PX,testX);
    %save(['data\data3\',num2str(ii),'Neu']);
%end


%plot
%{
plot(PD);
hold on
plot(ws(1,:),'r'); %% ws(c,t):= true theta of neuron c at time t.
legend('adaptive','grandtruth');
%}


%{
global TRIAL
global T0
global PHI1
global PHI2
global ETA
global GLOBAL
if GLOBAL
    global_sampling = 'TRUE';
else
    global_sampling = 'FALSE';
end
%}
%title(['phi1:',string(PHI1),';phi2:',string(PHI2),';eta:',string(ETA),';global:',string(GLOBAL),';t0(if global):',string(T0)])

%{
figure(2)
plot(testX(:,1),'r');hold on;plot(PX(:,1));
legend('grandtruth', 'adaptive');
figure(3);
plot(testX(:,2),'r');hold on; plot(PX(:,2));
legend('grandtruth', 'adaptive');
%}
%{
figure(2) %%problem with XS: cleared?
plot(testX(:,1),'r');hold on;plot(XS(:,1),'k');hold on;plot(XD(:,1),'m');hold on;plot(PX(:,1));
legend('grandtruth','static','dynamic','adaptive');
figure(3)
plot(testX(:,2),'r');hold on;plot(XS(:,2),'k');hold on;plot(XD(:,2),'m');hold on;plot(PX(:,2));
legend('grandtruth','static','dynamic','adaptive');
%}

% for i = 1:length(subNeu)
%     figure(i)
%     plot(PDsing(:,i));
%     hold on
%     plot(ws(i,:));
% end
% PX2 = dynamic_ParticleFilter_bychen2(testY(:,subNeu),testX(1,:),cov(testX),F,R,alpha,PDsing',1000,mag,pl);
% figure(i+1)
% subplot(2,1,1);plot(PX2(:,1));hold on;plot(testX(:,1),'r');hold on;plot(PXsing(:,1,1),'k');hold on;plot(PXsing(:,1,2),'g');
% subplot(2,1,2);plot(PX2(:,2));hold on;plot(testX(:,2),'r');hold on;plot(PXsing(:,2,1),'k');hold on;plot(PXsing(:,2,2),'g');
% 
% beta1 = mag*cos(PD2);
% beta2 = mag*sin(PD2);
% lambda2 = pl*exp(alpha+repmat(PX2(:,1),1,1).*beta1+repmat(PX2(:,2),1,1).*beta2);
% %plot
% plot(PD2);
% hold on
% plot(ws(1,:),'r');
% saveas(gcf,['data\Neu',num2str(length(subNeu)),'_PD']);
% figure(2)
% plot(lambda2);
% hold on
% plot(lambda(:,1),'r');
% saveas(gcf,['data\Neu',num2str(length(subNeu)),'_Lambda']);
% save(['data\Neu',num2str(length(subNeu)),'_data']);





