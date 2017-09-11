function [sr_t, m_t, v_t, g_t] = Adam(sr_pt, angle, m_pt, v_pt,tau, y_t, x_t, phi1, phi2, eta, pl, alpha_t, mag, beta1_t, beta2_t, neuron)
%% *_t: updated value; *_pt: previous timepoint value
%% m: estimation of first-order moment; v: estimation of second-order moment
%% y: the observation
%% x: estimated value according to x_(t-1)^i (i=1,2,...,N) and \sr_(t-1)
%% phi1: forget rate of gradient
%% phi2: forget rate of square gradient
%% eta: learning rate
%% sr,m,v,g: 1*1, y: 1*C, x: 1*2, phi,eta: 1*1, pl: 1*C, mag: 1*1, alpha_t,beta_t:1*C, neuron:1*1


g_t = gradient(sr_pt, angle, y_t, x_t, pl, alpha_t, mag, beta1_t, beta2_t, neuron);
m_t = phi1*m_pt + (1-phi1)*g_t;
v_t = phi2*v_pt + (1-phi2)*g_t.*g_t;
global T0
global GLOBAL
if tau <= T0 && GLOBAL
    m_t_hat = m_t/(1-phi1^tau);
    v_t_hat = v_t/(1-phi2^tau);
else
    m_t_hat = m_t;
    v_t_hat = v_t;
end
sr_t = sr_pt - eta * m_t_hat./(sqrt(v_t_hat)+eps(1));

%sr_t = sr_pt - eta * m_t;%/(sqrt(v_t)+eps(1));
