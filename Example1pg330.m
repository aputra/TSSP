%% Numerical solution of GPE using TSSP -- one dimensional
% REF: Bao et al. Jour of Comp Phys 187 (2003) 318-342
% Benchmark to Example 1
% -- Fig 1a: I believe there is some error in the paper, where the
%            computation is actually done up to t = 5 (not t = 4).

clear all; close all; clc;

vareps = 0.1;
kappa1 = 1.2649;
k = 0.01;
a = -16; b = 16;

% -- meshes
h = 1/16;               % h is the mesh size: h = (b-a)/M
M = round(1/h*(b-a));

x = a:h:b;
l = -M/2:(M/2-1);
j_id = 0:(M-1);
mu_l = 2*pi*l/(b-a);


% -- initial condition
% S0 = -log(exp(x)+exp(-x));
% psi_x = exp(-x.^2).*(1+0.2*cos(x/sqrt(vareps)).*cos(x/sqrt(vareps))).*exp(1i/vareps.*S0);
psi_x = 1/(pi)^(1/4)*exp(-x.^2/2);

Nk = 40/k;
idd = intersect(find(x<=3), find(x>=0));

psi2_cut = zeros(length(idd),Nk);
sigm_w = zeros(size(1:Nk));
x_com = zeros(size(1:Nk));

for kk = 1:Nk
    psi_star = exp(-1i* (x.^2/2 + kappa1*abs(psi_x).^2)...
                        *k/2/vareps)                            .*psi_x;

    psi_starHat = psi_star(1:M)*exp(-1i*(x(1:M)-a)'*mu_l);
    psi_doubleStar = 1/M * (exp(-1i*vareps*k*mu_l.^2/2).*psi_starHat) * exp(1i*mu_l'*(x(1:M)-a));
    
    psi_x(1:length(j_id)) = exp(-1i*(x(1:M).^2/2+kappa1*abs(psi_doubleStar).^2)*k/2/vareps).*psi_doubleStar;
    psi_x(M+1)=psi_x(1);

    psi2_cut(:,kk) = abs(psi_x(idd)).^2;
    psi2 = abs(psi_x).^2;
    
    % -- condensate width (see eq 4.2)
    x_av = sum(x.*psi2*h);
    x_com(kk) = x_av;           % snapshot of x_av during the TSSP
    sigma_sq = sum(((x-x_av).^2).*psi2*h);

    sigm_w(kk) = sqrt(sigma_sq);
    
%     figure(1);
%         plot(x(idd),psi2_cut(:,kk));
    
end

figure(2); plot((1:Nk)*k,sigm_w,'b-');
figure(3); plot(x(idd),psi2_cut(:,kk),'b+-');

% -- plot result as in fig 2
time_idx = round(20/k):round(30/k);
figure(4);surf(x(idd),time_idx*k,psi2_cut(:,time_idx).');view(0,90);shading interp