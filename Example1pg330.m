%% Numerical solution of GPE using TSSP
% Ref: Bao et al. Jour of Comp Phys 187 (2003) 318-342

clear all; close all; clc;

vareps = 0.025;
kappa1 = 1.2649;
a = -16; b = 16;

h = 1/64;          % h is the mesh size: h = (b-a)/M
M = round(1/h*(b-a));
k = 0.0025;

x = a:h:b;
l = -M/2:(M/2-1);
j_id = 0:(M-1);
mu_l = 2*pi*l/(b-a);

% S0 = -log(exp(x)+exp(-x));
% psi_x = exp(-x.^2).*(1+0.2*cos(x/sqrt(vareps)).*cos(x/sqrt(vareps))).*exp(1i/vareps.*S0);
psi_x = 1/(pi)^(1/4)*exp(-x.^2/2);

Nk = 5/k;
idd = intersect(find(x<=3), find(x>=0));

abc = zeros(length(idd),Nk);
sigm_w = zeros(size(1:Nk));
x_com = zeros(size(1:Nk));

% F(Nk) = struct('cdata',[],'colormap',[]);

for kk = 1:Nk
    psi_star = exp(-1i* (x.^2/2 + kappa1*abs(psi_x).^2)...
                        *k/2/vareps)                            .*psi_x;

    
    psi_starHat = psi_star(1:M)*exp(-1i*(x(1:M)-a)'*mu_l);
    psi_doubleStar = 1/M * (exp(-1i*vareps*k*mu_l.^2/2).*psi_starHat) * exp(1i*mu_l'*(x(1:M)-a));
    
    psi_x(1:length(j_id)) = exp(-1i*(x(1:M).^2/2+kappa1*abs(psi_doubleStar).^2)*k/2/vareps).*psi_doubleStar;
    psi_x(M+1)=psi_x(1);

    abc(:,kk) = abs(psi_x(idd)).^2;
    aa = abs(psi_x).^2;

    x_av = sum(x.*aa*h);
    x_com(kk) = x_av;
    sigma_sq = sum(((x).^2).*aa*h);

    sigm_w(kk) = sqrt(sigma_sq);
%     figure(1);
%         plot(x(idd),abc(:,kk));
    
end


figure(2); plot((1:Nk)*k,sigm_w,'b-');
figure(3); plot(x(idd),abc(:,kk),'b+-');