%% Numerical solution of GPE using TSSP
% Ref: Bao et al. Jour of Comp Phys 187 (2003) 318-342

clear all; close all; clc;

vareps = 0.1;
kappa1 = 1.2649;
a = -16; b =16;
M = 16*32;
h = 1/16;          % h is the mesh size: h = (b-a)/M
k = 0.02;

x = a:h:b;
% x = (a+h/2):h:(b-h/2);
l = -M/2:(M/2-1);
j_id = 0:(M-1);
mu_l = 2*pi*l/(b-a);

psi_x = 1/((pi*vareps)^(1/4))*exp(-x.^2/2/vareps);
Nk = 2000;
idd = intersect(find(x<8), find(x>-8));

abc = zeros(length(idd),Nk);
sigm_w = zeros(size(1:Nk));

for kk = 1:Nk
        psi_star = exp(-1i* (x.^2/2+kappa1*(psi_x).*(conj(psi_x)))...
                            *k/2/vareps)                            .*psi_x;

        psi_starHat = zeros(size(l));
        psi_doubleStar = zeros(size(j_id));

        for ll = 1:length(l)
            psi_starHat(ll) = sum(psi_star(1:M).*exp(-1i*mu_l(ll)*(x(1:M)-a)));
        end

        for jj = 1:length(j_id)
            psi_doubleStar(jj) = 1/M*sum(exp(-1i*vareps*k*mu_l.^2/2).*psi_starHat.*exp(1i*mu_l*(x(jj)-a)));
            psi_x(jj) = exp(-1i*(x(jj)^2/2+kappa1*abs(psi_doubleStar(jj))^2)*k/2/vareps)*psi_doubleStar(jj);
        end
        psi_x(M+1)=psi_x(1);
%         psi_x = exp(-1i*(x.^2/2+kappa1*abs(psi_doubleStar).^2)*k/2/vareps).*psi_doubleStar;
        abc(:,kk) = abs(psi_x(idd)).^2;
        aa = abs(psi_x).^2;
        
        x_av = sum(x.*aa*h);
        sigma_sq = sum(((x-x_av).^2).*aa*h);
% %         
% %         ia = max(aa);
% %         idx = intersect(find(aa<(ia/exp(1)+0.005)),find(aa>(ia/exp(1)-0.005)));
        sigm_w(kk) = sqrt(sigma_sq);
end

figure; plot((1:Nk)*k,sigm_w,'b-');
% figure;plot(x(201:313),abs(psi_x(201:313)).^2)