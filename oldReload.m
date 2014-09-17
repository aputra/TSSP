%% Numerical solution of GPE using TSSP - Two dimensional
% Ref: Bao et al. Jour of Comp Phys 187 (2003) 318-342

clear all; close all; clc;

vareps = 1;
gamma_y = 2;
kappa_2 = 0.1;
ax = -4; bx = 4;
ay = -4; by = 4;

h = 1/8;           % h is the mesh size: h = (b-a)/M
Mx = round(1/h*(bx-ax));
My = round(1/h*(by-ay));
k = 0.02;          % k is the time step interval

x = ax:h:bx;
y = ay:h:by;

lx = -Mx/2:(Mx/2-1);
ly = -My/2:(My/2-1);

j_idx = 0:(Mx-1);
j_idy = 0:(My-1);

mu_lx = 2*pi*lx/(bx-ax);
mu_ly = 2*pi*ly/(by-ay);

[yy,xx] = meshgrid(y,x);
% psi_xy = 1/sqrt(pi*vareps)*exp(-(xx.^2+yy.^2)/2/vareps);
S0 = cosh(sqrt(xx.^2 + 2*yy.^2));
psi_xy = gamma_y^(1/4)/sqrt(pi)*exp(-(xx.^2+gamma_y*yy.^2)/2/vareps).*exp(-1i/vareps*S0);

Nk = 4/k;

% % idd = intersect(find(x<=3), find(x>=0));
% % 
% % abc = zeros(length(idd),Nk);
sigmX = zeros(size(1:Nk));
sigmY = zeros(size(1:Nk));

x_com = zeros(size(1:Nk));
y_com = zeros(size(1:Nk));

psi_starHat = zeros(length(lx),length(ly));
psi_doubleStar = zeros(length(j_idx),length(j_idy));


for kk = 1:Nk
    psi_star = exp(-1i* (xx.^2/2 + gamma_y^2*yy.^2/2 + kappa_2*abs(psi_xy).^2)...
                                                       *k/2/vareps)             .*psi_xy;

    for ii = 1:length(lx)
        for jj = 1:length(ly)
            
            psi_starHat(ii,jj) = sum(sum(psi_star(1:Mx,1:My) .* ((exp(-1i*(x(1:Mx)-ax)*mu_lx(ii))).' * (exp(-1i*(y(1:My)-ay)*mu_ly(jj))))));
        end
    end
    
    
    for ii = 1:length(j_idx)
        for jj = 1:length(j_idy)
            psi_doubleStar(ii,jj) = 1/Mx * 1/My * sum(sum(((exp(-1i*vareps*k*mu_lx.^2/2+1i*mu_lx*(xx(ii,jj)-ax))).' * ...
                                                  (exp(-1i*vareps*k*mu_ly.^2/2+1i*mu_ly*(yy(ii,jj)-ay)))) .* psi_starHat));
                                              
            psi_xy(ii,jj) = exp(-1i*(x(ii)^2/2 + gamma_y^2*y(jj)^2/2+kappa_2*abs(psi_doubleStar(ii,jj)^2))*k/2/vareps)*psi_doubleStar(ii,jj);
        end
    end
    
    
    
    psi_xy(Mx+1,:)=psi_xy(1,:);
    psi_xy(:,My+1)=psi_xy(:,1);
    
    aa = abs(psi_xy).^2;

    x_av = sum(sum(xx.*aa*h*h));
    y_av = sum(sum(yy.*aa*h*h));
    x_com(kk) = x_av;
    y_com(kk) = y_av;
    
    sigmax_sq = sum(sum(((xx-x_av).^2).*aa*h*h));
    sigmay_sq = sum(sum(((yy-y_av).^2).*aa*h*h));
    
    sigmX(kk) = sqrt(sigmax_sq);
    sigmY(kk) = sqrt(sigmay_sq);
%     figure(1);
%         plot(x(idd),abc(:,kk));
    
end


figure(2); plot((1:Nk)*k,sigmX,'b-');
figure(2); hold on; plot((1:Nk)*k,sigmY,'b--');
% figure(3); plot(x(idd),abc(:,kk),'b+-'); 