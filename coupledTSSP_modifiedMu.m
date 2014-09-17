%% TSSP for coupled GPE with applications to rotating BECs
% Ref: Wang, J Comp and App Math 205 (2007) 88

clear all; close all; clc;

Omega = 0.7;

gamma_x = 1;
gamma_y = 1;

delta_t = 0.05;          % time step interval

ax = -12; bx = 12;    % xmin and xmax
cy = -12; dy = 12;    % ymin and ymax

h = 1/12;           % h is the mesh size: h = (bx-ax)/Mx = (dy-cy)/Ny
Mx = round((bx-ax)/h);
Ny = round((dy-cy)/h);

x_j = (ax+h):h:(bx-h);
y_k = (cy+h):h:(dy-h);

% % x_j = x_j + h/100; y_k = y_k + h/100;
% FFT summation indices for x and y spaces:
px = (-Mx/2+1):(Mx/2-1);
qy = (-Ny/2+1):(Ny/2-1);

j_idx = 0:(Mx-1);
m_idy = 0:(Ny-1);

% FFT reciprocal space
mu_p = 2*pi*px/(bx-ax);
la_q = 2*pi*qy/(dy-cy);

[yy,xx] = meshgrid(y_k,x_j);
% mv = 6;
% Initial wavefunction
psi_ho = (gamma_x*gamma_y)^(1/4)/sqrt(pi).*exp(-(gamma_x*xx.^2+gamma_y*yy.^2)/2);                     % xj * yk
psi_v = (xx+1i*yy)/sqrt(pi).*exp(-(xx.^2+yy.^2)/2);                     % xj * yk
psi_1 = ((1-Omega)*psi_ho + Omega*psi_v);

% psi_1 = psi_ho;


psi_1 = psi_1/sqrt(sum(sum(abs(psi_1).^2))*h*h);

tTotal = round(200/delta_t);
N1 = zeros(1,tTotal);

delta1 = zeros(1,tTotal);
Energy = zeros(1,tTotal);
beta2D = 100;


% V2D_1 = 1/2*(gamma_x^2*xx.^2 + gamma_y^2*yy.^2) + 4*exp(-((xx-1).^2+yy.^2)); % xj * yk
V2D = 1/2*(xx.^2+yy.^2); % xj * yk
idV0 = find(V2D==0);

for tt = 1:tTotal
    N1(tt) = sum(sum(abs(psi_1).^2))*h*h;
    mu_chem = -log(N1(tt))/2/delta_t;
    V2D_1 = V2D - mu_chem;
    % First step
    psi_1_step = sqrt(V2D_1.*exp(-delta_t*V2D_1)./ (V2D_1 + beta2D*abs(psi_1).^2 .* (1-exp(-delta_t*V2D_1)) )  ) .*psi_1;
    if mu_chem == 0 
        psi_1_step(idV0) = 1/sqrt(1+beta2D*delta_t*abs(psi_1(idV0)).^2).*psi_1(idV0);
    end
    
    % Second step
    psi_hat1 = FFT_x_sine(psi_1_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
    psi_1_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    % Third step
    psi_hat1 = FFT_y_sine(psi_1_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy    
    psi_1_stepNew = rotatingIFFT_y_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk

    % Fourth step
    
    V2D_1 = V2D - mu_chem;
    
    psi_1_step = sqrt(V2D_1.*exp(-delta_t*V2D_1)./ (V2D_1 + beta2D*abs(psi_1_stepNew).^2 .* (1-exp(-delta_t*V2D_1)) )  ) .*psi_1_stepNew;
    if mu_chem == 0
        psi_1_step(idV0) = 1/sqrt(1+beta2D*delta_t*abs(psi_1_stepNew(idV0)).^2).*psi_1_stepNew(idV0);
    end
    
    psi_1 = psi_1_step;
    N1(tt) = sum(sum(abs(psi_1).^2))*h*h;
    psi_1 = psi_1/sqrt(N1(tt));
    
    
%     figure(1); hold on; plot(tt*delta_t, step_diff, 'bx');
    figure(1); subplot(1,2,1); surf(y_k,x_j,abs(psi_1).^2); shading interp; view(0,90);
    figure(1); subplot(1,2,2); hold on; plot(tt*delta_t,mu_chem,'bx');
end

figure(3); plot((1:tTotal)*delta_t,N1);
figure(4); plot((1:tTotal)*delta_t,delta1);

% xrms = sqrt(sum(sum(abs(psi_1).^2 .* (xx.^2)))*h*h)
% yrms = sqrt(sum(sum(abs(psi_1).^2 .* (yy.^2)))*h*h)
% sqrt(delta1(tt))
[gx,gy] = gradient(psi_1,h);

% Energy = sum(sum(-conj(psi_1).*del2(psi_1,h) +V2D_1.*abs(psi_1).^2 + beta2D/2.*abs(psi_1).^4) ) *h*h
Energy = sum(sum(1/2*(abs(gx).^2+abs(gy).^2) +V2D_1.*abs(psi_1).^2 + beta2D/2.*abs(psi_1).^4) ) *h*h
mu_chem = Energy + sum(sum(beta2D/2.*abs(psi_1).^4) ) *h*h