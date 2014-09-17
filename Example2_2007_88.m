%% TSSP for coupled GPE with applications to rotating BECs
% Ref: Wang, J Comp and App Math 205 (2007) 88

clear all; close all; clc;

Omega = 0;

% % % gamma_x = 1;
% % % gamma_y = 1.5;

delta_t = 0.01;          % time step interval

ax = -10; bx = 10;    % xmin and xmax
cy = -10; dy = 10;    % ymin and ymax

h = 1/10;           % h is the mesh size: h = (bx-ax)/Mx = (dy-cy)/Ny
Mx = round((bx-ax)/h);
Ny = round((dy-cy)/h);

x_j = (ax+h):h:(bx-h);
y_k = (cy+h):h:(dy-h);

% FFT summation indices for x and y spaces:
px = (-Mx/2+1):(Mx/2-1);
qy = (-Ny/2+1):(Ny/2-1);

j_idx = 0:(Mx-1);
m_idy = 0:(Ny-1);

% FFT reciprocal space
mu_p = 2*pi*px/(bx-ax);
la_q = 2*pi*qy/(dy-cy);

[yy,xx] = meshgrid(y_k,x_j);

% Initial wavefunction
psi_ho = 1/sqrt(pi)*exp(-(xx.^2+yy.^2)/2);                     % xj * yk
psi_v = (xx+1i*yy)/sqrt(pi).*exp(-(xx.^2+yy.^2)/2);                     % xj * yk
psi_1 = ((1-Omega)*psi_ho + Omega*psi_v);
psi_1 = psi_1/sqrt(sum(sum(abs(psi_1).^2))*h*h);
% psi_1 = psi_v;

tTotal = round(3.5/delta_t);
N1 = zeros(1,tTotal);

delta1 = zeros(1,tTotal);

beta2D = 1622.9;
V2D_1 = 1/2*(xx.^2 + yy.^2); % xj * yk
idV0 = find(V2D_1==0);

for tt = 1:tTotal

    % First step
    psi_hat1 = FFT_x_sine(psi_1,Mx,Ny,x_j,mu_p,ax); % Ny * px
    psi_1_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    % Second step
    psi_hat1 = FFT_y_sine(psi_1_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy    
    psi_1_stepNew = rotatingIFFT_y_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk

    % Third step
    psi_1_step = sqrt(V2D_1.*exp(-2*delta_t*V2D_1)./ (V2D_1 + beta2D*abs(psi_1_stepNew).^2 .* (1-exp(-2*delta_t*V2D_1)) )  ) .*psi_1_stepNew;
    psi_1_step(idV0) = 1/sqrt(1+beta2D*2*delta_t*abs(psi_1_stepNew(idV0)).^2).*psi_1_stepNew(idV0);
    
    % Fourth step
    psi_hat1 = FFT_y_sine(psi_1_step,Mx,Ny,y_k,la_q,cy); % Mx * qy
    psi_1_step = rotatingIFFT_y_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk

    % Fifth step --> advancing from t_n to t_n+1 --> end of one loop
    psi_hat1 = FFT_x_sine(psi_1_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
    psi_1_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    psi_1 = psi_1_step.';
    N1(tt) = sum(sum(abs(psi_1).^2))*h*h;
    psi_1 = psi_1/sqrt(N1(tt));
    N1(tt) = sum(sum(abs(psi_1).^2))*h*h;
    delta1(tt) = sum(sum(abs(psi_1).^2 .* (xx.^2 + yy.^2)))*h*h;

%     figure(1); hold on; plot(tt*delta_t, step_diff, 'bx');
    figure(1); surf(y_k,x_j,abs(psi_1).^2); shading interp; view(0,90);
end

figure(3); plot((1:tTotal)*delta_t,N1);
figure(4); plot((1:tTotal)*delta_t,delta1);

save Ex2groundstate.mat psi_1;