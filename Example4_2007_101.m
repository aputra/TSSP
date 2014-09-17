%% TSSP for coupled GPE with applications to rotating BECs
% Numerical Simulations on Stationary States for 2 component BECs
% Ref: Wang, J Sci Comp (2009) 38: 149-163 -- 
% This program is tested to some extent, not overall runs match the
% reference.

clear all; close all; clc;

Omega = 0.6;

gamma_x = 1;
gamma_y = 1;

delta_t = 0.002;          % time step interval

ax = -12; bx = 12;    % xmin and xmax
cy = -12; dy = 12;    % ymin and ymax

h = 1/8;           % h is the mesh size: h = (bx-ax)/Mx = (dy-cy)/Ny
Mx = round((bx-ax)/h);
Ny = round((dy-cy)/h);

x_j = (ax+h):h:(bx-h) + h/400;
y_k = (cy+h):h:(dy-h) + h/400;

% FFT summation indices for x and y spaces:
px = (-Mx/2+1):(Mx/2-1);
qy = (-Ny/2+1):(Ny/2-1);

% j_idx = 0:(Mx-1);
% m_idy = 0:(Ny-1);

% FFT reciprocal space
mu_p = 2*pi*px/(bx-ax);
la_q = 2*pi*qy/(dy-cy);

[yy,xx] = meshgrid(y_k,x_j);

% Initial wavefunction
psi_ho = (gamma_x*gamma_y)^(1/4)/sqrt(pi).*exp(-(gamma_x*xx.^2+gamma_y*yy.^2)/2);                     % xj * yk
psi_v = (gamma_x*gamma_y)^(1/4)/sqrt(pi)*(xx+1i*yy).*exp(-(xx.^2+yy.^2)/2);                     % xj * yk
psi_cc = ((1-Omega)*psi_ho + Omega*psi_v);
% psi_1 = psi_v;
norm_fac = sum(sum(abs(psi_cc).^2)*h)*h;

psi_1 = psi_cc/sqrt(norm_fac*2);%*sqrt(10/11);
psi_2 = psi_cc/sqrt(norm_fac*2);%*sqrt(1/11);

tTotal = round(40/delta_t);
N1 = zeros(1,tTotal);
N2 = zeros(1,tTotal);

delta1 = zeros(1,tTotal);
delta2 = zeros(1,tTotal);

beta2D = 2000 * [1.0     0.7;
                 0.7   	 1.0];
            
% % V2D_1 = 1/2*(gamma_x^2*xx(1:Mx,1:Ny).^2 + gamma_y^2*yy(1:Mx,1:Ny).^2 ); % xj * yk
% % V2D_2 = 1/2*(gamma_x^2*xx(1:Mx,1:Ny).^2 + gamma_y^2*yy(1:Mx,1:Ny).^2 ); % xj * yk
            
V2D_1 = 1/2*(gamma_x^2*xx.^2 + gamma_y^2*yy.^2 ); % xj * yk
V2D_2 = 1/2*(gamma_x^2*xx.^2 + gamma_y^2*yy.^2 ); % xj * yk

for tt = 1:tTotal
    
    % First step
    psi_hat1 = FFT_x_sine(psi_1,Mx,Ny,x_j,mu_p,ax); % Ny * px
    psi_hat2 = FFT_x_sine(psi_2,Mx,Ny,x_j,mu_p,ax); % Ny * px
    
    psi_1_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    psi_2_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    % Second step
    psi_hat1 = FFT_y_sine(psi_1_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy
    psi_hat2 = FFT_y_sine(psi_2_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy
    
    psi_1_stepNew = rotatingIFFT_y_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    psi_2_stepNew = rotatingIFFT_y_sine(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    

    % Fourth step
    psi_1_step = exp(-delta_t*(V2D_1 + beta2D(1,1)*abs(psi_1_stepNew).^2 + beta2D(1,2)*abs(psi_2_stepNew).^2)).*psi_1_stepNew;   % xj * yk
    psi_2_step = exp(-delta_t*(V2D_2 + beta2D(2,1)*abs(psi_1_stepNew).^2 + beta2D(2,2)*abs(psi_2_stepNew).^2)).*psi_2_stepNew;   % xj * yk
    
    % Fifth step
    
    % Sixth step
    psi_hat1 = FFT_y_sine(psi_1_step,Mx,Ny,y_k,la_q,cy); % Mx * qy
    psi_hat2 = FFT_y_sine(psi_2_step,Mx,Ny,y_k,la_q,cy); % Mx * qy
    
    psi_1_step = rotatingIFFT_y_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    psi_2_step = rotatingIFFT_y_sine(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    
    % Seventh step --> advancing from t_n to t_n+1 --> end of one loop
    psi_hat1 = FFT_x_sine(psi_1_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
    psi_hat2 = FFT_x_sine(psi_2_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
    
    psi_1_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    psi_2_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    psi_1 = psi_1_step.';
    psi_2 = psi_2_step.';
    
    
    N1(tt) = sum(sum(abs(psi_1).^2)*h)*h;
    N2(tt) = sum(sum(abs(psi_2).^2)*h)*h;
    psi_1 = psi_1/sqrt(N1(tt)+N2(tt));%*sqrt(10/11);
    psi_2 = psi_2/sqrt(N1(tt)+N2(tt));%*sqrt(1/11);
    N1(tt) = sum(sum(abs(psi_1).^2)*h)*h;
    N2(tt) = sum(sum(abs(psi_2).^2)*h)*h;
    
    delta1(tt) = sum(sum(abs(psi_1).^2 .* (xx.^2 + yy.^2)))*h*h;
    delta2(tt) = sum(sum(abs(psi_2).^2 .* (xx.^2 + yy.^2)))*h*h;
    
    figure(1); subplot(1,2,1); surf(y_k,x_j,abs(psi_1).^2); shading interp; view(0,90);
    figure(1); subplot(1,2,2); surf(y_k,x_j,abs(psi_2).^2); shading interp; view(0,90);
end
% 
figure(3); plot((1:tTotal)*delta_t,N1);
hold on; plot((1:tTotal)*delta_t,N2,'b:');
hold on; plot((1:tTotal)*delta_t,N1+N2,'b--');


figure(4); plot((1:tTotal)*delta_t,delta1);
hold on; plot((1:tTotal)*delta_t,delta2,'b:');
hold on; plot((1:tTotal)*delta_t,delta1+delta2,'b--');