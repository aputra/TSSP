%% TSSP for coupled GPE with applications to rotating BECs
% Ref: Wang, J Comp and App Math 205 (2007) 88

clear all; close all; clc;

Omega = 0.8;
lambda = 0;

% % % gamma_x = 1;
% % % gamma_y = 1.5;

delta_t = 0.01;          % time step interval

ax = -8; bx = 8;    % xmin and xmax
cy = -8; dy = 8;    % ymin and ymax

h = 1/8;           % h is the mesh size: h = (bx-ax)/Mx = (dy-cy)/Ny
Mx = round((bx-ax)/h);
Ny = round((dy-cy)/h);

x_j = ax:h:bx;
y_k = cy:h:dy;

% FFT summation indices for x and y spaces:
px = -Mx/2:(Mx/2-1);
qy = -Ny/2:(Ny/2-1);

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
% psi_1 = psi_1/sqrt(sum(sum(abs(psi_1).^2))*h*h);
% % mu_chem = sqrt(1622.9/pi);
% % beta = 1622.9;
% % psi_1 = sqrt((mu_chem - 1/2*xx.^2 - 1/2*yy.^2)/beta);
% % psi_1(find(imag(psi_1)))=0;

psi_2 = zeros(size(xx));


tTotal = round(40/delta_t);
N1 = zeros(1,tTotal);
N2 = zeros(1,tTotal);

delta1 = zeros(1,tTotal);
delta2 = zeros(1,tTotal);

beta2D = 100 * [1    0;
                0    1];
V2D_1 = 1/2*(xx.^2 + yy.^2); % xj * yk
V2D_2 = 1/2*(xx.^2 + yy.^2); % xj * yk            
    
for tt = 1:tTotal
% %     V2D_1 = 1/2*(xx.^2 + yy.^2 + 10 + 0.499*(ff*cos(20.54*(tt-1)*delta_t) + gg*sin(20.54*(tt-1)*delta_t)) ); % xj * yk
% %     V2D_2 = 1/2*(xx.^2 + yy.^2 - 10 - 0.499*(ff*cos(20.54*(tt-1)*delta_t) + gg*sin(20.54*(tt-1)*delta_t)) ); % xj * yk
    
    % First step
    psi_hat1 = FFT_x_sine(psi_1,Mx,Ny,x_j,mu_p,ax); % Ny * px
%     psi_hat2 = FFT_x(psi_2,Mx,Ny,x_j,mu_p,ax); % Ny * px
    
    psi_1_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
%     psi_2_step = rotatingIFFT_x(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    % Second step
    psi_hat1 = FFT_y_sine(psi_1_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy
%     psi_hat2 = FFT_y(psi_2_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy
    
    psi_1_stepNew = rotatingIFFT_y_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
%     psi_2_step = rotatingIFFT_y(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    
    % Third step
%     psi_1_stepNew = cos(lambda*delta_t/2)*psi_1_step;% + 1i*sin(lambda*delta_t/2)*psi_2_step;    % xj * yk
%     psi_2_stepNew = 1i*sin(lambda*delta_t/2)*psi_1_step + cos(lambda*delta_t/2)*psi_2_step;    % xj * yk
    
    % Fourth step
    psi_1_stepNew = exp(-delta_t*(V2D_1 + beta2D(1,1)*abs(psi_1_stepNew).^2)).*psi_1_stepNew;   % xj * yk
%     psi_2_step = exp(-1i*delta_t*(V2D_2 + beta2D(2,1)*abs(psi_1_stepNew).^2 + beta2D(2,2)*abs(psi_2_stepNew).^2)).*psi_2_stepNew;   % xj * yk
    
    % Fifth step
%     psi_1_stepNew = cos(lambda*delta_t/2)*psi_1_step;% + 1i*sin(lambda*delta_t/2)*psi_2_step;    % xj * yk
%     psi_2_stepNew = 1i*sin(lambda*delta_t/2)*psi_1_step + cos(lambda*delta_t/2)*psi_2_step;    % xj * yk
    
    % Sixth step
    psi_hat1 = FFT_y_sine(psi_1_stepNew,Mx,Ny,y_k,la_q,cy); % Mx * qy
%     psi_hat2 = FFT_y(psi_2_stepNew,Mx,Ny,y_k,la_q,cy); % Mx * qy
    
    psi_1_step = rotatingIFFT_y_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
%     psi_2_step = rotatingIFFT_y(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    
    % Seventh step --> advancing from t_n to t_n+1 --> end of one loop
    psi_hat1 = FFT_x_sine(psi_1_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
%     psi_hat2 = FFT_x(psi_2_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
    
    psi_1_step = rotatingIFFT_x_sine(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
%     psi_2_step = rotatingIFFT_x(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
%     step_diff = sum(sum(abs(psi_1_step.'- psi_1).^2))*h*h;
    
    psi_1 = psi_1_step.';
%     psi_2 = psi_2_step.';
    
    N1(tt) = sum(sum(abs(psi_1).^2))*h*h;
%     N2(tt) = sum(sum(abs(psi_2).^2))*h*h;
    psi_1 = psi_1/sqrt(N1(tt));
    delta1(tt) = sum(sum(abs(psi_1).^2 .* (xx.^2 + yy.^2)))*h*h;
%     delta2(tt) = sum(sum(abs(psi_2).^2 .* (xx.^2 + yy.^2)))*h*h;
    
%     figure(1); hold on; plot(tt*delta_t, step_diff, 'bx');
    figure(1); surf(y_k,x_j,abs(psi_1).^2); shading interp; view(0,90);
end

figure(3); plot((1:tTotal)*delta_t,N1);
% hold on; plot((1:tTotal)*delta_t,N2,'b:');
% hold on; plot((1:tTotal)*delta_t,N1+N2,'b--');


figure(4); plot((1:tTotal)*delta_t,delta1);
% hold on; plot((1:tTotal)*delta_t,delta2,'b:');
% hold on; plot((1:tTotal)*delta_t,delta1+delta2,'b--');
% % 
% % 
% %     figure(1); subplot(1,2,1); surf(y_k,x_j,abs(psi_1).^2); shading interp; view(0,90); caxis([0 0.01])
% %     figure(1); subplot(1,2,2); surf(y_k,x_j,abs(psi_2).^2); shading interp; view(0,90); caxis([0 0.005])