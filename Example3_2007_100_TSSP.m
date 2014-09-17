%% TSSP for coupled GPE with applications to rotating BECs
% Ref: Wang, J Comp and App Math 205 (2007) 88-104
% Tested, result matched to Example 1 (Fig 1)

clear all; close all; clc;

Omega = 0.85;
lambda = -1;
load Ex3groundstate.mat;

gamma_x = 1;
gamma_y = 1;

delta_t = 0.002;          % time step interval

h = 1/10;           % h is the mesh size: h = (bx-ax)/Mx = (dy-cy)/Ny
ax = -10+h; bx = 10-h;    % xmin and xmax
cy = -10+h; dy = 10-h;    % ymin and ymax


Mx = round((bx-ax)/h);
Ny = round((dy-cy)/h);

x_j = ax:h:bx;
y_k = cy:h:dy;

% FFT summation indices for x and y spaces:
px = -Mx/2:(Mx/2-1);
qy = -Ny/2:(Ny/2-1);

% j_idx = 0:(Mx-1);
% m_idy = 0:(Ny-1);

% FFT reciprocal space
mu_p = 2*pi*px/(bx-ax);
la_q = 2*pi*qy/(dy-cy);

[yy,xx] = meshgrid(y_k,x_j);

% Initial wavefunction
% psi_1 = 1/sqrt(2*pi)*exp(-(xx.^2+yy.^2)/2);                     % xj * yk

psi_2 = zeros(size(psi_1));

tTotal = round(3.5/delta_t);
N1 = zeros(1,tTotal);
N2 = zeros(1,tTotal);

delta1 = zeros(1,tTotal);
delta2 = zeros(1,tTotal);
Lz = zeros(1,tTotal);

beta2D = 2000 * [1.03    1.0;
                 1.0     0.97];
            
% % V2D_1 = 1/2*(gamma_x^2*xx(1:Mx,1:Ny).^2 + gamma_y^2*yy(1:Mx,1:Ny).^2 ); % xj * yk
% % V2D_2 = 1/2*(gamma_x^2*xx(1:Mx,1:Ny).^2 + gamma_y^2*yy(1:Mx,1:Ny).^2 ); % xj * yk
            
% % V2D_1 = 1/2*(gamma_x^2*xx.^2 + gamma_y^2*yy.^2 ); % xj * yk
% % V2D_2 = 1/2*(gamma_x^2*xx.^2 + gamma_y^2*yy.^2 ); % xj * yk

for tt = 1:tTotal
    ff = xx.*yy;
    gg = 0;
    V2D_1 = 1/2*(xx.^2 + yy.^2) + 10 + 0.499*(ff*cos(20.54*(tt-1)*delta_t) + gg*sin(20.54*(tt-1)*delta_t)) ; % xj * yk
    V2D_2 = 1/2*(xx.^2 + yy.^2) - 10 - 0.499*(ff*cos(20.54*(tt-1)*delta_t) + gg*sin(20.54*(tt-1)*delta_t)) ; % xj * yk
    
    % First step
    psi_hat1 = FFT_x(psi_1,Mx,Ny,x_j,mu_p,ax); % Ny * px
    psi_hat2 = FFT_x(psi_2,Mx,Ny,x_j,mu_p,ax); % Ny * px
    
    psi_1_step = rotatingIFFT_x(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    psi_2_step = rotatingIFFT_x(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    % Second step
    psi_hat1 = FFT_y(psi_1_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy
    psi_hat2 = FFT_y(psi_2_step.',Mx,Ny,y_k,la_q,cy); % Mx * qy
    
    psi_1_step = rotatingIFFT_y(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    psi_2_step = rotatingIFFT_y(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    
    % Third step
    psi_1_stepNew = cos(lambda*delta_t/2)*psi_1_step + 1i*sin(lambda*delta_t/2)*psi_2_step;    % xj * yk
    psi_2_stepNew = 1i*sin(lambda*delta_t/2)*psi_1_step + cos(lambda*delta_t/2)*psi_2_step;    % xj * yk
    
    % Fourth step
    psi_1_step = exp(-1i*delta_t*(V2D_1 + beta2D(1,1)*abs(psi_1_stepNew).^2 + beta2D(1,2)*abs(psi_2_stepNew).^2)).*psi_1_stepNew;   % xj * yk
    psi_2_step = exp(-1i*delta_t*(V2D_2 + beta2D(2,1)*abs(psi_1_stepNew).^2 + beta2D(2,2)*abs(psi_2_stepNew).^2)).*psi_2_stepNew;   % xj * yk
    
    % Fifth step
    psi_1_stepNew = cos(lambda*delta_t/2)*psi_1_step + 1i*sin(lambda*delta_t/2)*psi_2_step;    % xj * yk
    psi_2_stepNew = 1i*sin(lambda*delta_t/2)*psi_1_step + cos(lambda*delta_t/2)*psi_2_step;    % xj * yk
    
    % Sixth step
    psi_hat1 = FFT_y(psi_1_stepNew,Mx,Ny,y_k,la_q,cy); % Mx * qy
    psi_hat2 = FFT_y(psi_2_stepNew,Mx,Ny,y_k,la_q,cy); % Mx * qy
    
    psi_1_step = rotatingIFFT_y(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    psi_2_step = rotatingIFFT_y(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,la_q,cy); % xj * yk
    
    % Seventh step --> advancing from t_n to t_n+1 --> end of one loop
    psi_hat1 = FFT_x(psi_1_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
    psi_hat2 = FFT_x(psi_2_step,Mx,Ny,x_j,mu_p,ax); % Ny * px
    
    psi_1_step = rotatingIFFT_x(Omega,delta_t,psi_hat1,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    psi_2_step = rotatingIFFT_x(Omega,delta_t,psi_hat2,Mx,Ny,x_j,y_k,mu_p,ax); % yk * xj
    
    psi_1 = psi_1_step.';
    psi_2 = psi_2_step.';
    
    [gy1,gx1] = gradient(psi_1,h);
    [gy2,gx2] = gradient(psi_2,h);
    
    L1 = sum(sum(conj(psi_1).*(-1i*xx.*gy1+1i*yy.*gx1)))*h*h;
    L2 = sum(sum(conj(psi_2).*(-1i*xx.*gy2+1i*yy.*gx2)))*h*h;
    Lz(tt) = L1 + L2;
    
    N1(tt) = sum(sum(abs(psi_1).^2))*h*h;
    N2(tt) = sum(sum(abs(psi_2).^2))*h*h;
    
    delta1(tt) = sum(sum(abs(psi_1).^2 .* (xx.^2 + yy.^2)))*h*h;
    delta2(tt) = sum(sum(abs(psi_2).^2 .* (xx.^2 + yy.^2)))*h*h;
    
    figure(1); subplot(1,2,1); surf(y_k,x_j,abs(psi_1).^2); shading interp; view(0,90);
    figure(1); subplot(1,2,2); surf(y_k,x_j,abs(psi_2).^2); shading interp; view(-30,30); title(tt*delta_t);
end

figure(3); plot((1:tTotal)*delta_t,N1);
hold on; plot((1:tTotal)*delta_t,N2,'b:');
hold on; plot((1:tTotal)*delta_t,N1+N2,'b--');
hold on; plot((1:tTotal)*delta_t,Lz,'go');

figure(4); plot((1:tTotal)*delta_t,delta1);
hold on; plot((1:tTotal)*delta_t,delta2,'b:');
hold on; plot((1:tTotal)*delta_t,delta1+delta2,'b--');