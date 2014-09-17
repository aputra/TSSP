%% GFDN for finding ground states of single condensate
% Ref: Wang's thesis -- using BEFD (Backward Euler finite difference).
% Not tested completely, program is slow... I hate running FDTD since the
% matrix is big and there is matrix division which is slow =(

clear all; close all; clc;

Omega = 0.7;

gamma_x = 1;
gamma_y = 1;

delta_t = 0.05;          % time step interval

ax = -6; bx = 6;    % xmin and xmax
cy = -6; dy = 6;    % ymin and ymax

h = 1/6;           % h is the mesh size: h = (bx-ax)/Mx = (dy-cy)/Ny
Mx = round((bx-ax)/h);
Ny = round((dy-cy)/h);

x_j = (ax+h):h:(bx-h);
y_k = (cy+h):h:(dy-h);

j_idx = 0:Mx;
m_idy = 0:Ny;

[yy,xx] = meshgrid(y_k,x_j);

% Initial wavefunction
psi_ho = (gamma_x*gamma_y)^(1/4)/sqrt(pi).*exp(-(gamma_x*xx.^2+gamma_y*yy.^2)/2);                     % xj * yk
psi_v = (gamma_x*gamma_y)^(1/4)/sqrt(pi)*(xx+1i*yy).*exp(-(xx.^2+yy.^2)/2);                     % xj * yk
psi_cc = ((1-Omega)*psi_ho + Omega*psi_v);
% psi_1 = psi_v;
norm_fac = sum(sum(abs(psi_cc).^2))*h*h;

psi_1 = psi_cc/sqrt(norm_fac);

tTotal = round(500/delta_t);
N1 = zeros(1,tTotal);

delta1 = zeros(1,tTotal);

beta2D = 100;
            
V2D_1 = 1/2*(gamma_x^2*xx.^2 + gamma_y^2*yy.^2 ); % xj * yk

% convert matrix to vector
aa = V2D_1.'; Vxy = aa(:);
bb = (-delta_t/2/h^2 - 1i*delta_t*Omega/2/h*y_k);
B_mat = repmat(bb,1,length(x_j)-1);
cc = (-delta_t/2/h^2 + 1i*delta_t*Omega/2/h*y_k);
C_mat = repmat(cc,1,length(x_j)-1);

dd = (-delta_t/2/h^2 + 1i*delta_t*Omega/2/h*x_j);
Drep = repmat(dd.',1,length(y_k)).';
D_mat = Drep(:);
D_mat((1:(Mx-1))*(Ny-1))=0;
ee = (-delta_t/2/h^2 - 1i*delta_t*Omega/2/h*x_j);
Erep = repmat(ee.',1,length(y_k)).';
E_mat = Erep(:);
E_mat((1:(Mx-1))*(Ny-1))=0;
ABC = diag(B_mat,length(y_k)) + diag(C_mat,-length(y_k)) + diag(D_mat(1:(length(D_mat)-1)),1) + diag(E_mat(1:(length(D_mat)-1)),-1);

% aa = (abs(psi_1).^2).'; psi_sqxy = aa(:);
aa = (psi_1).'; psi_xy = aa(:);

for tt = 1:tTotal
    FD_matrix = ABC + diag(1 + 2*delta_t/h^2 + delta_t*Vxy + delta_t*beta2D*abs(psi_xy).^2);
    psi_sqxy = FD_matrix\psi_xy;
%     FD_matrix = FD_matrix + diag(B_mat,length(y_k));
%     FD_matrix = FD_matrix + diag(C_mat,-length(y_k));
%     FD_matrix = FD_matrix + diag(D_mat(1:(length(D_mat)-1)),1);
%     FD_matrix = FD_matrix + diag(E_mat(1:(length(D_mat)-1)),-1);
    % First step
    
    N1(tt) = (sum(abs(psi_sqxy).^2))*h*h;
    psi_xy = psi_sqxy./sqrt(N1(tt));%*sqrt(10/11);
    
    psi_1 = vec2mat(psi_xy,Ny-1);
    N1(tt) = sum(sum(abs(psi_1).^2))*h*h;
    delta1(tt) = sum(sum(abs(psi_1).^2 .* (xx.^2 + yy.^2)))*h*h;
    figure(1); surf(y_k,x_j,abs(psi_1).^2*10); shading interp; view(0,90);
end