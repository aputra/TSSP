function [psi_hat] = FFT_x_sine(psi,Mx,Ny,x,mu,a)
    % matrix size of psi : Mx * Ny
    % matrix size of psi_hat: Ny * px
    
    psi_hat = 1/Mx * psi.' * exp(-1i* (x-a).' * mu);    
%     psi_hat = 1/Mx * psi(1:Mx,1:Ny).' * exp(-1i* (x(1:Mx)-a).' * mu);    
end