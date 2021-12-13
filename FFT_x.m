function [psi_hat] = FFT_x(psi,Mx,Ny,x,mu,a)
    % matrix size of psi : Mx * Ny
    % matrix size of psi_hat: Ny * px
    
    psi_hat = 1/Mx * psi(1:Mx,:).' * exp(-1i* (x(1:Mx)-a).' * mu);    
%     psi_hat = 1/Mx * psi(1:Mx,1:Ny).' * exp(-1i* (x(1:Mx)-a).' * mu);    

end

% TESTED CHANGES