function [psi_hat] = FFT_y_sine(psi,Mx,Ny,y,la,c)
    % matrix size of psi : Mx * Ny
    % matrix size of psi_hat: Mx * qy
    
    psi_hat = 1/Ny * psi * exp(-1i*(y-c).' * la);
%     psi_hat = 1/Ny * psi(1:Mx,1:Ny) * exp(-1i* (y(1:Ny)-c).' * la);
end