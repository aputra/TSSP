function [psi] = rotatingIFFT_x_mod(Omega,dt,psi_hat,Mx,Ny,x,y,mu,a)
    % matrix size of psi_hat: Ny * px
    % matrix size of psi: yk * xj
    
    PDEmat = exp(-dt/4* (repmat(mu.^2,Ny-1,1)+2*Omega*y.'*mu) );
    psi = (PDEmat.*psi_hat) * exp(1i* mu.'*(x-a));

end