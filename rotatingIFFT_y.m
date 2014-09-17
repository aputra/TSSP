function [psi] = rotatingIFFT_y(Omega,dt,psi_hat,Mx,Ny,x,y,la,c)
    % matrix size of psi_hat: Mx * qy
    % matrix size of psi: xj * yk
    
    PDEmat = exp(-1i*dt/4* (repmat(la.^2,Mx+1,1)-2*Omega*x.'*la) );
    psi = (PDEmat.*psi_hat) * exp(1i* la.'*(y-c));

end