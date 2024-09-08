function fval = KM_theory_S_K_with_surface_reflection_correction(X,nh,nsurr,L,R,T)%,Qsca)

% %Define the variables

S = X(1);
K = X(2);


%%
% surface reflectance correction
if nsurr==nh
    rc = 0;
    ri = 0;
else
    rc = ((real(nh)-real(nsurr)).^2+(imag(nh)-imag(nsurr)).^2)./((real(nh)+real(nsurr)).^2+(imag(nh)+imag(nsurr)).^2);
 
      n = nsurr./nh;
    
    
        fun = @(x) (sin(x)).*(cos(x)).*((abs((((((n.*n) - ((sin(x)).^2)).^0.5) - cos(x))./((((n.*n) - ((sin(x)).^2)).^0.5) + cos(x)))).^2) + ((abs((((n.*n).*cos(x)) - (((n.*n) - ((sin(x)).^2)).^0.5))./(((n.*n).*cos(x))+(((n.*n) - ((sin(x)).^2)).^0.5)))).^2));
        ri = integral(fun,0,(pi/2));
end   
   

%%

%Here R and T are the reflectance and transmittance values from Lumerical
%when air/matrix interface reflection is present

Rkm = (R-rc)./(1-rc-ri*(1-R));
Tkm = T*(1-ri*Rkm)./(1-rc);

fval(1,1) = Rkm - ((1-(real(sqrt(K./(K+2*S))))).*(1+(real(sqrt(K./(K+2*S))))).*(exp((real(sqrt(K.*(K+2*S))))*L)-exp(-(real(sqrt(K.*(K+2*S))))*L)))./((1+(real(sqrt(K./(K+2*S))))).^2.*exp((real(sqrt(K.*(K+2*S))))*L)-(1-(real(sqrt(K./(K+2*S))))).^2.*exp(-(real(sqrt(K.*(K+2*S))))*L));
fval(1,2) = Tkm - 4*(real(sqrt(K./(K+2*S))))./((1+(real(sqrt(K./(K+2*S))))).^2.*exp((real(sqrt(K.*(K+2*S))))*L)-(1-(real(sqrt(K./(K+2*S))))).^2.*exp(-(real(sqrt(K.*(K+2*S))))*L));

