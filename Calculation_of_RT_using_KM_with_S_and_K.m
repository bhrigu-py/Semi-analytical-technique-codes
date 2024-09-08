clc;
clear all;
close all


%%
%input data
data=xlsread('C:\Users\ADMIN\OneDrive - Indian Institute of Technology Bombay\Desktop\Experimental paper TiO2-PDMS\For paper\FDTD_SM_r0.2_f0.1_0.15_0.2.xlsx','r0.25');%('G:\data_files_without_substrate\RT comparison using new derived formula for coatings without substrate.xlsx','absrbng_diffuse');
lambda=data(:,11);

S = data(:,13);
K = data(:,14);

S = abs(S);
K = abs(K);

data1=xlsread('C:\Users\ADMIN\OneDrive - Indian Institute of Technology Bombay\Desktop\PhD topic\Mie code\Data_nk\PDMS_main_file.xlsx');
lmbda1=data1(:,1);
mh1=data1(:,2);
mh2=data1(:,3);
nh1=interp1(lmbda1,mh1,lambda);
nh2=interp1(lmbda1,mh2,lambda);
nh = nh1+1i*nh2; %matrix index
nsurr = 1.0; %index of surrounding
L=300; %thickness of the film
nb = 1.0;

%%
%diffuse R and T
for i=1:length(K)
    if K(i)>0
        gama_g(i)=real(sqrt(K(i)./(K(i)+2*S(i))));
        A_g(i)=real(sqrt(K(i).*(K(i)+2*S(i))));
        if A_g(i)*L<10^5
            R1(i)=real(((1-gama_g(i)).*(1+gama_g(i)).*(exp(A_g(i)*L)-exp(-A_g(i)*L)))./((1+gama_g(i)).^2.*exp(A_g(i)*L)-(1-gama_g(i)).^2.*exp(-A_g(i)*L)));
            T1(i)=4*gama_g(i)./((1+gama_g(i)).^2.*exp(A_g(i)*L)-(1-gama_g(i)).^2.*exp(-A_g(i)*L));
        else
             R1(i) = real((K(i)+S(i)-sqrt(K(i).^2+2*K(i).*S(i)))./S(i));
             T1(i) = 0;
        end
    else
        R1(i) = S(i)*L/(S(i)*L+1);
        T1(i) = 1-R1(i);
    end
end

% a = 1+K./S;
% b = (a.^2-1).^(1/2);
% 
% R1 = 1./(a+b.*coth(b.*S*L));
% T1 = b./(a.*sinh(b.*S*L)+b.*cosh(b.*S*L));
       
% if A_g*L>10^5
%    R1 = real((K+S-sqrt(K.^2+2*K.*S))./S);
% end
     
    
   %% surface reflectance correction

   for i=1:length(nh)
    if nsurr==nh%nh(i)
        rc(i) = 0;
        ri(i) = 0;
    else
        rc(i) = ((real(nh(i))-real(nsurr)).^2+(imag(nh(i))-imag(nsurr)).^2)./((real(nh(i))+real(nsurr)).^2+(imag(nh(i))+imag(nsurr)).^2);
      
            n(i) = nsurr./nh(i);
        
            fun = @(x) (sin(x)).*(cos(x)).*((abs((((((n(i).*n(i)) - ((sin(x)).^2)).^0.5) - cos(x))./((((n(i).*n(i)) - ((sin(x)).^2)).^0.5) + cos(x)))).^2) + ((abs((((n(i).*n(i)).*cos(x)) - (((n(i).*n(i)) - ((sin(x)).^2)).^0.5))./(((n(i).*n(i)).*cos(x))+(((n(i).*n(i)) - ((sin(x)).^2)).^0.5)))).^2));
            ri(i) = integral(fun,0,(pi/2));

    end
    
 %introducing surface reflection correction
    R_KM(i) = rc(i)+(1-rc(i)).*(1-ri(i)).*R1(i)./(1-ri(i).*R1(i));
    T_KM(i) = (1-rc(i)).*T1(i)./(1-ri(i).*R1(i));
   end

% 
%     if nsurr==nh
%         rc = 0;
%         ri = 0;
%     else
%         rc = ((real(nh)-real(nsurr)).^2+(imag(nh)-imag(nsurr)).^2)./((real(nh)+real(nsurr)).^2+(imag(nh)+imag(nsurr)).^2);
%       
%             n = nsurr./nh;
%         
%             fun = @(x) (sin(x)).*(cos(x)).*((abs((((((n.*n) - ((sin(x)).^2)).^0.5) - cos(x))./((((n.*n) - ((sin(x)).^2)).^0.5) + cos(x)))).^2) + ((abs((((n.*n).*cos(x)) - (((n.*n) - ((sin(x)).^2)).^0.5))./(((n.*n).*cos(x))+(((n.*n) - ((sin(x)).^2)).^0.5)))).^2));
%             ri(i) = integral(fun,0,(pi/2));
% 
%     end
%     
%  %introducing surface reflection correction
%     R_KM = rc+(1-rc).*(1-ri).*R1./(1-ri.*R1);
%     T_KM = (1-rc).*T1./(1-ri.*R1);
  
   R_KM = R_KM';
   T_KM = T_KM';
%% Reflection and transmission in the presence of substrate with reflectance Rb
%Kubelka 1948
%Rb is calculated using fresnel reflection that is reflectance at the
%interface of matrix and substrate. If no substrate then Rb is simply
%reflectance at matrix and air.
%Rb = ((real(nh)-nb).^2+(imag(nh)).^2)./((real(nh)+nb).^2+(imag(nh)).^2);
%Rb = ((real(nh(i))-real(nb)).^2+(imag(nh(i))-imag(nb)).^2)./((real(nh(i))+real(nb)).^2+(imag(nh(i))+imag(nb)).^2);
Rb = ((real(nh)-real(nb)).^2+(imag(nh)-imag(nb)).^2)./((real(nh)+real(nb)).^2+(imag(nh)+imag(nb)).^2);
R = R_KM + (T_KM.^2.*Rb)./(1-R_KM.*Rb);
%T = T_KM;
if Rb==0
    T = T_KM;
else
    T = ((R-R_KM).*(1./Rb - R_KM)).^0.5;
end
A = 1-R-T;

N = [R, T, A];
    


    