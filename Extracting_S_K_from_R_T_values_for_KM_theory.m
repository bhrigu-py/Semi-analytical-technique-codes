%Solving KM equations to caluclate S and K using reflectivity and
%transmissivity data
clc;
close all;
clear all

data=xlsread('C:\Users\ADMIN\OneDrive - Indian Institute of Technology Bombay\Desktop\Experimental paper TiO2-PDMS\For paper\FDTD_SM_r0.2_f0.1_0.15_0.2.xlsx','newSigmaTiO2_PDMS_poly_f15');%('G:\data_files_without_substrate\RT comparison using new derived formula for coatings without substrate.xlsx','absrbng_diffuse');
lambda=data(1:200,6);
R_total = data(1:200,7);
T_total = data(1:200,8);

%%


data1=xlsread('C:\Users\ADMIN\OneDrive - Indian Institute of Technology Bombay\Desktop\PhD topic\Mie code\Data_nk\PDMS_main_file.xlsx');
lmbda1=data1(:,1);
mh1=data1(:,2);
mh2=data1(:,3);
nh1=interp1(lmbda1,mh1,lambda);
nh2=interp1(lmbda1,mh2,lambda);
nh = nh1+1i*nh2; %matrix index
nsurr = 1.4; %index of surrounding
L=10; %thickness of the film
nb = 1.4;%refractive index of substrate, this code works when nb = nh, if nh is imaginary, nb = real(nh). that is Rb should not be high.
%If Rb is high (when imag(nh) is very large, then substrate correction in reflection obtained from Lumerical is required before using this code. Or this correction can be added in the
%function KM_theory_S_K_with_surface_reflection_correction(x,nh,nsurr,L,R_total(i),T_total(i)),x0)

%%
%with surface correction
x0 = [0.001;0.001]; %initial guess

for i = 1:length(R_total)
    if (1 - (R_total(i)+T_total(i)))<0.04
    S(i) = R_total(i)/(L*(1-R_total(i)));
    K(i) = 0;
    else
    xsolve{i} =fsolve(@(x) KM_theory_S_K_with_surface_reflection_correction(x,nh(i),nsurr,L,R_total(i),T_total(i)),x0);
    S(i) = xsolve{i}(1);
    K(i) = xsolve{i}(2);
    end
end
%%
% for j = 1:length(R_total)
% S(j) = xsolve{j}(1);
% end
S = abs(S');
% for j = 1:length(R_total)
% K(j) = xsolve{j}(2);
% end
K = abs(K');

figure
plot(lambda,S)
hold on
plot(lambda,K)

xlabel('Wavelength (nm)');
legend('S','K')