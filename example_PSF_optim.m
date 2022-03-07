% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen
% author: Yiming Li & Shuang Fu
% email: liym2019@sustech.edu.cn & fus2020@mail.sustech.edu.cn
% date: 2022.3.2
% Tested with Matlab 2019a
%% This example first optimizes a pupil using 21 zernike modes,then projects it on the DM 
% and continues to optimize using DM influence functions
clear;close all;clc;
addpath('zernike','chirp_z','optimize_psf_crlb_zernike','opt_realDM_psf_crlb','utils_funcs')
%% define the parameters used
paraSim.NA = 1.35;                      % numerical aperture of obj             
paraSim.refmed = 1.405;                 % refractive index of sample medium
paraSim.refcov = 1.524;                 % refractive index of converslip
paraSim.refimm = 1.405;                 % refractive index of immersion oil
paraSim.lambda = 660;                   % wavelength of emission
paraSim.objStage0 = -600;              % nm, distance from coverslip to nominal focal plane, when RI match,this is useless
paraSim.zemit0 = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0);% reference emitter z position, nm, distance of molecule to coverslip
paraSim. pixelSizeX = 108;              % nm, pixel size of the image
paraSim. pixelSizeY = 108;              % nm, pixel size of the image
paraSim.Npupil = 64;                    % sampling at the pupil plane
paraSim.show = false;                   % plot the CRLB,pupil,psf model
paraSim.alpha = 0;                      % the coefficient of regularization term to confine the spatial extent of the PSF
% total 21 zernike aberration orders, n,m,amplitude(nm)
paraSim.aberrations = [2,-2,0.0; 2, 2,0.0; 3,-1,0.0; 3,1,0.0; 4, 0,0.0;...
                       3,-3,0.0; 3, 3,0.0; 4,-2,0.0; 4,2,0.0; 5,-1,0.0;...
                       5, 1,0.0; 6, 0,0.0; 4,-4,0.0; 4,4,0.0; 5,-3,0.0;...
                       5, 3,0.0; 6,-2,0.0; 6, 2,0.0; 7,1,0.0; 7,-1,0.0;...
                       8, 0,0.0];
% parameters to compute the CRLB and psf model
Nmol = 31;                              % number of discrete z positions 
Npixels = 71;                           % number of pixels of the image
Nphotons = 2000 +0*10000*rand(1,Nmol);  % photons for each PSF
bg = 20 +0*10*rand(1,Nmol);             % background for each PSF
paraSim.Nphotons = Nphotons;
paraSim.bg = bg;
paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;
paraSim.xemit = ones(1,Nmol)*0;                %nm x positions of each PSF
paraSim.yemit = ones(1,Nmol)*0;                %nm y positions of each PSF
paraSim.zemit = linspace(-3000,3000,Nmol)*1;     %nm z positions of each PSF
paraSim.objStage = linspace(-1000,1000,Nmol)*0;%nm objective positions of each PSF
%% define the initial pupil using Zernike modes
aberrations_start = [2, 2,0.0; 2, 2,60.0; 3,-1,0.0; 3,1,0.0; 4, 0,0.0;...
                     3,-3,0.0; 3, 3,0.0; 4,-2,0.0; 4,2,0.0; 5,-1,0.0;...
                     5, 1,0.0; 6, 0,0.0; 4,-4,0.0; 4,4,0.0; 5,-3,0.0;...
                     5, 3,0.0; 6,-2,0.0; 6, 2,0.0; 7,1,0.0; 7,-1,0.0;...
                     8, 0,0.0];
% aberrations_start(:,3) = 200*(2*rand([21,1])-1);

% compute and show the crlb, PSFs, pupil of these Zernike coefficients
paraSim.show = true;disp('Zernike start crlb')
[y0,grad0] = sum_crlb_at_z(aberrations_start(:,3),paraSim);
paraSim.show = false;
pupil_start = Zernike_construct_pupil(aberrations_start, paraSim);
figure;imagesc(pupil_start);title('Zernike start pupil')
%% optimize Zernike coefficients to minimize average CRLB_3D (take miniutes to hours)
f = @(zernike_coefs)sum_crlb_at_z(zernike_coefs,paraSim);
% fminicon
lb = -paraSim.lambda*ones(21,1);
ub = paraSim.lambda*ones(21,1);
options = optimoptions(@fmincon,'algorithm','interior-point','display','iter-detailed',...
'diagnostics','on','SpecifyObjectiveGradient',true,'UseParallel',true,...
'MaxIterations',400,'OptimalityTolerance',1e-6);
[optimized_coefs,fval,exitflag,output,lambda,grad,hess]=...
    fmincon(f, aberrations_start(:,3),[],[],[],[],lb,ub,[],options);
aberrations_end =aberrations_start;
aberrations_end(:,3)=optimized_coefs;

% compute and show the crlb, PSFs, pupil of Zernike optimization
paraSim.show = true;disp('Zernike final crlb')
[y0,grad0] = sum_crlb_at_z(aberrations_end(:,3),paraSim);
paraSim.show = false;
pupil_end = Zernike_construct_pupil(aberrations_end, paraSim);
figure;imagesc(pupil_end);title('Zernike final pupil')
%% load DM influence functions, project the Zernike pupil on DM
% (xxx1875 corresponds to 1.5 NA, xxx1687.5 corresponds to 1.35NA)
paraSim.dm_voltage = zeros(140,1);
load('./data/17BW023#111_20210609_170825DM_influence_function_and_u_flat1687.5.mat');
smooth=true;
paraSim.dm_basis = dm_basis_resize(DM_impulse,paraSim.Npupil,smooth);

% look the DM basis
tmp=zeros(paraSim.Npupil,paraSim.Npupil,140);
for i=1:140
tmp(:,:,i)=paraSim.dm_basis(i,:,:);
end
imageslicer(tmp);

% project the optimized zernike pupil on the DM
pupil = Zernike_construct_pupil(aberrations_end, paraSim);
dm_coefs_ls = realDMbasis_decompose_ls(pupil,paraSim.dm_basis);
reconstruction = DM_construct_pupil(dm_coefs_ls, paraSim);
figure;
subplot(1,3,1);
imagesc(pupil);axis square;title('Zernike pupil');colorbar
subplot(1,3,2);
imagesc(reconstruction);axis square;title('DM pupil');colorbar
subplot(1,3,3);
imagesc(reconstruction-pupil);axis square;title('error');colorbar
disp(['DM charges for each actuator']);
disp(dm_coefs_ls');
disp(['rse: ',num2str( 100* norm(reconstruction-pupil,2)/norm(pupil,2)),'%'] )
%% compute and show the crlb, PSFs with the projected pupil
dm_voltage_start = dm_coefs_ls;

% confine the voltage between the applicable range [-1,1]
dm_voltage_start(dm_voltage_start<-1)=-0.99999;
dm_voltage_start(dm_voltage_start>1)=0.99999;

paraSim.show = true;disp('DM start crlb')
[y0,grad0] = realDM_sum_crlb_at_z(dm_voltage_start,paraSim);
paraSim.show = false;
pupil_start = DM_construct_pupil(dm_voltage_start, paraSim);
figure;imagesc(pupil_start);title('DM start pupil')
%% optimize the voltages of DM actuators to minimize average CRLB_3D (take hours to days)
f = @(dm_voltage)realDM_sum_crlb_at_z(dm_voltage,paraSim);
lb = -1*ones(140,1);
ub = 1*ones(140,1);
fminconoptions = optimoptions(@fmincon,'algorithm','interior-point',...
    'display','iter-detailed','diagnostics','on','UseParallel',true,...
    'SpecifyObjectiveGradient',true,'CheckGradients',false,...
    'MaxIterations',50,'OptimalityTolerance',1e-2);
[optimized_coefs,fval,exitflag,output,lambda,grad,hess]=...
    fmincon(f, dm_voltage_start,[],[],[],[],lb,ub,[],fminconoptions);
dm_voltage_end = optimized_coefs;

% compute and show the crlb, PSFs, pupil phase with the optimized DM voltages
paraSim.show = true;disp('DM final crlb');
[y0,grad0] = realDM_sum_crlb_at_z(dm_voltage_end,paraSim);
paraSim.show = false;
pupil_end = DM_construct_pupil(dm_voltage_end, paraSim);
figure;imagesc(pupil_end);title('DM final pupil')
%% output the control voltage(0,1) loaded on the DM
% Membrane DMs, like the Boston Multi-DM, exhibit a quadratic response of 
% the displacement of one actuator with respect to the applied voltage.
% To have a better linear response you can define a new control variable.
% The relationship here is u = 2*(v.^2)-1
v_flat = load('./data/17BW023#111_FLAT_MAP_COMMAND.txt');
u_flat = 2*(v_flat.^2)-1;
u_final = dm_voltage_end+u_flat;
u_final( u_final<-1 ) = -1;
u_final( u_final>1 ) = 1;

v_final = sqrt((u_final+1)/2);

fileID = fopen(['./output/',filesep,date,'-DMcontrol.txt'],'w');
fprintf(fileID,'%12.8f\n',v_final);
fclose(fileID);
