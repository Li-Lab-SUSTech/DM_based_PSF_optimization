function [loss,dlossdzer,hessian] =  realDM_sum_crlb_at_z(dm_voltage,parameters)
parameters.dm_voltage=dm_voltage;

% common term for calculate PSF and PSF derivatives
NA = parameters.NA;
refmed = parameters.refmed;
refcov = parameters.refcov;
refimm = parameters.refimm;
lambda = parameters.lambda;
Npupil = parameters.Npupil;
sizeX = parameters.sizeX;
sizeY = parameters.sizeY;
pixelSizeX = parameters.pixelSizeX;
pixelSizeY = parameters.pixelSizeY;
xrange = pixelSizeX*sizeX/2;
yrange = pixelSizeY*sizeY/2;

% pupil radius (in diffraction units) and pupil coordinate sampling
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

% calculation of relevant Fresnel-coefficients for the interfaces
% between the medium and the cover slip and between the cover slip
% and the immersion fluid
CosThetaMed = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
CosThetaCov = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refcov^2);
CosThetaImm = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refimm^2);
%need check again, compared to original equation, FresnelPmedcov is
%multipiled by refmed
FresnelPmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaCov+refcov*CosThetaMed);
FresnelSmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaMed+refcov*CosThetaCov);
FresnelPcovimm = 2*refcov*CosThetaCov./(refcov*CosThetaImm+refimm*CosThetaCov);
FresnelScovimm = 2*refcov*CosThetaCov./(refcov*CosThetaCov+refimm*CosThetaImm);
FresnelP = FresnelPmedcov.*FresnelPcovimm;
FresnelS = FresnelSmedcov.*FresnelScovimm;

% Apoidization for sine condition
% apoid = sqrt(CosThetaImm)./CosThetaMed;
apoid = 1 ./ sqrt(CosThetaMed);
% definition aperture
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
Amplitude = ApertureMask.*apoid;


% setting of vectorial functions
Phi = atan2(YPupil,XPupil);
CosPhi = cos(Phi);
SinPhi = sin(Phi);
CosTheta = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
SinTheta = sqrt(1-CosTheta.^2);

pvec{1} = FresnelP.*CosTheta.*CosPhi;
pvec{2} = FresnelP.*CosTheta.*SinPhi;
pvec{3} = -FresnelP.*SinTheta;
svec{1} = -FresnelS.*SinPhi;
svec{2} = FresnelS.*CosPhi;
svec{3} = 0;

PolarizationVector = cell(2,3);
for jtel = 1:3
  PolarizationVector{1,jtel} = CosPhi.*pvec{jtel}-SinPhi.*svec{jtel};
  PolarizationVector{2,jtel} = SinPhi.*pvec{jtel}+CosPhi.*svec{jtel};
end

wavevector = cell(1,3);
wavevector{1} = (2*pi*NA/lambda)*XPupil;
wavevector{2} = (2*pi*NA/lambda)*YPupil;
wavevector{3} = (2*pi*refimm/lambda)*CosThetaImm;
wavevectorzmed = (2*pi*refmed/lambda)*CosThetaMed; 

%%
ImageSizex = xrange*NA/lambda;
ImageSizey = yrange*NA/lambda;

[Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,sizeX);
[Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,sizeY);

parameters.XPupil = XPupil;
parameters.PolarizationVector = PolarizationVector;
parameters.Amplitude = Amplitude;
parameters.wavevector = wavevector;
parameters.wavevectorzmed = wavevectorzmed;

parameters.Ax = Ax;
parameters.Bx = Bx;
parameters.Dx = Dx;
parameters.Ay = Ay;
parameters.By = By;
parameters.Dy = Dy;

% calculate intensity normalization function using the PSFs at focus
% position without any aberration
FieldMatrix = cell(2,3);
for itel = 1:2
  for jtel = 1:3
    PupilFunction = Amplitude.*PolarizationVector{itel,jtel};
    IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
    FieldMatrix{itel,jtel} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
  end
end

intFocus = zeros(sizeX,sizeY);
for jtel = 1:3
    for itel = 1:2
        intFocus = intFocus + (1/3)*abs(FieldMatrix{itel,jtel}).^2;
    end
end
normIntensity = sum(intFocus(:));
parameters.normIntensity = normIntensity;

%% build theta_GT
Nz= parameters.Nmol;
num_dmpixels = size(parameters.dm_voltage,1);
shared = [ones(1,num_dmpixels) 0 0 0 0 0];  % 1 is shared parameters between z slices, 0 is free parameters between z slices, only consider  [x, y, z, I, bg]
sumShared = sum(shared);
numparams = (num_dmpixels+5) * Nz - sumShared * (Nz-1);
parameters.numparams = numparams;
parameters.num_dmpixels = num_dmpixels;
thetainit_GT = zeros(numparams,1);

bg0_GT = zeros(1,Nz);
Nph_GT = zeros(1,Nz);
x0_GT = zeros(1,Nz);
y0_GT = zeros(1,Nz);
z0_GT = zeros(1,Nz);

% center of mass with nm unit
ImageSizex = parameters.pixelSizeX*parameters.sizeX/2;
ImageSizey = parameters.pixelSizeY*parameters.sizeY/2;

DxImage = 2*ImageSizex/parameters.sizeX;
DyImage = 2*ImageSizey/parameters.sizeY;
ximagelin = -ImageSizex+DxImage/2:DxImage:ImageSizex;
yimagelin = -ImageSizey+DyImage/2:DyImage:ImageSizey;
[YImage,XImage] = meshgrid(yimagelin,ximagelin);
for i = 1:Nz
    bg0_GT(i) = parameters.bg(i);
    Nph_GT(i) = parameters.Nphotons(i);
end
x0_GT = parameters.xemit;
y0_GT = parameters.yemit;
z0_GT = parameters.zemit;


allTheta_GT = zeros(num_dmpixels+5,Nz);
allTheta_GT(num_dmpixels+1,:)=x0_GT';
allTheta_GT(num_dmpixels+2,:)=y0_GT';
allTheta_GT(num_dmpixels+3,:)=z0_GT';
allTheta_GT(num_dmpixels+4,:)=Nph_GT';
allTheta_GT(num_dmpixels+5,:)=bg0_GT';
allTheta_GT(1:num_dmpixels,:) = repmat(parameters.dm_voltage,[1 Nz]);

% for
map = zeros(numparams,3);
n=1;
for i = 1:num_dmpixels+5
    if shared(i)==1
        map(n,1)= 1;
        map(n,2)=i;
        map(n,3)=0;
        n = n+1;
    elseif shared(i)==0
        for j = 1:Nz
            map(n,1)=0;
            map(n,2)=i;
            map(n,3)=j;
            n = n+1;
        end
    end
end

parameters.map = map;
parameters.zemitStack = zeros(Nz,1)'; % move emitter
parameters.objStageStack = zeros(Nz,1)'; %move objStage
parameters.ztype = 'emitter';
parameters.map = map;
parameters.sizeZ = length(parameters.zemit);

for i = 1:numparams
    if map(i,1)==1
        thetainit_GT(i)= mean(allTheta_GT(map(i,2),:));
    elseif map(i,1)==0
        thetainit_GT(i) = allTheta_GT(map(i,2),map(i,3));
    end
end
%% 
[newDudt,dudxyzibgddm,model,PSF,PSFderivatives] =  realDM_Derivatives_for_crlb_optim(thetainit_GT,parameters,map,shared,0);

NV = parameters.numparams;
t2 = 1./model; % [1].Smith, C.S., et al., Fast, single-molecule localization that achieves theoretically minimum uncertainty. Nature Methods, 2010. 7(5): p. 373-375.
hessian = zeros(NV,NV);
for m = 1:NV
    if map(m,1)~=2
        mm = map(m,2); 
        mmShared = map(m,1); 
        mmChannel = map (m,3); 
%         if map(m,2)<=21
%             Mnormterm = normLamda*newTheta(map(m,2));
%         else
%             Mnormterm = 0;
%         end
        temp1 = squeeze(newDudt(:,:,:,mm)); 
        for n = m:NV
%             if map(n,2)<=21
%                 Nnormterm = normLamda*newTheta(map(n,2));
%             else
%                 Nnormterm = 0;
%             end
            if map(n,1)~=2
                nn = map(n,2);
                nnShared = map(n,1);
                nnChannel = map(n,3);
                temp2 = squeeze(newDudt(:,:,:,nn));
                
                if mmShared ==1 &&nnShared ==1 
                    temp = t2.*temp1.*temp2;
                    hessian((m-1)*NV+n) = sum(temp(:));
                    hessian((n-1)*NV+m) = hessian((m-1)*NV+n);
                    
                elseif mmShared ==1&&nnShared==0 
                    temp = t2(:,:,nnChannel).*temp1(:,:,nnChannel).*temp2(:,:,nnChannel);
                    hessian((m-1)*NV+n) = sum(temp(:));
                    hessian((n-1)*NV+m) = hessian((m-1)*NV+n);
                 
                elseif nnShared ==1&&mmShared==0
                    temp = t2(:,:,mmChannel).*temp1(:,:,mmChannel).*temp2(:,:,mmChannel);
                    hessian((m-1)*NV+n) = sum(temp(:));
                    hessian((n-1)*NV+m) = hessian((m-1)*NV+n);
                    
                elseif mmShared ==0&&nnShared==0
                    if mmChannel == nnChannel
                        temp = t2(:,:,mmChannel).*temp1(:,:,mmChannel).*temp2(:,:,mmChannel);
                        hessian((m-1)*NV+n) = sum(temp(:));
                        hessian((n-1)*NV+m) = hessian((m-1)*NV+n);
                    end
                    
                end
            end
            
        end
    end
end
% CRLB = sqrt(diag(inv(hessian)));
% rcond(hessian)

%% loss_v1 
% % tmp_1 = sqrt(1./diag(hessian)); % crlb开根号
% tmp_1 = 1./diag(hessian); % crlb不开根号
% ind_start = numAberrations+1;
% ind_end = Nz*3+numAberrations;
% sum_crlb=sum(tmp_1(ind_start:ind_end));
% 
% % gradient for dloss/dzernike
% dlossdzer = zeros(numAberrations,1);
% 
% for jz=1:Nz
%     for jpos=1:3
%         
%         dIppdzer = zeros(numAberrations,1);
%         for jzer = 1:numAberrations
%             tmp = -1./(model(:,:,jz).^2) .*...
%                 (newDudt(:,:,jz,numAberrations+jpos).^2).*...
%                 newDudt(:,:,jz,jzer) + 2./model(:,:,jz).*...
%                 newDudt(:,:,jz,numAberrations+jpos).*...
%                 dudxyzibgdzer(:,:,jz,jpos,jzer);
%             dIppdzer(jzer) = sum(tmp(:));
%         end
%         
% %         % 
% %         dlossdzer = dlossdzer-1/2/...
% %             (hessian(numAberrations+(jpos-1)*Nz+jz,...
% %             numAberrations+(jpos-1)*Nz+jz)^(3/2))*dIppdzer;
%         % 
%         dlossdzer = dlossdzer-1/...
%             (hessian(numAberrations+(jpos-1)*Nz+jz,...
%             numAberrations+(jpos-1)*Nz+jz)^2)*dIppdzer;
%     end
% end

%% loss_v2 
loss = 0;
mask = [1,0,0,0,0;
        0,1,0,0,0;
        0,0,1,0,0;
        0,0,0,0,0
        0,0,0,0,0];

% x_crlb=zeros([1,Nz]);
% y_crlb=zeros([1,Nz]);
% z_crlb=zeros([1,Nz]);
subfishers = zeros([5,5,Nz]);
for jz=1:Nz
    for i=1:5
        for k =1:5
            subfishers(i,k,jz) = hessian(num_dmpixels+jz+(i-1)*Nz,...
                                     num_dmpixels+jz+(k-1)*Nz);
        end
    end
    weighted_fisher = inv( subfishers(:,:,jz) ).*mask;
    loss = loss+1/Nz*trace(weighted_fisher);
    
%     crlb_tmp =  diag( weighted_fisher ) ;
%     x_crlb(jz) = crlb_tmp(1);
%     y_crlb(jz) = crlb_tmp(2);
%     z_crlb(jz) = crlb_tmp(3);
end
% sum_crlb = sum(x_crlb(:))+sum(y_crlb(:))+sum(z_crlb(:));

% dfisher/dzernike
dfisherddm = zeros([5,5,Nz,num_dmpixels]);
for jzer=1:num_dmpixels
   for jz = 1:Nz
      for i=1:5
          for k=1:5
              tmp1 = newDudt(:,:,jz,num_dmpixels+i).*...
                     newDudt(:,:,jz,num_dmpixels+k).*( -1./ (model(:,:,jz).^2) ).*...
                     newDudt(:,:,jz,jzer);
              tmp2 = 1./model(:,:,jz).*newDudt(:,:,jz,num_dmpixels+i).*...
                     dudxyzibgddm(:,:,jz,k,jzer);
              tmp3 = 1./model(:,:,jz).*newDudt(:,:,jz,num_dmpixels+k).*...
                     dudxyzibgddm(:,:,jz,i,jzer);
              tmp = tmp1+tmp2+tmp3;
              dfisherddm(i,k,jz,jzer)=sum(tmp(:));
          end
      end
   end
end

% dloss/dzer
dlossdzer = zeros(num_dmpixels,1);
for jzer=1:num_dmpixels
    for jz=1:Nz
        tmp_inv = inv(subfishers(:,:,jz));
        dlossdzer(jzer)=dlossdzer(jzer)+...
            1/Nz*trace( mask.*( -tmp_inv*dfisherddm(:,:,jz,jzer)*tmp_inv ) );
    end
end

%% regularization term
if parameters.alpha~=0
    alpha = parameters.alpha;
    x = linspace(-(sizeX-1)/2,(sizeX-1)/2,sizeX);
    [x1,y1] = meshgrid(x,x);
    mask_reg = x1.^2+y1.^2;
    mask_reg = mask_reg*alpha;
    mask_reg = repmat(mask_reg,1,1,Nz);
    loss = loss + sum(mask_reg.*PSF,'all')/Nz;

    for jzer=1:num_dmpixels
        dlossdzer(jzer)=dlossdzer(jzer)+sum(PSFderivatives(:,:,:,jzer).*mask_reg,'all')/Nz;
    end
end
%% plot CRLB and model
if parameters.show
    x_crlb=zeros([1,Nz]);
    y_crlb=zeros([1,Nz]);
    z_crlb=zeros([1,Nz]);
    for j=1:Nz
        fisher_tmp = zeros(5);
        for i=1:5
            for k =1:5
                fisher_tmp(i,k) = hessian(num_dmpixels+j+(i-1)*Nz,...
                                         num_dmpixels+j+(k-1)*Nz);
            end
        end
        sqrt_crlb_tmp = sqrt(diag(inv(fisher_tmp)));
%         sqrt_crlb_tmp = diag(inv(fisher_tmp));
        x_crlb(j) = sqrt_crlb_tmp(1);
        y_crlb(j) = sqrt_crlb_tmp(2);
        z_crlb(j) = sqrt_crlb_tmp(3);
    end

    figure;plot(parameters.zemit,x_crlb,'b','displayname','CRLB^{1/2}_{y}');hold on;
    plot(parameters.zemit,y_crlb,'r','displayname','CRLB^{1/2}_{x}');hold on;
    plot(parameters.zemit,z_crlb,'g','displayname','CRLB^{1/2}_{z}');hold on;
    title(['CRLB^{1/2} at different z positions, ','photons:',...
        num2str(parameters.Nphotons(1)),',bg:',num2str(parameters.bg(1))])
    xlabel('z(nm)')
    ylabel('CRLB^{1/2}(nm)')
    legend

    disp(['average crlb_3D is: ',num2str( (sum(x_crlb.^2)+sum(y_crlb.^2)+sum(z_crlb.^2))/Nz )]);
    
    imageslicer(permute(model,[2 1 3]));% imageslicer 有一个转置的操作
end
