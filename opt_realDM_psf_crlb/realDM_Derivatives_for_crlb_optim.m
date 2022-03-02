function [dudt,dudxyzibgddm,model,PSF,PSFderivatives] =  realDM_Derivatives_for_crlb_optim(theta,parameters,map,shared,sigma)

% map[shared parameterID channelID]

lambda = parameters.lambda;
zemit0 = parameters.zemit0;
objStage0= parameters.objStage0; 
zemitStack = parameters.zemitStack;
objStageStack = parameters.objStageStack;
sizeX = parameters.sizeX;
sizeY = parameters.sizeY; 
sizeZ = parameters.sizeZ;
XPupil=parameters.XPupil;
PolarizationVector=parameters.PolarizationVector;
Amplitude=parameters.Amplitude;
wavevector=parameters.wavevector;
wavevectorzmed=parameters.wavevectorzmed;

dmbasis=parameters.dm_basis;

Ax=parameters.Ax;
Bx=parameters.Bx;
Dx=parameters.Dx;
Ay=parameters.Ay;
By=parameters.By;
Dy=parameters.Dy;
normIntensity = parameters.normIntensity;

numparams = parameters.numparams;
num_dmpixels = parameters.num_dmpixels;
ztype = parameters.ztype;

parameterID = num_dmpixels+5;
channelID = parameters.sizeZ;

repTheta = zeros(parameterID,channelID);

for i = 1:numparams
    if map(i,1)==0
        repTheta(map(i,2),map(i,3)) = theta(i);
    elseif map(i,1)==1||map(i,1)==2
        repTheta(map(i,2),1:channelID) = theta(i)*ones(1,channelID);
    end
end


% % calculation aberration function(zernike)
% Waberration = zeros(size(XPupil));
% zernikecoefs = squeeze(repTheta(1:numAberrations,1));% seems to work only when aberration is shared??
% zernikecoefs = normfac.*zernikecoefs;
% for j = 1:numel(zernikecoefs)
%   Waberration = Waberration+zernikecoefs(j)*squeeze(allzernikes(j,:,:));  
% end
% calculation aberration function(DM_basis)
Waberration = zeros(size(XPupil));
dm_coefs = squeeze(repTheta(1:num_dmpixels,1));% seems to work only when aberration is shared??
for j = 1:numel(dm_coefs)
  Waberration = Waberration+dm_coefs(j)*squeeze(dmbasis(j,:,:));  
end
PhaseFactor = exp(1i*Waberration);
% PhaseFactor = exp(2*pi*1i*Waberration/lambda);

numders = 3+num_dmpixels;

xemit = repTheta(num_dmpixels+1,:);
yemit = repTheta(num_dmpixels+2,:);
zemit = repTheta(num_dmpixels+3,:);
if strcmp(ztype,'stage')
    objStage = repTheta(num_dmpixels+3,:)'+objStageStack(:)+objStage0(:);
    zemit = zemitStack(:)+zemit0(:);
elseif strcmp(ztype,'emitter')
    objStage = objStageStack(:)+objStage0(:);
    zemit = repTheta(num_dmpixels+3,:)'+zemitStack(:)+zemit0(:);
end

%% dfield/dtheta
FieldMatrix = cell(2,3,sizeZ);
FieldMatrixDerivatives = cell(2,3,sizeZ,numders);

dudxyzddm_field = cell(2,3,sizeZ,3,num_dmpixels);

for jz = 1:sizeZ
    % xyz induced phase contribution
    Wxyz= (-1*xemit(jz))*wavevector{1}+(-1*yemit(jz))*wavevector{2}+(zemit(jz))*wavevectorzmed;
    
    PositionPhaseMask = exp(1i*(Wxyz+(objStage(jz))*wavevector{3}));
    
    for itel = 1:2
        for jtel = 1:3
            % pupil functions and FT to matrix elements
            PupilMatrix = Amplitude.*PhaseFactor.*PolarizationVector{itel,jtel};
            PupilFunction = PositionPhaseMask.*PupilMatrix;
            IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
            FieldMatrix{itel,jtel,jz} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
            
            % pupil functions for xy-derivatives and FT to matrix elements
            for jder = 1:2
                if shared(jder+num_dmpixels)~=2
                    PupilFunction_xy = -1i*wavevector{jder}.*PositionPhaseMask.*PupilMatrix;
                    IntermediateImage = transpose(cztfunc(PupilFunction_xy,Ay,By,Dy));
                    FieldMatrixDerivatives{itel,jtel,jz,num_dmpixels+jder} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
                end
            end
            
            % pupil functions for z-derivative and FT to matrix elements
            if shared(3+num_dmpixels)~=2
                if strcmp(ztype,'stage') % find stage position
                    PupilFunction_z = 1i*wavevector{3}.*PositionPhaseMask.*PupilMatrix;
                elseif strcmp(ztype,'emitter') % find emitter position
                    PupilFunction_z = 1i*wavevectorzmed.*PositionPhaseMask.*PupilMatrix;
                end
                IntermediateImage = transpose(cztfunc(PupilFunction_z,Ay,By,Dy));
                FieldMatrixDerivatives{itel,jtel,jz,num_dmpixels+3} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
            end
            
            % pupil functions for Zernike mode-derivative and FT to matrix elements
            %        n = 1;
            for jzer = 1:num_dmpixels
                if shared(jzer)~=2
                    PupilFunction_zer = 1i*squeeze(dmbasis(jzer,:,:)).*PositionPhaseMask.*PupilMatrix;
                    IntermediateImage = transpose(cztfunc(PupilFunction_zer,Ay,By,Dy));
                    FieldMatrixDerivatives{itel,jtel,jz,jzer} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
                end
            end
            
            %% dfield/(dxyzddm)
            for jpos = 1:3
                if jpos~=3
                    for jzer = 1:num_dmpixels
                        tmp = wavevector{jpos}.*squeeze(dmbasis(jzer,:,:)).*PositionPhaseMask.*PupilMatrix;
                        IntermediateImage = transpose(cztfunc(tmp,Ay,By,Dy));
                        dudxyzddm_field{itel,jtel,jz,jpos,jzer}=transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
                    end
                else
                    for jzer = 1:num_dmpixels
                        if strcmp(ztype,'stage')
                            tmp = -1*wavevector{3}.*squeeze(dmbasis(jzer,:,:)).*PositionPhaseMask.*PupilMatrix;
                        elseif strcmp(ztype,'emitter') 
                            tmp = -1*wavevectorzmed.*squeeze(dmbasis(jzer,:,:)).*PositionPhaseMask.*PupilMatrix;
                        end
                        IntermediateImage = transpose(cztfunc(tmp,Ay,By,Dy));
                        dudxyzddm_field{itel,jtel,jz,jpos,jzer}=transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
                    end                    
                end
            end
            %%
        end
    end
end

%% dpsf/dtheta
PSF = zeros(sizeX,sizeY,sizeZ);
PSFderivatives = zeros(sizeX,sizeY,sizeZ,numders);

dPSFdxyzddm = zeros(sizeX,sizeY,sizeZ,3,num_dmpixels);
for jz = 1:sizeZ
    for jtel = 1:3
        for itel = 1:2
            PSF(:,:,jz) = PSF(:,:,jz) + (1/3)*abs(FieldMatrix{itel,jtel,jz}).^2;
            for jder = 1:numders
                if shared(jder)~=2
                    PSFderivatives(:,:,jz,jder) = PSFderivatives(:,:,jz,jder) +...
                        (2/3)*real(conj(FieldMatrix{itel,jtel,jz}).*...
                        FieldMatrixDerivatives{itel,jtel,jz,jder});
                    % 1/3*real(conj(FieldMatrixDerivatives{itel,jtel,jz,jder}).*
                    % FieldMatrix{itel,jtel,jz}+conj(FieldMatrix{itel,jtel,jz}).*
                    % FieldMatrixDerivatives{itel,jtel,jz,jder});
                end
            end
            %% du/(dxyzdzer)
            for jpos=1:3
                for jzer = 1:num_dmpixels
                    dPSFdxyzddm(:,:,jz,jpos,jzer) = dPSFdxyzddm(:,:,jz,jpos,jzer)+...
                        (2/3)*real( conj(FieldMatrixDerivatives{itel,jtel,jz,jzer}).*...
                        FieldMatrixDerivatives{itel,jtel,jz,num_dmpixels+jpos}+...
                        conj(FieldMatrix{itel,jtel,jz}).*...
                        dudxyzddm_field{itel,jtel,jz,jpos,jzer});
%                     dudxyzdzer(:,:,jz,jpos,jzer) = dudxyzdzer(:,:,jz,jpos,jzer)+...
%                         (1/3)*real( conj(dudxyzdzer_field{itel,jtel,jz,jpos,jzer}).*...
%                         FieldMatrix{itel,jtel,jz}+...
%                         conj(FieldMatrixDerivatives{itel,jtel,jz,numAberrations+jpos}).*...
%                         FieldMatrixDerivatives{itel,jtel,jz,jzer}+...
%                         conj(FieldMatrixDerivatives{itel,jtel,jz,jzer}).*...
%                         FieldMatrixDerivatives{itel,jtel,jz,numAberrations+jpos}+...
%                         conj(FieldMatrix{itel,jtel,jz}).*dudxyzdzer_field{itel,jtel,jz,jpos,jzer});                    
                end
            end
            %%
        end
    end
end


PSF = PSF/normIntensity;
PSFderivatives = PSFderivatives/normIntensity;
dPSFdxyzddm =dPSFdxyzddm / normIntensity;

if sigma>0
    h = sigma;  %% otf 
    PSF = convn(PSF,h,'same');
    for jder = 1:numders
        PSFderivatives(:,:,:,jder) = convn(PSFderivatives(:,:,:,jder),h,'same');
    end
end


Nph = repTheta(num_dmpixels+4,:);
bg = repTheta(num_dmpixels+5,:);

model = PSF;
for i = 1:sizeZ
    model(:,:,i) = Nph(i)*PSF(:,:,i)+bg(i);
end

dudt = zeros(sizeX,sizeY,sizeZ,num_dmpixels+5);
dudxyzibgddm = zeros(sizeX,sizeY,sizeZ,5,num_dmpixels);
for i = 1:sizeZ
%     for jp = 1:2
%         if shared(numAberrations+jp)~=2
%             dudt(:,:,i,numAberrations+jp) = Nph(i)*PSFderivatives(:,:,i,numAberrations+jp);
%         end
%     end
%     for jp = 1:numAberrations
%         if shared(jp)~=2
%             dudt(:,:,i,jp) = Nph(i)*PSFderivatives(:,:,i,jp);
%         end  
%     end
%     if shared(numAberrations+3)~=2
%         dudt(:,:,i,numAberrations+3) = Nph(i)*PSFderivatives(:,:,i,numAberrations+3);
%     end
    dudt(:,:,i,1:num_dmpixels+3)=Nph(i)*PSFderivatives(:,:,i,:);
    dudxyzibgddm(:,:,i,1:3,:) = Nph(i)*dPSFdxyzddm(:,:,i,:,:);
    
end

dudt(:,:,:,num_dmpixels+4) = PSF;
dudt(:,:,:,num_dmpixels+5) = ones(size(PSF));

dudxyzibgddm(:,:,:,4,:) = PSFderivatives(:,:,:,1:num_dmpixels);
dudxyzibgddm(:,:,:,5,:) = zeros(size(dudt(:,:,:,1:num_dmpixels)));




