function dm_coefs_ls = realDMbasis_decompose_ls(pupil,dm_basis)
%% mask should be the same as in the paraSim
PupilSize = 1.0;
DxyPupil = 2*PupilSize/size(pupil,1);
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);

%% least square
pupil(isnan(pupil))=0;   
ApertureMask(isnan(ApertureMask)) = 0;


[dimx,dimy] = size(pupil);

dmbasis_reshaped = [];
for i =1:size(dm_basis,1)
    dmbasis_reshaped(:,i)=reshape(squeeze(dm_basis(i,:,:)).*ApertureMask, dimx*dimy, 1);
end

pupil_reshaped = reshape(pupil,dimx*dimy,1);

% dm_coefs_ls = (dmbasis_reshaped'*dmbasis_reshaped)^(-1)*(dmbasis_reshaped'*pupil_reshaped);


%% ·ÇÂúÖÈ£¬ÐèÒª svd
[u,s,v]=svd(dmbasis_reshaped,'econ');

dm_coefs_ls = v*s^(-1)*u'*pupil_reshaped;

end