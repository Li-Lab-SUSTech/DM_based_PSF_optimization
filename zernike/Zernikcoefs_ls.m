function [zernike_coefs_ls] = Zernikcoefs_ls(pupil,orders_input)
%% mask should be the same as in the paraSim
PupilSize = 1.0;
DxyPupil = 2*PupilSize/size(pupil,1);
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);


%% least square
pupil(isnan(pupil))=0;   
ApertureMask(isnan(ApertureMask)) = 0;

orders = orders_input(:,1:2);
allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));

[dimx,dimy] = size(pupil);

zernike_reshaped = [];
for i =1:size(orders,1)
    zernike_reshaped(:,i)=reshape(squeeze(allzernikes(i,:,:)).*ApertureMask, dimx*dimy, 1);
end



pupil_reshaped = reshape(pupil,dimx*dimy,1);

% zernike_coefs_ls = (zernike_reshaped'*zernike_reshaped)^(-1)*(zernike_reshaped'*pupil_reshaped);

%% 非满秩，需要 svd
[u,s,v]=svd(zernike_reshaped,'econ');

zernike_coefs_ls = v*s^(-1)*u'*pupil_reshaped./normfac;

end