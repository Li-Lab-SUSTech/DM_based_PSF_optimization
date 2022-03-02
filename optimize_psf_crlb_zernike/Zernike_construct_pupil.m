function pupil = Zernike_construct_pupil(aberrations_orders, paraSim)
%% mask should be the same as in the paraSim
PupilSize = 1.0;
DxyPupil = 2*PupilSize/paraSim.Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
%%
Waberration = zeros(size(XYPupil));
orders = aberrations_orders;

zernikecoefs = orders(:,3);
normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));  
zernikecoefs = normfac.*zernikecoefs; 
allzernikes = get_zernikefunctions(orders,XPupil,YPupil);


for j = 1:numel(zernikecoefs)
  Waberration = Waberration+zernikecoefs(j)*squeeze(allzernikes(j,:,:));  
end


Waberration = Waberration.*ApertureMask;
pupil = Waberration/paraSim.lambda*2*pi;  % rad