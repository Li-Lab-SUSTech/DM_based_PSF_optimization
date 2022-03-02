function pupil_phase = DM_construct_pupil(dm_voltage_start, paraSim)
%% mask should be the same as in the paraSim
PupilSize = 1.0;
DxyPupil = 2*PupilSize/paraSim.Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);
%%
pupil_phase = zeros(size(XPupil));
for j = 1:numel(dm_voltage_start)
  pupil_phase = pupil_phase+dm_voltage_start(j)*squeeze(paraSim.dm_basis(j,:,:));  
end
pupil_phase = pupil_phase.*ApertureMask;  % rad