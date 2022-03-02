function dm_basis=dm_basis_resize(dm_impulse,Nsize,smooth)
dm_basis_filtered=zeros(size(dm_impulse));
dm_basis = zeros(140,Nsize,Nsize);

%% mask should be the same as in the paraSim
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Nsize;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);

%%
for j=1:140
    if smooth
%         freq=fftshift(fft2(squeeze(dm_basis(j,:,:))));
%         freq_filtered = freq.*filter;
%         basis_filtered=ifft2(ifftshift(freq_filtered));
%         dm_basis(j,:,:) = real(basis_filtered).*ApertureMask;
        basis_filtered = medfilt2(squeeze(dm_impulse(j,:,:)),[5,5]);
        dm_basis_filtered(j,:,:) = real(basis_filtered);
    else
        dm_basis_filtered(j,:,:) = dm_impulse(j,:,:);
    end
    
    dm_basis(j,:,:)=imresize( squeeze(dm_basis_filtered(j,:,:)),[Nsize,Nsize] ).*ApertureMask;

end

%% look the dependent basis()
% tmp=zeros(Nsize,Nsize,140);
% for i=1:140
% tmp(:,:,i)=dm_basis(i,:,:);
% end
% imageslicer(tmp);
% disp('original dm_basis')
% 
% tmp_basis=tmp(:,:,64);
% figure;mesh(tmp_basis)
% freq=fftshift(fft2(tmp_basis));
% figure;mesh(abs(freq))
% 
% X = -Nsize/2+1:1:Nsize/2;
% [filtx,filty] = meshgrid(X,X);
% filter = double((filtx.^2+filty.^2)<21^2);
% figure;mesh(filter);
% freq_filt = freq.*filter;
% basis_filtered=ifft2(ifftshift(freq_filt));
% figure;mesh(real(basis_filtered))
% % 
% % % figure;mesh(abs(basis_filtered))
% 
% 
% tmp_basis=tmp(:,:,6);
% figure;mesh(tmp_basis)
% tmp_filtered = medfilt2(tmp_basis,[3,3]);
% figure;mesh(tmp_filtered)
