function [Z] = Zernikcoefs_ls(p,orders_input,x,y,mask)
%% for loop, inner product for every zernike polynomials
% orders = orders_input(:,1:2);
% allzernikes = get_zernikefunctions(orders,x,y);
% for i = 1:length(orders(:,2))
%     Z(i) = 2*(orders(i,1)+1)./pi/( 1+double(orders(i,2)==0)).* ...
%            sum(p.*conj(squeeze(allzernikes(i,:,:))),'all');
% 
% end



%% least square
p(isnan(p))=0;   
mask(isnan(mask)) = 0;

orders = orders_input(:,1:2);
allzernikes = get_zernikefunctions(orders,x,y);

[dimx,dimy] = size(p);

zernike_reshaped = [];
for i =1:size(orders,1)
    zernike_reshaped(:,i)=reshape(squeeze(allzernikes(i,:,:)).*mask, dimx*dimy, 1);
end



pupil_reshaped = reshape(p,dimx*dimy,1);

Z = (zernike_reshaped'*zernike_reshaped)^(-1)*(zernike_reshaped'*pupil_reshaped);

end