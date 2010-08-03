function [corr_im,spobj,spref]=speckleinterferometry(obj,start_pos_o,end_pos_o,ref,start_pos_r,end_pos_r)
%SPECKLEINTERFEROMETRY Speckle interferometry.

spobj=zeros(end_pos_o(1:2)-start_pos_o(1:2)+[1,1]);
spref=zeros(end_pos_r(1:2)-start_pos_r(1:2)+[1,1]);
if numel(spobj)~=numel(spref)
  error('obj and ref must have the same sizes.')
end
for k=start_pos_o(3):end_pos_o(3)
  im=double(fits_read_image_subset(obj,...
    [start_pos_o(1),start_pos_o(2),k],...
    [end_pos_o(1),end_pos_o(2),k]));
  im=im/mean(im(:));
  sp=fftn(im);
  sp=sp/mean(abs(sp(:)));
  spobj=spobj+abs(sp).^2;
end
spobj=spobj/(end_pos_o(3)-start_pos_o(3)+1);
for k=start_pos_r(3):end_pos_r(3)
  im=double(fits_read_image_subset(ref,...
    [start_pos_r(1),start_pos_r(2),k],...
    [end_pos_r(1),end_pos_r(2),k]));
  im=im/mean(im(:));
  sp=fftn(im);
  sp=sp/mean(abs(sp(:)));
  spref=spref+abs(sp).^2;
end
spref=spref/(end_pos_r(3)-start_pos_r(3)+1);
sp_im=spobj./spref;
corr_im=fftshift(ifftn(sp_im));
return