function imre=itsa(filename,start_pos,end_pos,NUMIT,imre)
%ITSA Iterative shift-and-add.

if nargin==4
  if numel(NUMIT)>1
    imre=NUMIT;
    NUMIT=1;
  else
    imre=ssa(filename,start_pos,end_pos);
  end
end
x_m=floor(0.5*(start_pos(2)+end_pos(2)));
y_m=floor(0.5*(start_pos(1)+end_pos(1)));
for l=1:NUMIT
  H=conj(fftn(imre));
  imre_new=zeros(size(imre));
  for k=start_pos(3):end_pos(3)
    im=double(fits_read_image_subset(filename,...
      [start_pos(1),start_pos(2),k],...
      [end_pos(1),end_pos(2),k]));
    corr_im=fftshift(ifftn(fftn(im).*H));
    [~,x]=max(max(corr_im,[],1));
    [~,y]=max(max(corr_im,[],2));
    imre_new=imre_new+shift(im,x_m-x,y_m-y);
  end
  imre=imre_new/(end_pos(3)-start_pos(3)+1);
end
return