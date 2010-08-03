function [imre,imle]=ssa(filename,start_pos,end_pos)
%SSA Simple shift-and-add.

imle=zeros(end_pos(1:2)-start_pos(1:2)+[1,1]);
imre=imle;
for k=start_pos(3):end_pos(3)
  im=double(fits_read_image_subset(filename,...
    [start_pos(1),start_pos(2),k],...
    [end_pos(1),end_pos(2),k]));
  imle=imle+im;
end
imle=imle/(end_pos(3)-start_pos(3)+1);
x_m=floor(0.5*(start_pos(2)+end_pos(2)));
y_m=floor(0.5*(start_pos(1)+end_pos(1)));
for k=start_pos(3):end_pos(3)
  im=double(fits_read_image_subset(filename,...
    [start_pos(1),start_pos(2),k],...
    [end_pos(1),end_pos(2),k]));
  [~,x]=max(max(im,[],1));
  [~,y]=max(max(im,[],2));
  imre=imre+shift(im,x_m-x,y_m-y);
end
imre=imre/(end_pos(3)-start_pos(3)+1);
end