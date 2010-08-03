function [imle,imre,obj]=ssademo
imle=zeros(256);
obj=zeros(256);
obj(128,128)=100;
obj(100,100)=120;
obj(100,131)=200;
imre=zeros(511);
for n=1:1000
    imse=abs(ifft2(fft2(obj).*fft2(ifftshift(psfsse(1e-1,5e-2,5,256)))));
    imse=imse/sum(sum(imse));
    imle=imle+imse;
    h=zeros(256);
    [~,d1]=max(max(imse,[],2));
    [~,d2]=max(max(imse,[],1));
    h(257-d1,257-d2)=1;
    imre=conv2(imse,h)+(imre);
end
h=zeros(256);
[~,d1]=max(max(obj,[],2));
[~,d2]=max(max(obj,[],1));
h(257-d1,257-d2)=1;
obj=conv2(obj,h);
figure('Name','object')
imagesc(obj);
axis equal;
axis off;
h=zeros(256);
[~,d1]=max(max(imle,[],2));
[~,d2]=max(max(imle,[],1));
h(257-d1,257-d2)=1;
imle=conv2(imle,h);
figure('Name','long-exposure')
imagesc(imle);
axis equal;
axis off;
figure('Name','ssa reconstructed')
imagesc(imre);
axis equal;
axis off;
return