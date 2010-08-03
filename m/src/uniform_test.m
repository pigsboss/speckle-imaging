function [i_obs,i_cal]=uniform_test(i_obj,M,K)
u_obj=ifft2(ifftshift(sqrt(i_obj)));
i_obs=zeros(size(i_obj));
u_cal=zeros(size(i_obj));
for k=1:K
	[X,Y]=meshgrid(0:7);
	[XI,YI]=meshgrid((0:size(i_obj,1)-1)*7/size(i_obj,1));
	real_add=interp2(X,Y,M*(0.5-rand(8)),XI,YI);
	imag_add=interp2(X,Y,M*(0.5-rand(8)),XI,YI);
	i_obs_k=abs(fftshift(fft2(u_obj.*((real_add)+1i*(imag_add))))).^2;
	i_obs=i_obs+i_obs_k;
	u_obs_k=sqrt(ifft2(ifftshift(i_obs_k)));
	u_cal=u_cal+u_obs_k;
end
i_cal=abs(fftshift(fft2(u_cal))).^2;
figure
imagesc(i_obs);
figure
imagesc(i_cal);
endfunction