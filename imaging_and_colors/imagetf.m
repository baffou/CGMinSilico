function imagetf(itf)

imagebw(-sqrt(abs(fftshift(fft2(itf)))))
set(gca,'YDir','normal')
colormap(gca,'Gray')
