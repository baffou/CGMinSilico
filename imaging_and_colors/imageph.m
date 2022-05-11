function imageph(W)

imagesc(W*1e9)
set(gca,'DataAspectRatio',[1 1 1])
colormap(gca,phase1024)
cb = colorbar;
ylabel(cb,'[nm]');
title('OPD image [nm]')
set(gca,'YDir','normal')