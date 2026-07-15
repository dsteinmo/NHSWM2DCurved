%clear;
%plotstart;
%load out0009567.mat

figure(1); clf; 
set(gcf,'renderer','painters');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [10 8],'PaperPosition',[0 0 10 8]);
colormap(darkjet);

pf2d(N,x,y,eta); shading interp; colorbar;
xlabel('x (m)'); ylabel('y (m)');
drawnow;
print('-dpng','-r300','donutpenin_nospureddy.png');
