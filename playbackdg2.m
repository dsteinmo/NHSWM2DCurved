%plotstart;

stuff = ls('out*.mat');
filelist = strread(stuff,'%s');
filelist = sort(filelist); %put into proper time ordering

numfiles = length(filelist);

for ii=1:numfiles
   filenamestr = filelist{ii};
   load(filenamestr);
   figure(2); colormap(darkjet); PlotField2D(N,x,y,eta);
   colorbar; caxis([-3.2 3.2]); title(['t=' num2str(time)]); view([0 90]);
   axis tight;
   drawnow; 

   print('-dpng', ['frame' sprintf('%07d',ii) '.png']);
end
   
