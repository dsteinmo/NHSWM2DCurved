%data on belize:~/RAID/curvedpeninsula/

clear;
close all;

plotstart;

stuff = ls('out*.mat');
filelist = strread(stuff,'%s');
filelist = sort(filelist); %put into proper time ordering

numfiles = length(filelist);

maxcont = 3.2;
mincont = -3.2;

figure(1); clf;
set(gcf,'renderer','painters');
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
        'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
          'DefaultAxesFontWeight','bold');
      
colormap(darkjet);

numpans = 8;

%ascii codes for abcd...
panlabels = 97:1:97+(numpans-1);

panHnds = zeros(numpans,1);
panTimes = zeros(numpans,1);

for j=1:length(filelist)
    disp(filelist{j});
    load(filelist{j},'eta','time');
    time(end)/3600
    
    eta(eta>maxcont) = maxcont;
    eta(eta<mincont) = mincont;
    
    panHnds(j) = subplot(2,4,j);
    pf2d(N,x,y,eta); shading interp;
    %caxis([-2 2]);
    caxis([mincont maxcont]);
    axis off;
    %xlabel('x (km)');
    %ylabel('y (km)');
    text(0.85*max(x(:)),0.85*max(y(:)),sprintf('(%c)',panlabels(j)));
    panTimes(j) =time;
    hold on;
    PlotDomain2D_black;
    
end
E=colorbar('Location','Southoutside');


% Specify some parameters for the plot
x0   = .135;  % spacing between and around figures (inches)
y0   = .8;    % offset from the bottom (inches)
w    = 1.58; % size of each subfigure (w x w inches)
y0cb = 0.25; % offset of colourbar from botom (inches)

%should make this part not-hard-coded as well...

set(panHnds(1),'Units','inches','Position',[x0       y0+w+x0  w w]);
set(panHnds(2),'Units','inches','Position',[w+2*x0   y0+w+x0  w w]);
set(panHnds(3),'Units','inches','Position',[2*w+3*x0 y0+w+x0  w w]);
set(panHnds(4),'Units','inches','Position',[3*w+4*x0 y0+w+x0  w w]);

set(panHnds(5),'Units','inches','Position',[x0        y0       w w]);
set(panHnds(6),'Units','inches','Position',[w+2*x0    y0       w w]);
set(panHnds(7),'Units','inches','Position',[2*w+3*x0  y0       w w]);
set(panHnds(8),'Units','inches','Position',[3*w+4*x0  y0       w w]);

set(E,'Units','inches','Position',[x0+0.25*w y0cb 1.5*w+2*w+3*x0 y0cb]);


%total width, height
W=5*x0 + 4*w;      % 5 spacing and four panels wide
H=2*w + 2*x0 + y0; % four spacing, fpir panels and one y0 in height
fprintf('Figure is %.2fin wide, %.2fin tall\n',W,H);
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [W H],'PaperPosition',[0 0 W H]);
set(gcf, 'Renderer', 'Painters');

filename = 'curvedpenin_conts_thesis.png';
filenamepdf = [filename(1:end-4) '.pdf'];
filenameeps = [filename(1:end-4) '.eps'];


print('-dpng','-r900',filename);
system(['convert ' filename ' ' filenameeps]);
system(['epstopdf ' filenameeps]);
system(['rm ' filenameeps]);
system(['cp ' filenamepdf ' ~/work/phdthesis/figures/']);
