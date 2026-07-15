%H=10;
g=0.0025*9.81;

%plotstart;

stuff = ls('out*.mat');
filelist = strread(stuff,'%s');
filelist = sort(filelist); %put into proper time ordering

numfiles = length(filelist);

APEint = zeros(numfiles,1);
KEint = zeros(numfiles,1);
TEint = zeros(numfiles,1);
times = zeros(numfiles,1);

for ii=1:numfiles
   filenamestr = filelist{ii};
   load(filenamestr);

   APE = 0.5*g*eta.^2;
   KE = 0.5*(H+eta).*(u.^2+v.^2);

   APEint(ii) = dgint(APE,V,J);
   KEint(ii) = dgint(KE,V,J);
   TEint(ii) = APEint(ii)+KEint(ii);
   times(ii) = time;
end
figure(2); clf;
plot(times,KEint,times,APEint,times,TEint);
legend('KE','APE','KE+APE');
xlabel('time (s)'); ylabel('Area-Integrated Energy');
drawnow;

print('-dpng','energyplot.png');
