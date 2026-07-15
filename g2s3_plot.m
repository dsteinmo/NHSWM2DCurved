clear;
% Driver script for solving the 2D vacuum Euler's equations 
Globals2D;

% Order of polynomial approximation (N) 
N = 6;



%fluxtype = 'Roe'; %Roe not modified for SW yet.
%should also code an HLLC solver at some point.
fluxtype = 'HLL';
%fluxtype = 'LF';

%InitialSolution = @IsentropicVortexIC2D;
BCSolution = @ForwardStepBC2D;
ExactSolution = [];

% Read in Mesh
load('donut_peninsula.mat');
% Initialize solver and construct grid and metric
StartUp2D;


%return;



%Derek: initialize analytic boundary data and make a spline re-construction
%Note: the spline data can be built in mesh generation
rmax= 8345.2;
rmin=1000;
dtheta = pi/30;  %specify spacing at the bdry.
theta  = (-pi:dtheta:(pi-dtheta))';

rmax = rmax - 4e3*sech((theta - pi/2)/(pi/10)).^2;

%outer nodes
node   = [rmax.*cos(theta) rmax.*sin(theta)];
node = [node; node(1,:)];



     
% %inner nodes  %work in progress
% %might be okay to "cheat" like this, as long as you 
% %close off each loop. be careful though damn it!
% nodesinner = [rmin*cos(theta) rmin*sin(theta)];
% nodesinner = [nodesinner; nodesinner(1,:)];
% 
% %node = [node; rmin*cos(theta) rmin*sin(theta)];
% node = [node; nodesinner];

%arclength parameterization of boundary.
[xs,ys, t] = ParametricSpline(node(:,1),node(:,2));

figure(1);
clf;
set(gcf,'Renderer','Painters');
set(gcf,'PaperUnits','inches','PaperSize',[8 8],'PaperPosition',[0 0 8 8]);
%set(gcf,'DefaultLegendFontSize',16);
tt = linspace(min(t),max(t),20*length(t));
plot(node(:,1),node(:,2),'*b',ppval(xs,tt),ppval(ys,tt),'-r',[0 0],[0 0], '-k', 'linewidth',2);
h = legend('Raw boundary data','Spline Interpolant','Triangular Mesh','Location','SouthEast');
set(h,'FontSize',20);
hold on;
PlotMesh2D;
axis(1.0e+03 *[-1.2428    1.2594    3.7979    5.75]);
axis off;
%print('-dpng','-r900','spinevslinear.png');

%hard code elements to curve (again, can be identified in mesh generation)
%perhaps by looking for re-entrant corners.
k = [307 156 311 306]';
%flag only those faces that lie on the boundary
f=zeros(size(k));
for jj=1:length(k)
    f(jj) = find(BCType(k(jj),:)==Wall);
end
figure(2); clf;
set(gcf,'Renderer','Painters');
set(gcf,'PaperUnits','inches','PaperSize',[8 4],'PaperPosition',[0 0 8 4]);
subplot(1,2,1);
plot(x(Fmask(:,1),k(1)),y(Fmask(:,1),k(1)),'-k');
hold on;
plot(x(Fmask(:,2),k(1)),y(Fmask(:,2),k(1)),'-k');
plot(x(Fmask(:,3),k(1)),y(Fmask(:,3),k(1)),'-k','linewidth',2);
plot(x(Fmask(:,1),k(2)),y(Fmask(:,1),k(2)),'-k');
plot(x(Fmask(:,2),k(2)),y(Fmask(:,2),k(2)),'-k','linewidth',2);
plot(x(Fmask(:,3),k(2)),y(Fmask(:,3),k(2)),'-k');
%plot(VX(EToV(k(1),[1 2 3 1])),VY(EToV(k(1),[1 2 3 1])),'-k'); hold on;
%plot(VX(EToV(k(2),[1 2 3 1])),VY(EToV(k(2),[1 2 3 1])),'-k');
%plot(VX(EToV([k(1) k(2)],[1 2 3 1])),VY(EToV([k(1) k(2)],[1 2 3 1])),'-k');
hold on;
plot(x(:,k(1)),y(:,k(1)),'.','Color',[.5 0 .5],'markersize',15);
plot(x(:,k(2)),y(:,k(2)),'.','Color',[.5 0 .5],'markersize',15);
%plot(x(:,[k(1) k(2)]),y(:,[k(1) k(2)]),'.');
axis off;
 
% adjust curved elements if Baccording to simulation type
% hard wired, cylinder (centered at (0,0) with radius .5)
%[k,f] = find(BCType==Cyl);
curved = [];
if(~isempty(k))
 cylfaces = [k,f];
 curved = sort(unique(k));
 
 %MakeCylinder2D(cylfaces, .5, 0, 0);
 MakeCurvedEdges_derek(cylfaces,xs,ys,t);
 
 % turn cylinders into walls
 ids = find(BCType==Cyl);  BCType(ids) = Wall;
end
straight = setdiff(1:K, curved);
BuildBCMaps2D

subplot(1,2,2);
hold off;
plot(x(Fmask(:,1),k(1)),y(Fmask(:,1),k(1)),'-k');
hold on;
plot(x(Fmask(:,2),k(1)),y(Fmask(:,2),k(1)),'-k');
plot(x(Fmask(:,3),k(1)),y(Fmask(:,3),k(1)),'-k','linewidth',2);
plot(x(:,k(1)),y(:,k(1)),'.','Color',[.5 0 .5],'markersize',15);

plot(x(Fmask(:,1),k(2)),y(Fmask(:,1),k(2)),'-k');
plot(x(Fmask(:,2),k(2)),y(Fmask(:,2),k(2)),'-k','linewidth',2);
plot(x(Fmask(:,3),k(2)),y(Fmask(:,3),k(2)),'-k');
plot(x(:,k(2)),y(:,k(2)),'.','Color',[.5 0 .5],'markersize',15);
axis off;

print('-dpng','-r900','deformandblend.png');
return;

%set physical parameters.
g=0.0025*9.81;

Lx = max(x(:));
Hbar = 12.8;
H = Hbar*ones(Np,K);
eta0 = (Hbar/4)*(1/Lx)*x;

fRot = 7.8828e-5; %set f-plane rotation

%Q(:,:,1) = H+eta0;
% compute initial condition (time=0)
%Q = feval(InitialSolution, x, y, 0);
Q = zeros(Np,K,4);
Q(:,:,1) = H+eta0;
Q(:,:,2) = 0*Q(:,:,2);
Q(:,:,3) = 0*Q(:,:,3);
Q(:,:,4) = 0*Q(:,:,4);

%close all;

% Solve Problem
FinalTime = 3600*24*3;
NumOuts=200;
[Q] = CurvedEuler2D(Q, FinalTime, NumOuts, BCSolution, fluxtype,g,H,fRot); 
