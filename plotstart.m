%TestAll_belize;
Globals2D;

% Order of polynomial approximation (N) 
N = 8;



%fluxtype = 'Roe'; %Roe not modified for SW yet.
%should also code an HLLC solver at some point.
fluxtype = 'HLL';
%fluxtype = 'LF';

%InitialSolution = @IsentropicVortexIC2D;
BCSolution = @ForwardStepBC2D;
ExactSolution = [];

% Read in Mesh
load('donut_peninsula.mat');


%return;

% Initialize solver and construct grid and metric
StartUp2D;

%Derek: initialize analytic boundary data and make a spline re-construction
%Note: the spline data can be built in mesh generation
rmax= 8345.2;
rmin=1000;
dtheta = pi/80;  %specify spacing at the bdry.
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

%figure(1);
%clf;

%PlotMesh2D;
%hold on;
%tt = linspace(min(t),max(t),20*length(t));
%plot(node(:,1),node(:,2),'.g',ppval(xs,tt),ppval(ys,tt),'-r');
%drawnow;

%hard code elements to curve (again, can be identified in mesh generation)
%perhaps by looking for re-entrant corners.
k = [1017 1018 1172 1173 1025 815]';
f = [2 3 3 3 3 2]';

 
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


%set physical parameters.
g=0.0025*9.81;

Lx = max(x(:));
Hbar = 12.8;
H = Hbar*ones(Np,K);
eta0 = (Hbar/4)*(1/Lx)*x;
%eta0 = 1*exp(-((x-.35*Lx)/1e3).^2-((y-.35*Lx)/1e3).^2);
% 
% figure(5);
% pf2d(N,x,y,eta0); shading interp; colorbar; drawnow;
% return;

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

return;
[Q] = CurvedEuler2D(Q, FinalTime, NumOuts, BCSolution, fluxtype,g,H,fRot); 
