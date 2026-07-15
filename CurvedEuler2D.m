function Q = CurvedEuler2D(Q, FinalTime, NumOuts, ExactSolutionBC, fluxtype,g,H,fRot)

% function Q = CurvedEuler2D(Q, FinalTime, ExactSolution, ExactSolutionBC, fluxtype)
% Purpose  : Integrate 2D Euler equations using a 3rd order SSP RK

Globals2D;

%save nodal versions of normal from being overwritten by Gauss nodes.
nxnod = nx; nynod = ny;

% build cubature information
CubatureOrder = floor(2*(N+1)*3/2); 
cub = CubatureVolumeMesh2D(CubatureOrder);

% build Gauss node data
NGauss = floor((N+1)*2); gauss = GaussFaceMesh2D(NGauss);

%compute Hx & Hy. Not sure if you want to do it with cubature or what
%cH = cub.V*H;
Hx=0*H;
Hy=0*H;

disp('constructing Laplacian...');
[negLap,MM] = CurvedPoissonIPDG2D();
NHCoeff = (H(1,1)^2)/6;  %assumed to be scalar here
NHOP = -NHCoeff*negLap -MM;
clear negLap;
[ll,uu,pp,qq]= lu(NHOP);
clear NHOP;
Abc = CurvedPoissonIPDGbc2D();
disp('done');
%return;

% compute initial timestep
gamma = 1.4;
dt = EulerDT2D(Q, g); tstep = 1; time = 0;
rhsQ = 0*Q; resQ = 0*Q;
dtinit = dt;

OutInterval = FinalTime/NumOuts;
OutTimeCounter =0;

Filt = Filter2D(N,0.9*N,4);

%output at initial time
eta = Q(:,:,1)-H;
u = Q(:,:,2)./Q(:,:,1);
v = Q(:,:,3)./Q(:,:,1);
save(sprintf('out%07d.mat',tstep),'eta','u','v','N','K','x','y','time','H');
disp('outputting at t=0');

% outer time step loop 
figure(1); colormap(darkjet);
while (time<FinalTime)
  
  if(time+dt>FinalTime)
    dt = FinalTime-time;
  end
  
  % 3rd order SSP Runge-Kutta
  rhsQ  = CurvedEulerRHS2D(Q, time, ExactSolutionBC, fluxtype,g,Hx,Hy,fRot);
  for n=1:3
       rhsQ(:,:,n) = Filt*rhsQ(:,:,n);
  end
  Q1 = Q + dt*rhsQ;
  
  rhsQ  = CurvedEulerRHS2D(Q1, time, ExactSolutionBC, fluxtype,g,Hx,Hy,fRot);
  for n=1:3
       rhsQ(:,:,n) = Filt*rhsQ(:,:,n);
  end
  Q2 = (3*Q + Q1 + dt*rhsQ)/4;
  
  rhsQ  = CurvedEulerRHS2D(Q2, time, ExactSolutionBC, fluxtype,g,Hx,Hy,fRot);
  for n=1:3
       rhsQ(:,:,n) = Filt*rhsQ(:,:,n);
  end
  Q = (Q + 2*Q2 + 2*dt*rhsQ)/3;
   


  %compute NH correction, uses most recent momentum flux (note that my spectral codes
  %use momentum flux from previous time-step)
  a1 = rhsQ(:,:,2); a2 = rhsQ(:,:,3);
%   %evaluate at cubature nodes
%   ca1 = cub.V*a1; ca2 = cub.V*a2; 
%   %first step is getting the divergence of vector (a1,a2)
%   ddr = (cub.Dr')*(cub.W.*(cub.rx.*ca1 + cub.ry.*ca2));
%   dds = (cub.Ds')*(cub.W.*(cub.sx.*ca1 + cub.sy.*ca2));
%   diva = ddr+dds; %then interpolate to gauss nodes and do surface terms
%   %shouldn't do this: 
%   diva = MM\diva(:);
%   diva = reshape(diva,Np,K);

  %Note: since operator is linear, can do this with purely nodal DG:
  diva = Div2D(a1,a2);  %compute weak divergence
  dax = zeros(3*Nfp,K); 
  day = zeros(3*Nfp,K);
  dax(:) = (a1(vmapM)-a1(vmapP));  %form field differences at interfaces
  day(:) = (a2(vmapM)-a2(vmapP));
  
  fluxa = (nxnod.*dax + nynod.*day)/2.0; %compute central flux
  preRHS = -(diva - LIFT*(Fscale.*fluxa));  %rhs = - div a

  
  
  %deal with BC's
  %nodal code does:
  %qW = zeros(Nfp*Nfaces, K);
  %qW(mapW) = f*(-hv(vmapW).*nx(mapW) + hu(vmapW).*ny(mapW))./NHCoeff(vmapW);
  %interpolate hu & hv to gauss nodes (note: uses _interior_ points too to
  %get a good approximation on boundary)
  
  %ghu = gauss.interp*Q(:,:,2);
  %ghv = gauss.interp*Q(:,:,3);
  
  %zbc = zeros(gauss.NGauss*Nfaces*K,1);
  
  %the inhomogenous BC's seem to create numerical boundary layers at the
  %walls in plots of z, while homogenous bc's do not.
  %zbc(gauss.mapW) = fRot*(-ghv(gauss.mapW).*gauss.nx(gauss.mapW) + ... 
  %                     ghu(gauss.mapW).*gauss.ny(gauss.mapW))./NHCoeff;

  %zbc = Abc*zbc;
  
  myRHS = MM*(preRHS(:)); %homogenous bc.
%  myRHS = MM*(preRHS(:)+zbc); %junk at boundary happens regardless of bc sign
                              %think homogenous bc is the way to go.
                              
  z = qq*(uu\(ll\(pp*myRHS)));
  z = reshape(z,Np,K);
  
%   if mod(tstep,50) == 0 || tstep == 1
%       figure(1); clf;
%       %pf2d(N,x,y,Q(:,:,1)-H); colorbar;
%       pf2d(N,x,y,z); shading interp; colorbar;
%       title(['t=' num2str(time)]); drawnow;
%       disp('plotted');
%   end
  
  %Compute DG gradient of z
  [zx,zy] = Grad2D(z);  %first, get weak gradient
  %then handle fluxes
  
  %form field differences, %do i need to modified + traces for Dirichlet
  %conditions on u.n, or can i get away with this?
  dz = zeros(3*Nfp,K);
  dz(:) = (z(vmapM)-z(vmapP));
  
  %compute central fluxes
  fluxzx = nxnod.*dz/2.0; 
  fluxzy = nynod.*dz/2.0;   %Thought: may need to make sure zero-Dirichlet is imposed by modifying `+' traces here to re-inforce Dirichlet
                         %...i think you tried that and it was an epic fail.
  
  NHcorrx = NHCoeff.*(zx - LIFT*(Fscale.*fluxzx)); %compute NH correction to u&v eqns
  NHcorry = NHCoeff.*(zy - LIFT*(Fscale.*fluxzy));
  
  %time-step NH corrections.
  Q(:,:,2) = Q(:,:,2) + dt*NHcorrx;
  Q(:,:,3) = Q(:,:,3) + dt*NHcorry;

   
  %if it's time to output, then output!
  if OutTimeCounter >= OutInterval
      disp(['outputting at t=' num2str(time)]);
      eta = Q(:,:,1)-H;
      u = Q(:,:,2)./Q(:,:,1);
      v = Q(:,:,3)./Q(:,:,1);
      
    
      figure(1);
      subplot(3,1,1);
      pf2d(N,x,y,rhsQ(:,:,1));
      colormap(darkjet); shading interp; colorbar;
      subplot(3,1,2);
      pf2d(N,x,y,rhsQ(:,:,2));
      colormap(darkjet); shading interp; colorbar;
      subplot(3,1,3);
      pf2d(N,x,y,rhsQ(:,:,3));
      colormap(darkjet); shading interp; colorbar;

      save(sprintf('out%07d.mat',tstep),'eta','u','v','z','N','K','x','y','time','H');
      return;
    
    OutTimeCounter =0;
  end

  % Increment time and compute new timestep
  time = time+dt;
  OutTimeCounter = OutTimeCounter+dt;
  dt = EulerDT2D(Q, g);

  if dt <= 0.1*dtinit
       disp('time-step is getting puny. probably approaching an instability. Terminating...');
       return;
  end
  
  tstep = tstep+1;

  if mod(tstep,1000) == 0
      disp(['tstep = ' num2str(tstep)]);
  end

  if any(isnan(Q(:,:,1)))
     disp('NaNs found, terminating...');
  end
 
end;

%output at final time
eta = Q(:,:,1)-H;
u = Q(:,:,2)./Q(:,:,1);
v = Q(:,:,3)./Q(:,:,1);
save(sprintf('out%07d.mat',tstep),'eta','u','v','z','N','K','x','y','time','H');
disp(['outputting at t=' num2str(time)]);
return;
