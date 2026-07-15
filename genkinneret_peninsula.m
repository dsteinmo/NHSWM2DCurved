function [Nv,VX,VY,K,EToV,BCType,node,edge] = genkinneret_peninsula()
    % Geometry
    rmax= 8345.2;
    rmin=1000;
    %hdata.hmax=1000;
    dtheta = pi/30;  %specify spacing at the bdry.
    theta  = (-pi:dtheta:(pi-dtheta))';

    rmax = rmax - 4e3*sech((theta - pi/2)/(pi/10)).^2;
    
    %outer nodes
    node   = [rmax.*cos(theta) rmax.*sin(theta)];
    aa=length(node);
    
    %inner nodes
    node = [node; rmin*cos(theta) rmin*sin(theta)];
    
    %outer edges
    edge = [(1:aa-1)' (2:aa)'];
    edge = [edge; aa 1];
    
    %inner edges
    edge = [edge; (aa+1:length(node)-1)' (aa+2:length(node))'];
    edge = [edge; length(node) aa+1];
    
    

    % Make mesh
    %[Vert,EToV] = mesh2d(node,edge,hdata);
     [Vert,EToV] = mesh2d(node,edge);
    
    axis on;
    grid on; 
    xlabel('x (m)'); ylabel('y (m)');
    axis tight;
    dimp = size(Vert); dimt = size(EToV);
    Nv = dimp(1); K = dimt(1);

    disp(['Mesh generated. ' num2str(Nv) ' vertices on ' num2str(K) ' elements.']);

    %stuff below for DG
    VX = Vert(:,1); VY = Vert(:,2);

    %keyboard;

    % Reorder elements to ensure counter clockwise orientation
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    EToV(i,:) = EToV(i,[1 3 2]);
    %done reordering

    % Build connectivity matrix
    EToE = tiConnect2D(EToV);

    %find all boundary nodes (nodesOuter). Ain't nothin to it, but to do it.
    tol = 1e-8;
    edgenum = findedge(Vert,node,edge,tol); %this func is part of mesh2d.
    kk=1;
    for jj = 1:length(edgenum)
        %if edgenum of point jj is nonzero, then it lies on the boundary,
        %so put it in list of boundary pts.
        if edgenum(jj) ~= 0  
            nodesOuter(kk) = jj;
            kk=kk+1;
        end        
    end

    hold on;
    plot(VX(nodesOuter),VY(nodesOuter),'.r');
    drawnow;
    hold off;

    %allocate BCType table.
    BCType = 0*EToE;

    %Insert the correct BC flags for boundaries
    Wall=3;
    %BCType = CorrectBCTable_derek(EToV,BCType,nodesOuter,Wall,K);
    BCType = CorrectBCTable_derek(EToV,VX,VY,node,edge,BCType,Wall);
    %Need to do this to make vertex arrays consistent with main scripts.
    VX = VX';
    VY = VY';
    
    save('donut_peninsula.mat','Nv','VX','VY','K','EToV','BCType','node','edge')
end
