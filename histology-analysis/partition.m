function localpart = partition(IDX, Ibinary, I, bprops, npartr, npartq)

global path groupnm outnm

global imsave ext

global stain vartype tissue

% Extract boundary properties (from boundaries.m)
ibound = bprops.ibound;    obound = bprops.obound;    centpt = bprops.centpt;

% Store copy of original boundaries
iboundi = ibound;    oboundi = flipud(obound);

% ----- Pseudo-radial streamline calculation ------------------------------
% Solve Laplace equation [del^2(phi) = div(grad(phi)) = 0] on tissue region to identify smooth trajectory between curves
% - Invert tissue mask and extract dimensions
phi0 = imcomplement(Ibinary);    [Nx0,Ny0] = size(phi0);

% - Build masks of boundaries to impose as boundary conditions
phi0i = false(Nx0,Ny0);    phi0o = false(Nx0,Ny0);
for II = 1:length(iboundi(:,1));    phi0i(iboundi(II,1),iboundi(II,2)) = 1;    end
for II = 1:length(oboundi(:,1));    phi0o(oboundi(II,1),oboundi(II,2)) = 1;    end

% - Find indices of all rows/columns of tissue mask
[r,c] = find(phi0);    rmin = min(r)-1;    rmax = max(r)+1;    cmin = min(c)-1;    cmax = max(c)+1;

% - Crop image to reduce degrees of freedom
phiR = phi0(rmin:rmax,cmin:cmax);    [Nx,Ny] = size(phiR);

% - Build mask of outer contour to impose as boundary condition
phiO = false(Nx,Ny);
for II = 1:length(oboundi(:,1));    phiO(oboundi(II,1)-rmin+1,oboundi(II,2)-cmin+1) = 1;    end

% - Extract indices of pixels inside/outside of region
U =  find(phiR);    % inside
W = find(~phiR);    % outside

% - Construct 4-point stencil for Laplacian solution (north, south, east, west)
UN = U - 1;    UE = U + Nx;    US = U + 1;    UW = U - Nx;

% - Initialize weights for each term in stencil
V = ones(size(U));

% - Sparse linear system coefficients for pixels INSIDE of region
IN = [U   U   1.00*V
    U  UN  -0.25*V
    U  UE  -0.25*V
    U  US  -0.25*V
    U  UW  -0.25*V ];

% - Sparse linear system coefficients for pixels OUTSIDE of region
OUT = [W  W  1.00*ones(size(W))];

% - Sparse linear system coefficients for ALL pixels in image
ALL = [IN; OUT];

% - Construct sparse matrix of coefficients
A = sparse(ALL(:,1), ALL(:,2), ALL(:,3));

% - Impose boundary conditions on OUTER boundary only
B = zeros(Nx,Ny);    B = B(:);    B(phiO(:)) = 1;

% - Solve sparse linear system
phiR = A\B;    phiR = reshape(phiR,[Nx,Ny]);

% - Project solution onto original-sized image;  phi is the solution to the Laplace equation
phi = zeros(Nx0,Ny0);    phi(rmin:rmax,cmin:cmax) = phiR;

% - Enforce the correct boundary/outside conditions
phi(phi0 == 0) = NaN;    phi(phi0i == 1) = 0;    phi(phi0o == 1) = 1;


% Compute components of the field TANGENT to the trajectory between inner and outer boundaries
[phiX,phiY] = gradient(phi);    Tx = phiX./sqrt(phiX.^2 + phiY.^2);    Ty = phiY./sqrt(phiX.^2 + phiY.^2);


% ---- Refine local boundaries :: CIRCUMFERENTIAL VARIATION ---------------
if strcmpi(vartype,'circ');
    
    % Initialize storage arrays for segment endpoints and line-type check
    ipart0 = zeros(npartq,2);
    
    % Define number of points in each partition
    intsize = round(linspace(1,length(iboundi(:,1)),(npartq+1)));
    
    % Divide inner contour based on partition size
    for II = 1:npartq;       ipart0(II,:) = iboundi(intsize(II),:);      end
    
    % Compute streamlines of Laplace solution over region
    % - Define grid over image
    [X,Y] = meshgrid(linspace(1,Ny0,Ny0),linspace(1,Nx0,Nx0));
    
    % - Define 2pixel-wide neighborhood
    nbh = 2.*[0,1; 0,-1; 1,0; -1,0; 1,1; 1,-1; -1,1; -1,-1];
    
    % - Compute streamlines at each inner partition points
    sxp   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(1,2)),(ipart0(:,1)+nbh(1,1)),[0.1,1e5]);
    sxn   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(2,2)),(ipart0(:,1)+nbh(2,1)),[0.1,1e5]);
    syp   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(3,2)),(ipart0(:,1)+nbh(3,1)),[0.1,1e5]);
    syn   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(4,2)),(ipart0(:,1)+nbh(4,1)),[0.1,1e5]);
    sxpyp = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(5,2)),(ipart0(:,1)+nbh(5,1)),[0.1,1e5]);
    sxpyn = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(6,2)),(ipart0(:,1)+nbh(6,1)),[0.1,1e5]);
    sxnyp = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(7,2)),(ipart0(:,1)+nbh(7,1)),[0.1,1e5]);
    sxnyn = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(8,2)),(ipart0(:,1)+nbh(8,1)),[0.1,1e5]);
    
    % - Concatenate all streamlines
    SL = cat(1,sxp,sxn,syp,syn,sxpyp,sxpyn,sxnyp,sxnyn);    [SLx,SLy] = size(SL);
    
    % Initialize storage for inner and outer boundary points
    ipart = NaN(SLy,2);    iidx = NaN(SLy,1);
    opart = NaN(SLy,2);    oidx = NaN(SLy,1);
    
    % Find closest points to successful streamlines on the inner and outer contours
    for II = 1:SLy;
        
        % Initialize logical array to keep streamlines
        keep = true(SLx,1);
        
        % Streamlines with NaNs did not traverse across region
        for JJ = 1:SLx;
            if any(isnan(SL{JJ,II}(2,:)))
                keep(JJ,1) = false;
            end
        end
        
        % Initialize index arrays for inner and outer boundaries
        kidx = find(keep)';    idxI = NaN(length(kidx),1);    idxO = NaN(length(kidx),1);
        
        % Compute pointwise distance between streamlines and inner/outer boundary
        for JJ = 1:length(kidx);
            fidx = find(isfinite(SL{kidx(JJ),II}(:,1)));
            
            pchkI = fliplr(SL{kidx(JJ),II}(1,:));    pchkO = fliplr(SL{kidx(JJ),II}(fidx(end),:));
            
            dI = pdist2(pchkI,iboundi);    dImin = min(dI);
            dO = pdist2(pchkO,oboundi);    dOmin = min(dO);
            
            chkI = find(dI == dImin)';    idxI(JJ,1) = chkI(1);
            chkO = find(dO == dOmin)';    idxO(JJ,1) = chkO(1);
        end
        
        % Assign inner/outer boundary point based on minimum distance to streamlines
        [~,sidxI] = min(abs(round(idxI - mean(idxI))));
        [~,sidxO] = min(abs(round(idxO - mean(idxO))));
        
        % Store partitions on inner/outer boudary
        iidx(II,1)  = idxI(sidxI);                oidx(II,1) = idxO(sidxO);
        ipart(II,:) = iboundi(idxI(sidxI),:);    opart(II,:) = oboundi(idxO(sidxO),:);
    end
    
    % Update index arrays to close partitions
    iidx = cat(1,iidx,iidx(1));    oidx = cat(1,oidx,oidx(1));
    
    %     % --------------------
    %     % Plot streamline paths on Laplace solution
    %     figure; imagesc(phi)
    %     hold on;
    %     for II  = 1:length(sxp);
    %         plot(sxp{II}(:,1),sxp{II}(:,2),'r');      plot(sxn{II}(:,1),sxn{II}(:,2),'g');      plot(syp{II}(:,1),syp{II}(:,2),'b');      plot(syn{II}(:,1),syn{II}(:,2),'w')
    %         plot(sxpyp{II}(:,1),sxpyp{II}(:,2),'c');  plot(sxpyn{II}(:,1),sxpyn{II}(:,2),'m');  plot(sxnyp{II}(:,1),sxnyp{II}(:,2),'y');  plot(sxnyn{II}(:,1),sxnyn{II}(:,2),'k')
    %     end
    %
    %     plot(ipart(:,2),ipart(:,1),'*');    plot(opart(:,2),opart(:,1),'*')
    %     for II = 1:npartq;    plot([ipart(II,2) opart(II,2)],[ipart(II,1) opart(II,1)],'w');    end
    %     % --------------------
    
    % Re-define boundaries as cell arrays for variational analysis
    ibound = {ibound};    obound = {obound};
    
    h = waitbar(0,'Building CIRCUMFERENTIAL Partitions...','Name',outnm{IDX,:});
    
    % Extract boundaries for each partition
    lbound = cell(npartq,1);    rbound = cell(npartq,1);
    for II = 1:npartq;
        
        % Update progress bar
        waitbar(II/npartq)
        
        % Isolate inner boundary from partitioned section
        ichk = [iidx(II,1) iidx(II+1,1)];
        if ichk(2) < ichk(1)
            ibound{II,1} = iboundi([ichk(1):end,1:ichk(2)],:);
        else
            ibound{II,1} = iboundi(ichk(1):ichk(2),:);
        end
        
        % Isolate inner boundary from partitioned section
        ochk = [oidx(II,1) oidx(II+1,1)];
        if ochk(2) < ochk(1)
            obound{II,1} = oboundi([ochk(1):end,1:ochk(2)],:);
        else
            obound{II,1} = oboundi(ochk(1):ochk(2),:);
        end
        
        % Construct mask of local partition
        % - Polygon defining local area
        ypoly = cat(1,ibound{II,1}(1,1),ibound{II,1}(end,1),obound{II,1}(end,1),obound{II,1}(1,1));
        xpoly = cat(1,ibound{II,1}(1,2),ibound{II,1}(end,2),obound{II,1}(end,2),obound{II,1}(1,2));
        
        % - Create mask from polygon
        Imsk = inpoly(cat(2,X(:),Y(:)),cat(2,xpoly,ypoly));    Imsk = reshape(Imsk,size(X));

        % - Number of objects in mask image
        cc = bwconncomp(Imsk,8);    num = cc.NumObjects;
        
        % - Make sure all corners are in object boundary by updating polygon and mask, if needed
        if num > 1
            ypoly = cat(1,ibound{II,1}(:,1),flipud(obound{II,1}(:,1)));
            xpoly = cat(1,ibound{II,1}(:,2),flipud(obound{II,1}(:,2)));
            Imsk = inpoly(cat(2,X(:),Y(:)),cat(2,xpoly,ypoly));    Imsk = reshape(Imsk,size(X));
        end
        
        % Extract object boundary
        bound = bwtraceboundary(Imsk,[ypoly(1), xpoly(1)], 'N', 8);
        
        [~,iidx1] = ismember(ibound{II,1}(1,:),bound,'rows');    [~,iidx2] = ismember(ibound{II,1}(end,:),bound,'rows');
        [~,oidx1] = ismember(obound{II,1}(1,:),bound,'rows');    [~,oidx2] = ismember(obound{II,1}(end,:),bound,'rows');
        
        % Extract and define left and right bounds
        if II == 1;
            lbound{II,1} = bound(iidx1:oidx1,:);    rbound{II,1} = bound(oidx2:iidx2,:);
        elseif II < npartq;
            lbound{II,1} = rbound{II-1,1};          rbound{II,1} = bound(oidx2:iidx2,:);
        elseif II == npartq;
            lbound{II,1} = rbound{II-1,1};          rbound{II,1} = lbound{1,1};
        end
        
    end
    
    close(h)
    
    % ---- Refine local boundaries :: RADIAL VARIATION ------------------------
elseif strcmpi(vartype,'rad')
    
    % Re-define boundaries as cell arrays for variational analysis
    ibound = {ibound};    obound = {obound};
    
    % Management of left and right boundaries
    lbound = cell(1,npartr);    rbound = cell(1,npartr);
    
    % Determine pointwise radial step size and disk size
    rsz = 1/npartr;
    
    h = waitbar(0,'Building RADIAL Partitions...','Name',outnm{IDX,:});
    
    % Reconstruct image from Laplace solution
    for II = 1:npartr;
        
        % Update progress bar
        waitbar(II/npartr)
        
        % Initialize current radius level
        rlevel = II*rsz;
        
        % Threshold Laplace solution and compute skeleton to identify level
        Ir = abs(phi - rlevel) < 1e-2;    Ir = bwmorph(Ir,'thin',Inf);
        
        % Get pixel coordinates of radius level
        rb = find(Ir);    [rbx,rby] = ind2sub([Nx0,Ny0],rb);     [rbx,rby] = points2contour(rbx,rby,1,'cw');    rmid = cat(2,rbx',rby');
        
        % Re-order contour based on centroid location
        for JJ = centpt(1)-0.5:Ny0;
            if Ir(centpt(2)-0.5,JJ) > 0 && Ir(centpt(2)-0.5,JJ-1) == 0;    mpt = [JJ centpt(2)-0.5];    end
        end
        
        % Find location of start point in boundary array and shift accordingly
        [~,idx] = ismember(fliplr(mpt),rmid,'rows');    rmid = circshift(rmid,[-idx+1,0]);
        
        % Update inner/outer boundaries
        obound{1,II} = rmid;      if II < npartr;  ibound{1,II+1} = obound{1,II};  end
        
    end
    
    close(h)
    
    % ---- Refine local boundaries :: CIRCUMFERENTIAL & RADIAL VARIATION ------
elseif strcmpi(vartype,'both');
    
    % First, perform RADIAL partitioning ----------------------------------
    iboundr = {iboundi};    oboundr = {oboundi};
    
    % Determine pointwise radial step size and disk size
    rsz = 1/npartr;
    
    h = waitbar(0,'Building RADIAL Partitions...','Name',outnm{IDX,:});
    
    % Reconstruct image from Laplace solution
    for II = 1:npartr;
        
        % Update progress bar
        waitbar(II/npartr)
        
        % Initialize current radius level
        rlevel = II*rsz;
        
        % Threshold Laplace solution and compute skeleton to identify level
        Ir = abs(phi - rlevel) < 1e-2;    Ir = bwmorph(Ir,'thin',Inf);
        
        % Get pixel coordinates of radius level
        rb = find(Ir);    [rbx,rby] = ind2sub([Nx0,Ny0],rb);     [rbx,rby] = points2contour(rbx,rby,1,'cw');    rmid = cat(2,rbx',rby');
        
        % Re-order contour based on centroid location
        for JJ = centpt(1)-0.5:Ny0;
            if Ir(centpt(2)-0.5,JJ) > 0 && Ir(centpt(2)-0.5,JJ-1) == 0;    mpt = [JJ centpt(2)-0.5];    end
        end
        
        % Find location of start point in boundary array and shift accordingly
        [~,idx] = ismember(fliplr(mpt),rmid,'rows');    rmid = circshift(rmid,[-idx+1,0]);
        
        % Update inner/outer boundaries
        oboundr{1,II} = rmid;      if II < npartr;  iboundr{1,II+1} = oboundr{1,II};  end
        
    end
    
    close(h)
    
    
    % Next, perform CIRCUMFERENTIAL partitioning --------------------------
    % Initialize storage arrays for segment endpoints and line-type check
    ipart0 = zeros(npartq,2);
    
    % Define number of points in each partition
    intsize = round(linspace(1,length(iboundi(:,1)),(npartq+1)));
    
    % Divide inner contour based on partition size
    for II = 1:npartq;       ipart0(II,:) = iboundi(intsize(II),:);      end
    
    % Compute streamlines of Laplace solution over region
    % - Define grid over image
    [X,Y] = meshgrid(linspace(1,Ny0,Ny0),linspace(1,Nx0,Nx0));
    
    % - Define 2pixel-wide neighborhood
    nbh = 2.*[0,1; 0,-1; 1,0; -1,0; 1,1; 1,-1; -1,1; -1,-1];
    
    % - Compute streamlines at each inner partition points
    sxp   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(1,2)),(ipart0(:,1)+nbh(1,1)),[0.1,1e5]);
    sxn   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(2,2)),(ipart0(:,1)+nbh(2,1)),[0.1,1e5]);
    syp   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(3,2)),(ipart0(:,1)+nbh(3,1)),[0.1,1e5]);
    syn   = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(4,2)),(ipart0(:,1)+nbh(4,1)),[0.1,1e5]);
    sxpyp = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(5,2)),(ipart0(:,1)+nbh(5,1)),[0.1,1e5]);
    sxpyn = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(6,2)),(ipart0(:,1)+nbh(6,1)),[0.1,1e5]);
    sxnyp = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(7,2)),(ipart0(:,1)+nbh(7,1)),[0.1,1e5]);
    sxnyn = stream2(X,Y,Tx,Ty,(ipart0(:,2)+nbh(8,2)),(ipart0(:,1)+nbh(8,1)),[0.1,1e5]);
    
    % - Concatenate all streamlines
    SL = cat(1,sxp,sxn,syp,syn,sxpyp,sxpyn,sxnyp,sxnyn);    [SLx,SLy] = size(SL);
    
    % Initialize storage for inner and outer boundary points
    ipart = NaN(SLy,2);    iidx = NaN(SLy,1);
    opart = NaN(SLy,2);    oidx = NaN(SLy,1);
    
    % Find closest points to successful streamlines on the inner and outer contours
    for II = 1:SLy;
        
        % Initialize logical array to keep streamlines
        keep = true(SLx,1);
        
        % Streamlines with NaNs did not traverse across region
        for JJ = 1:SLx;
            if any(isnan(SL{JJ,II}(2,:)))
                keep(JJ,1) = false;
            end
        end
        
        % Initialize index arrays for inner and outer boundaries
        kidx = find(keep)';    idxI = NaN(length(kidx),1);    idxO = NaN(length(kidx),1);
        
        % Compute pointwise distance between streamlines and inner/outer boundary
        for JJ = 1:length(kidx);
            fidx = find(isfinite(SL{kidx(JJ),II}(:,1)));
            
            pchkI = fliplr(SL{kidx(JJ),II}(1,:));    pchkO = fliplr(SL{kidx(JJ),II}(fidx(end),:));
            
            dI = pdist2(pchkI,iboundi);    dImin = min(dI);
            dO = pdist2(pchkO,oboundi);    dOmin = min(dO);
            
            chkI = find(dI == dImin)';    idxI(JJ,1) = chkI(1);
            chkO = find(dO == dOmin)';    idxO(JJ,1) = chkO(1);
        end
        
        % Assign inner/outer boundary point based on minimum distance to streamlines
        [~,sidxI] = min(abs(round(idxI - mean(idxI))));
        [~,sidxO] = min(abs(round(idxO - mean(idxO))));
        
        % Store partitions on inner/outer boudary
        iidx(II,1)  = idxI(sidxI);                oidx(II,1) = idxO(sidxO);
        ipart(II,:) = iboundi(idxI(sidxI),:);    opart(II,:) = oboundi(idxO(sidxO),:);
    end
    
    % Update index arrays to close partitions
    iidx = cat(1,iidx,iidx(1));    oidx = cat(1,oidx,oidx(1));
    
    %     % --------------------
    %     % Plot streamline paths on Laplace solution
    %     figure; imagesc(phi)
    %     hold on;
    %     for II  = 1:length(sxp);
    %         plot(sxp{II}(:,1),sxp{II}(:,2),'r');      plot(sxn{II}(:,1),sxn{II}(:,2),'g');      plot(syp{II}(:,1),syp{II}(:,2),'b');      plot(syn{II}(:,1),syn{II}(:,2),'w')
    %         plot(sxpyp{II}(:,1),sxpyp{II}(:,2),'c');  plot(sxpyn{II}(:,1),sxpyn{II}(:,2),'m');  plot(sxnyp{II}(:,1),sxnyp{II}(:,2),'y');  plot(sxnyn{II}(:,1),sxnyn{II}(:,2),'k')
    %     end
    %
    %     plot(ipart(:,2),ipart(:,1),'*');    plot(opart(:,2),opart(:,1),'*')
    %     for II = 1:npartq;    plot([ipart(II,2) opart(II,2)],[ipart(II,1) opart(II,1)],'w');    end
    %     % --------------------
    
    % Re-define boundaries as cell arrays for variational analysis
    iboundq = {iboundi};    oboundq = {oboundi};
    
    h = waitbar(0,'Building CIRCUMFERENTIAL Partitions...','Name',outnm{IDX,:});
    
    % Extract boundaries for each partition
    lboundq = cell(npartq,1);    rboundq = cell(npartq,1);
    for II = 1:npartq;
        
        % Update progress bar
        waitbar(II/npartq)
        
        % Isolate inner boundary from partitioned section
        ichk = [iidx(II,1) iidx(II+1,1)];
        if ichk(2) < ichk(1)
            iboundq{II,1} = iboundi([ichk(1):end,1:ichk(2)],:);
        else
            iboundq{II,1} = iboundi(ichk(1):ichk(2),:);
        end
        
        % Isolate inner boundary from partitioned section
        ochk = [oidx(II,1) oidx(II+1,1)];
        if ochk(2) < ochk(1)
            oboundq{II,1} = oboundi([ochk(1):end,1:ochk(2)],:);
        else
            oboundq{II,1} = oboundi(ochk(1):ochk(2),:);
        end
        
        % Construct mask of local partition
        % - Polygon defining local area
        ypoly = cat(1,iboundq{II,1}(1,1),iboundq{II,1}(end,1),oboundq{II,1}(end,1),oboundq{II,1}(1,1));
        xpoly = cat(1,iboundq{II,1}(1,2),iboundq{II,1}(end,2),oboundq{II,1}(end,2),oboundq{II,1}(1,2));
        
        % - Create mask from polygon
        Imsk = inpoly(cat(2,X(:),Y(:)),cat(2,xpoly,ypoly));    Imsk = reshape(Imsk,size(X));
        
        % - Number of objects in mask image
        cc = bwconncomp(Imsk,8);    num = cc.NumObjects;
        
        % - Make sure all corners are in object boundary by updating polygon and mask, if needed
        if num > 1
            ypoly = cat(1,iboundq{II,1}(:,1),flipud(oboundq{II,1}(:,1)));
            xpoly = cat(1,iboundq{II,1}(:,2),flipud(oboundq{II,1}(:,2)));
            Imsk = inpoly(cat(2,X(:),Y(:)),cat(2,xpoly,ypoly));    Imsk = reshape(Imsk,size(X));
        end
        
        % Extract object boundary
        bound = bwtraceboundary(Imsk,[ypoly(1), xpoly(1)], 'N', 8);
        
        [~,iidx1] = ismember(iboundq{II,1}(1,:),bound,'rows');    [~,iidx2] = ismember(iboundq{II,1}(end,:),bound,'rows');
        [~,oidx1] = ismember(oboundq{II,1}(1,:),bound,'rows');    [~,oidx2] = ismember(oboundq{II,1}(end,:),bound,'rows');
        
        % Extract and define left and right bounds
        if II == 1;
            lboundq{II,1} = bound(iidx1:oidx1,:);     rboundq{II,1} = bound(oidx2:iidx2,:);
        elseif II < npartq;
            lboundq{II,1} = rboundq{II-1,1};          rboundq{II,1} = bound(oidx2:iidx2,:);
        elseif II == npartq;
            lboundq{II,1} = rboundq{II-1,1};          rboundq{II,1} = lboundq{1,1};
        end
        
    end
    
    close(h)
    
    % Finally, combine RADIAL and CIRCUMFERENTIAL partitions --------------
    % Initialize storage for each partition
    ibound = cell(npartq,npartr);    obound = cell(npartq,npartr);    lbound = cell(npartq,npartr);    rbound = cell(npartq,npartr);
    lidxq  = cell(npartq,npartr);    ridxq  = cell(npartq,npartr);    lidxr  = cell(npartq,npartr);    ridxr  = cell(npartq,npartr);
    
    h = waitbar(0,'Combining Partitions...','Name',outnm{IDX,:});
    
    % Extract boundaries for each partition
    for KK = 1:npartr;
        for II = 1:npartq;
            
            % Update progress bar
            waitbar((npartq*(KK-1)+II)/(npartq*npartr));
            
            % Intersection of left / radial boundaries
            d = pdist2(lboundq{II,1},oboundr{1,KK},'euclidean');    idx = find(d == min(min(d)));
            [lidxq{II,KK},lidxr{II,KK}] = ind2sub(size(d),idx(1));
            
            % Intersection of right / radial boundaries
            d = pdist2(rboundq{II,1},oboundr{1,KK},'euclidean');    idx = find(d == min(min(d)));
            [ridxq{II,KK},ridxr{II,KK}] = ind2sub(size(d),idx(1));
            
            % Extract inner boundary from partitioned section
            if KK == 1;
                ibound{II,KK} = iboundq{II,KK};
            else
                ibound{II,KK} = obound{II,KK-1};
            end
            
            % Extract left/right boundaries from partitioned section
            if KK == 1
                if II == 1;
                    lbound{II,KK} = flipud(lboundq{II,1}(1:lidxq{II,KK},:));
                    rbound{II,KK} = rboundq{II,1}(ridxq{II,KK}:end,:);
                elseif II == npartq;
                    lbound{II,KK} = lboundq{II,1}(lidxq{II,KK}:end,:);
                    rbound{II,KK} = flipud(rboundq{II,1}(1:ridxq{II,KK},:));
                else
                    lbound{II,KK} = lboundq{II,1}(lidxq{II,KK}:end,:);
                    rbound{II,KK} = rboundq{II,1}(ridxq{II,KK}:end,:);
                end
            else
                if II == 1;
                    lbound{II,KK} = flipud(lboundq{II,1}(lidxq{II,KK-1}:lidxq{II,KK},:));
                    rbound{II,KK} = rboundq{II,1}(ridxq{II,KK}:ridxq{II,KK-1},:);
                elseif II == npartq;
                    lbound{II,KK} = lboundq{II,1}(lidxq{II,KK}:lidxq{II,KK-1},:);
                    rbound{II,KK} = flipud(rboundq{II,1}(ridxq{II,KK-1}:ridxq{II,KK},:));
                else
                    lbound{II,KK} = lboundq{II,1}(lidxq{II,KK}:lidxq{II,KK-1},:);
                    rbound{II,KK} = rboundq{II,1}(ridxq{II,KK}:ridxq{II,KK-1},:);
                end
            end
            
            % Extract outer boundary from partitioned section
            if lidxr{II,KK} > ridxr{II,KK}
                obound{II,KK} = oboundr{1,KK}(ridxr{II,KK}:lidxr{II,KK},:);
            else
                obound{II,KK} = oboundr{1,KK}([ridxr{II,KK}:end,1:lidxr{II,KK}],:);
            end
            
        end
    end
    
    close(h)
    
end

% Compile boundaries from each local partition
if strcmpi(vartype,'rad');     npartq = 1;    end
if strcmpi(vartype,'circ');    npartr = 1;    end

localpart = cell(npartq,npartr);

if strcmpi(vartype,'circ') || strcmpi(vartype,'rad');
    for KK = 1:npartr;
        for II = 1:npartq;
            localpart{II,KK} = cat(1, ibound{II,KK}, rbound{II,KK}, flipud(obound{II,KK}), lbound{II,KK});
        end
    end
elseif strcmpi(vartype,'both');
    for KK = 1:npartr;
        for II = 1:npartq;
            if KK == 1;
                localpart{II,KK} = cat(1, ibound{II,KK}, flipud(rbound{II,KK}), obound{II,KK}, lbound{II,KK});
            else
                localpart{II,KK} = cat(1, flipud(ibound{II,KK}), flipud(rbound{II,KK}), obound{II,KK}, lbound{II,KK});
            end
        end
    end
end

if strcmpi(tissue,'Myocardial Infarction         ');
    laplace_eqn.X  = X;     laplace_eqn.Y  = Y;
    laplace_eqn.Tx = Tx;    laplace_eqn.Ty = Ty;
end

% Save .mat file with partitions and thickness
% - Extract name of image
fname = outnm{IDX,:};

% - Compile path to save image
savepath = strcat(path,groupnm,'/',fname,'/');

% - Create new directory, if needed
if ~isdir(savepath);    mkdir(savepath);    end

% - Update filename for saving
if strcmpi(stain,'IF');    mfname = strrep(fname,rednm,mergenm);    else    mfname = fname;    end

% - Save .mat file
if strcmpi(tissue,'Myocardial Infarction         ');
    save(strcat(path,groupnm,'/',fname,'/',mfname,'_partition_',char(vartype),'.mat'),'localpart','laplace_eqn');
else
    save(strcat(path,groupnm,'/',fname,'/',mfname,'_partition_',char(vartype),'.mat'),'localpart');
end

% Plot overlay of partitions (and indices) on original image
% - Plot x,y coordinates of current partition
figure; imshow(I); hold on
for II = 1:npartq;
    for KK = 1:npartr;
        plot(localpart{II,KK}(:,2),localpart{II,KK}(:,1),'LineWidth',3,'Color',[0.0,0.0,0.0]);
    end
end

% - Generate labels for each partition location
for II = 1:npartq;
    for KK = 1:npartr;
        if strcmpi(vartype,'rad');
            midxi = ibound{II,KK}(:,2) == round(Nx0/2);    imin = min(ibound{II,KK}(midxi,1));
            midxo = obound{II,KK}(:,2) == round(Nx0/2);    omin = min(obound{II,KK}(midxo,1));
            text(round(Nx0/2),(imin+omin)/2,strcat('(',num2str(II),',',num2str(KK),')'),'Color','k','Fontsize',12,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
        elseif strcmpi(vartype,'circ') || strcmpi(vartype,'both');
            xmid = mean(localpart{II,KK}(:,2));    ymid = mean(localpart{II,KK}(:,1));
            text(xmid,ymid,strcat('(',num2str(II),',',num2str(KK),')'),'Color','k','Fontsize',12,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
end

% Extract color data from figure and write partitioned image to directory
if imsave;
    set(gcf,'color','w');    f = getframe(gcf);
    imwrite(f.cdata,strcat(path,groupnm,'/',fname,'/',mfname,'_partition_',char(vartype),'.tif'),'tif');
    close all force
end

fprintf('...Done! \n')


% _________________________ nested functions ______________________________
% -------------------------------------------------------------------------
    function [Xout,Yout,varargout] = points2contour(Xin,Yin,P,direction,varargin)
        
        %points2contour
        %Tristan Ursell
        %Sept 2013
        %
        %[Xout,Yout]=points2contour(Xin,Yin,P,direction)
        %[Xout,Yout]=points2contour(Xin,Yin,P,direction,dlim)
        %[Xout,Yout,orphans]=points2contour(Xin,Yin,P,direction,dlim)
        %[Xout,Yout,orphans,indout]=points2contour(Xin,Yin,P,direction,dlim)
        %
        %Given any list of 2D points (Xin,Yin), construct a singly connected
        %nearest-neighbor path in either the 'cw' or 'ccw' directions.  The code
        %has been written to handle square and hexagon grid points, as well as any
        %non-grid arrangement of points.
        %
        %'P' sets the point to begin looking for the contour from the original
        %ordering of (Xin,Yin), and 'direction' sets the direction of the contour,
        %with options 'cw' and 'ccw', specifying clockwise and counter-clockwise,
        %respectively.
        %
        %The optional input parameter 'dlim' sets a distance limit, if the distance
        %between a point and all other points is greater than or equal to 'dlim',
        %the point is left out of the contour.
        %
        %The optional output 'orphans' gives the indices of the original (Xin,Yin)
        %points that were not included in the contour.
        %
        %The optional output 'indout' is the order of indices that produces
        %Xin(indout)=Xout and Yin(indout)=Yout.
        %
        %There are many (Inf) situations where there is no unique mapping of points
        %into a connected contour -- e.g. any time there are more than 2 nearest
        %neighbor points, or in situations where the nearest neighbor matrix is
        %non-symmetric.  Picking a different P will result in a different contour.
        %Likewise, in cases where one point is far from its neighbors, it may be
        %orphaned, and only connected into the path at the end, giving strange
        %results.
        %
        %The input points can be of any numerical class.
        %
        %Note that this will *not* necessarily form the shortest path between all
        %the points -- that is the NP-Hard Traveling Salesman Problem, for which
        %there is no deterministic solution.  This will, however, find the shortest
        %path for points with a symmetric nearest neighbor matrix.
        %
        %see also: bwtraceboundary
        
        %check to make sure the vectors are the same length
        if length(Xin) ~= length(Yin);    error('Input vectors must be the same length.');    end
        
        %check to make sure point list is long enough
        if length(Xin) < 2;    error('The point list must have more than two elements.');    end
        
        %check distance limit
        if ~isempty(varargin)
            dlim = varargin{1};
            if dlim <= 0;    error('The distance limit parameter must be greater than zero.');    end
        else
            dlim = -1;
        end
        
        %check direction input
        if and(~strcmp(direction,'cw'),~strcmp(direction,'ccw'));    error(['Direction input: ' direction ' is not valid, must be either "cw" or "ccw".']);    end
        
        %check to make sure P is in the right range
        P = round(P);     npts = length(Xin);
        if or(P < 1,P > npts);    error('The starting point P is out of range.');    end
        
        %adjust input vectors for starting point
        if size(Xin,1) == 1
            Xin = circshift(Xin,[0,1-P]);
            Yin = circshift(Yin,[0,1-P]);
        else
            Xin = circshift(Xin,[1-P,0]);
            Yin = circshift(Yin,[1-P,0]);
        end
        
        %find distances between all points
        D = zeros(npts,npts);
        for q1 = 1:npts;     D(q1,:) = sqrt((Xin(q1)-Xin).^2 + (Yin(q1)-Yin).^2);    end
        
        %max distance and avoid self-connections
        maxD = max(D(:));    D = D + eye(npts)*maxD;
        
        %apply distance contraint by removing bad points and starting over
        if dlim > 0
            
            D(D >= dlim) = -1;
            
            %find bad points
            bad_pts = sum(D,1) == -npts;
            orphans = find(bad_pts);
            
            %check starting point
            if sum(orphans == P) > 0;    error('The starting point index is a distance outlier, choose a new starting point.');    end
            
            %get number of good points
            Xin = Xin(~bad_pts);    Yin = Yin(~bad_pts);    npts = length(Xin);
            
            %find distances between all points
            D = zeros(npts,npts);
            for q1 = 1:npts;    D(q1,:) = sqrt((Xin(q1)-Xin).^2 + (Yin(q1)-Yin).^2);    end
            
            %max distance and avoid self-connections
            maxD = max(D(:));    D=D+eye(npts)*maxD;
        else
            orphans = [];    bad_pts = zeros(size(Xin));
        end
        
        %tracking vector (has this original index been put into the ordered list?)
        track_vec = zeros(1,npts);
        
        %construct directed graph
        Xout = zeros(1,npts);       Xout(1) = Xin(1);
        Yout = zeros(1,npts);       Yout(1) = Yin(1);
        indout0 = zeros(1,npts);    indout0(1) = 1;
        
        p_now = 1;    track_vec(p_now) = 1;
        for q1 = 2:npts
            %get current row of distance matrix and remove used points
            curr_vec = D(p_now,:);    curr_vec(track_vec==1)=maxD;
            
            %find index of closest non-assigned point
            p_temp = find(curr_vec == min(curr_vec),1,'first');
            
            %reassign point
            Xout(q1) = Xin(p_temp);    Yout(q1) = Yin(p_temp);
            
            %move index
            p_now = p_temp;
            
            %update tracking
            track_vec(p_now) = 1;
            
            %update index vector
            indout0(q1) = p_now;
        end
        
        %undo the circshift
        temp1 = find(~bad_pts);    indout = circshift(temp1(indout0),[P,0]);
        
        %%%%%%% SET CONTOUR DIRECTION %%%%%%%%%%%%
        %contour direction is a *global* feature that cannot be determined until
        %all the points have been sequentially ordered.
        
        %calculate tangent vectors
        tan_vec = zeros(npts,3);
        for q1 = 1:npts
            if q1 == npts
                tan_vec(q1,:) = [Xout(1)-Xout(q1),Yout(1)-Yout(q1),0];
                tan_vec(q1,:) = tan_vec(q1,:)/norm(tan_vec(q1,:));
            else
                tan_vec(q1,:) = [Xout(q1+1)-Xout(q1),Yout(q1+1)-Yout(q1),0];
                tan_vec(q1,:) = tan_vec(q1,:)/norm(tan_vec(q1,:));
            end
        end
        
        %determine direction of contour
        local_cross = zeros(1,npts);
        for q1 = 1:npts
            if q1 == npts
                cross1 = cross(tan_vec(q1,:),tan_vec(1,:));
            else
                cross1 = cross(tan_vec(q1,:),tan_vec(q1+1,:));
            end
            local_cross(q1) = asin(cross1(3));
        end
        
        %figure out current direction
        if sum(local_cross) < 0;    curr_dir = 'cw';    else    curr_dir = 'ccw';    end
        
        %set direction of the contour
        if and(strcmp(curr_dir,'cw'),strcmp(direction,'ccw'))
            Xout = fliplr(Xout);    Yout = fliplr(Yout);
        end
        
        %varargout
        if nargout == 3;    varargout{1} = orphans;    end
        if nargout == 4;    varargout{1} = orphans;    varargout{2} = indout;    end
        
    end


    function [cn,on] = inpoly(p,node,edge,TOL)
        
        %  INPOLY: Point-in-polygon testing.
        %
        % Determine whether a series of points lie within the bounds of a polygon
        % in the 2D plane. General non-convex, multiply-connected polygonal
        % regions can be handled.
        %
        % SHORT SYNTAX:
        %
        %   in = inpoly(p,node);
        %
        %   p   : The points to be tested as an Nx2 array [x1 y1; x2 y2; etc].
        %   node: The vertices of the polygon as an Mx2 array [X1 Y1; X2 Y2; etc].
        %         The standard syntax assumes that the vertices are specified in
        %         consecutive order.
        %
        %   in  : An Nx1 logical array with IN(i) = TRUE if P(i,:) lies within the
        %         region.
        %
        % LONG SYNTAX:
        %
        %  [in,on] = inpoly(p,node,edge);
        %
        %  edge: An Mx2 array of polygon edges, specified as connections between
        %        the vertices in NODE: [n1 n2; n3 n4; etc]. The vertices in NODE
        %        do not need to be specified in connsecutive order when using the
        %        extended syntax.
        %
        %  on  : An Nx1 logical array with ON(i) = TRUE if P(i,:) lies on a
        %        polygon edge. (A tolerance is used to deal with numerical
        %        precision, so that points within a distance of
        %        eps^0.8*norm(node(:),inf) from a polygon edge are considered "on"
        %        the edge.
        %
        % EXAMPLE:
        %
        %   polydemo;       % Will run a few examples
        %
        % See also INPOLYGON
        
        % The algorithm is based on the crossing number test, which counts the
        % number of times a line that extends from each point past the right-most
        % region of the polygon intersects with a polygon edge. Points with odd
        % counts are inside. A simple implementation of this method requires each
        % wall intersection be checked for each point, resulting in an O(N*M)
        % operation count.
        %
        % This implementation does better in 2 ways:
        %
        %   1. The test points are sorted by y-value and a binary search is used to
        %      find the first point in the list that has a chance of intersecting
        %      with a given wall. The sorted list is also used to determine when we
        %      have reached the last point in the list that has a chance of
        %      intersection. This means that in general only a small portion of
        %      points are checked for each wall, rather than the whole set.
        %
        %   2. The intersection test is simplified by first checking against the
        %      bounding box for a given wall segment. Checking against the bbox is
        %      an inexpensive alternative to the full intersection test and allows
        %      us to take a number of shortcuts, minimising the number of times the
        %      full test needs to be done.
        %
        %   Darren Engwirda: 2005-2007
        %   Email          : d_engwirda@hotmail.com
        %   Last updated   : 23/11/2007 with MATLAB 7.0
        %
        % Problems or suggestions? Email me.
        
        % ERROR CHECKING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin<4
            TOL = 1.0e-12;
            if nargin<3
                edge = [];
                if nargin<2
                    error('Insufficient inputs');
                end
            end
        end
        nnode = size(node,1);
        if isempty(edge)                                                           % Build edge if not passed
            edge = [(1:nnode-1)' (2:nnode)'; nnode 1];
        end
        if size(p,2)~=2
            error('P must be an Nx2 array.');
        end
        if size(node,2)~=2
            error('NODE must be an Mx2 array.');
        end
        if size(edge,2)~=2
            error('EDGE must be an Mx2 array.');
        end
        if max(edge(:))>nnode || any(edge(:)<1)
            error('Invalid EDGE.');
        end
        
        % PRE-PROCESSING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n  = size(p,1);
        nc = size(edge,1);
        
        % Choose the direction with the biggest range as the "y-coordinate" for the
        % test. This should ensure that the sorting is done along the best
        % direction for long and skinny problems wrt either the x or y axes.
        dxy = max(p,[],1)-min(p,[],1);
        if dxy(1)>dxy(2)
            % Flip co-ords if x range is bigger
            p = p(:,[2,1]);
            node = node(:,[2,1]);
        end
        tol = TOL*min(dxy);
        
        % Sort test points by y-value
        [y,i] = sort(p(:,2));
        x = p(i,1);
        
        % MAIN LOOP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cn = false(n,1);     % Because we're dealing with mod(cn,2) we don't have
        % to actually increment the crossing number, we can
        % just flip a logical at each intersection (faster!)
        on = cn;
        for k = 1:nc         % Loop through edges
            
            % Nodes in current edge
            n1 = edge(k,1);
            n2 = edge(k,2);
            
            % Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
            %           - also get xmin = min(x1,x2), xmax = max(x1,x2)
            y1 = node(n1,2);
            y2 = node(n2,2);
            if y1<y2
                x1 = node(n1,1);
                x2 = node(n2,1);
            else
                yt = y1;
                y1 = y2;
                y2 = yt;
                x1 = node(n2,1);
                x2 = node(n1,1);
            end
            if x1>x2
                xmin = x2;
                xmax = x1;
            else
                xmin = x1;
                xmax = x2;
            end
            
            % Binary search to find first point with y<=y1 for current edge
            if y(1)>=y1
                start = 1;
            elseif y(n)<y1
                start = n+1;
            else
                lower = 1;
                upper = n;
                for j = 1:n
                    start = round(0.5*(lower+upper));
                    if y(start)<y1
                        lower = start;
                    elseif y(start-1)<y1
                        break;
                    else
                        upper = start;
                    end
                end
            end
            
            % Loop through points
            for j = start:n
                % Check the bounding-box for the edge before doing the intersection
                % test. Take shortcuts wherever possible!
                
                YY = y(j);   % Do the array look-up once & make a temp scalar
                if YY<=y2
                    XX = x(j);   % Do the array look-up once & make a temp scalar
                    if XX>=xmin
                        if XX<=xmax
                            
                            % Check if we're "on" the edge
                            on(j) = on(j) || (abs((y2-YY)*(x1-XX)-(y1-YY)*(x2-XX))<tol);
                            
                            % Do the actual intersection test
                            if (YY<y2) && ((y2-y1)*(XX-x1)<(YY-y1)*(x2-x1))
                                cn(j) = ~cn(j);
                            end
                            
                        end
                    elseif YY<y2   % Deal with points exactly at vertices
                        % Has to cross edge
                        cn(j) = ~cn(j);
                    end
                else
                    % Due to the sorting, no points with >y
                    % value need to be checked
                    break
                end
            end
            
        end
        
        % Re-index to undo the sorting
        cn(i) = cn|on;
        on(i) = on;
        
    end      % inpoly()


end