function thick = thickness(IDX, Ibinary, bprops)

global path groupnm outnm

global stain imsave ext

% Make sure that all libraries needed for fast marching have been compiled and added correctly
if ~exist(strcat('histology-analysis/fm/perform_front_propagation_2d.mexw64'),'file');
    
    % Change current directory
    cd histology-analysis/fm/
    
    % Setup C++ compiler, compile mex files, and reset directory
    mex -setup c++
    mex ./mex/perform_front_propagation_2d.cpp ./mex/perform_front_propagation_2d_mex.cpp ./mex/fib.cpp
    cd ../..
    
end

% Extract boundary properties (from boundaries.m)
ibound = bprops.ibound;    obound = bprops.obound;    centpt = bprops.centpt;

% ----- Thickness calculation ----------------------------------------------
% Invert tissue mask and extract dimensions
phi0 = imcomplement(Ibinary);    [Nx0,Ny0] = size(phi0);

% Create mask to restrict propagation outside of tissue geometry (inside: +Inf, outside: -Inf)
M = zeros(Nx0,Ny0) - Inf;    M(phi0 == 1) = +Inf;    options.constraint_map = M;


% Compute OUTWARD trajectory length, L0
% - Initialize front propagation for fast marching
start_points = ibound;    W = ones(Nx0,Ny0);

% - Execute fast marching
[L0,~,~] = perform_fast_marching(W, start_points', options);


% Compute INWARD trajectory length, L1
% - Initialize front propagation for fast marching
start_points = obound;    W = ones(Nx0,Ny0);

% - Execute fast marching
[L1,~,~] = perform_fast_marching(W, start_points', options);


% Evaluate thickness as W = L0 + L1 (Yezzi and Prince 2003, Rocha et al. 2007)
W  = Nx0.*(L0 + L1);    % [pixels]
Wb =  L0./(L0 + L1);    % Normalized: 0 - inner, 1 - outer boundary

% Compute skeleton of tissue section to approximate midline
sz = 5e-2;
Wmid = abs(Wb - 0.5) < sz;    Wmid = bwmorph(Wmid,'thin',Inf);

% Check for branching midline
Wmid = bwmorph(Wmid,'spur',Inf);

% branch = bwmorph(Wmid,'branchpoints');    branch = sum(branch(:));
% while branch
%     if abs(sz - 1e-2) < 1e-6;    sz = sz/10;    else    sz = sz - 1e-2;    end
%     Wmid = abs(Wb - 0.5) < sz;                  Wmid = bwmorph(Wmid,'thin',Inf);
%     branch = bwmorph(Wmid,'branchpoints');    branch = sum(branch(:));
% end

% Get pixel coordinates of midline
Widx = find(Wmid);    [Wx,Wy] = ind2sub([Nx0,Ny0],Widx);    Wbnd = cat(2,Wx,Wy);


% Order points extracted from midline
[Wx,Wy] = points2contour(round(Wbnd(:,1)),round(Wbnd(:,2)),1,'cw');    Wbnd = cat(2,Wx',Wy');


% Smooth midline boundary
smoothval = 20;    Wbsm = smooth2(cat(1,Wbnd(end-(smoothval/2)+1:end,:),Wbnd,Wbnd(1:(smoothval/2),:)),smoothval,0);

% Reinitialize storage for image
clear Wmid;    Wmid = false(Nx0,Ny0);

% Plot white contours on black background
for II = 1:length(Wbsm(:,1));    Wmid(round(Wbsm(II,1)),round(Wbsm(II,2))) = 1;    end

% Connect smoothed ends of inner/outer contours
Wmid = linept(Wmid,round(Wbsm(1,1)),round(Wbsm(1,2)),round(Wbsm(end,1)),round(Wbsm(end,2)));


% Get pixel coordinates of smoothed midline
[Wx,Wy] = points2contour(round(Wbsm(:,1)),round(Wbsm(:,2)),1,'cw');    Wbnd = cat(2,Wx',Wy');    Wbnd = flipud(Wbnd);


% Re-order theta orientation to have 0deg on +x-axis
% - Find start of midline contour
for II = centpt(1)-0.5:Ny0;
    if Wmid(centpt(2)-0.5,II) > 0 && Wmid(centpt(2)-0.5,II-1) == 0;    mpt = [II centpt(2)-0.5];    end
end

% - Find location of start in boundary array and shift accordingly
[~,idx] = ismember(fliplr(mpt),Wbnd,'rows');    Wbnd = circshift(Wbnd,[-idx+1,0]);   


% Extract smoothed midline thickness
Wavg = Wmid.*W;

% Extract final re-ordered thickness values
thavg = NaN(length(Wbnd(:,1)),1);
for II = 1:length(Wbnd(:,1));    thavg(II,1) = Wavg(Wbnd(II,1),Wbnd(II,2));    end
thavg = cat(1,thavg,thavg(1));

% ----- File/Image saving -------------------------------------------------
% Extract name of image
fname = outnm{IDX,:};

% Compile path to save image
savepath = strcat(path,groupnm,'/',fname,'/');

% Create new directory, if needed
if ~isdir(savepath);    mkdir(savepath);    end

% Update filename for saving
if strcmpi(stain,'IF');    mfname = strrep(fname,rednm,mergenm);    else    mfname = fname;    end

% Define circumferential coordinate array
theta = linspace(0,360,length(thavg));

if imsave
    % Set parameters for figue sizing
    scrsz = get(0,'ScreenSize');
    figpos = [0.25*scrsz(3) 0.04*scrsz(4) 0.92*scrsz(4) (Nx0/Ny0)*0.92*scrsz(4)];
    
    % Save thickness map figure
    figure('Position',figpos);  f = imagesc(W);  axis('image');  set(f,'AlphaData',isfinite(W));
    hold on;
    plot(ibound(:,2),ibound(:,1),'k','LineWidth',2);  plot(obound(:,2),obound(:,1),'k','LineWidth',2);
    plot(Wbnd(:,2),Wbnd(:,1),'k--','LineWidth',2);    plot(mpt(1),mpt(2),'ro','MarkerFaceColor','r','MarkerSize',13);
    set(gca,'visible','off');  colorbar;  set(gca,'FontSize',16)
    
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',strcat(path,groupnm,'/',fname,'/',mfname,'_thickness-map.',ext));
    close(gcf);
    
    % Plot and save thickness variation line plot
    figure; plot(theta,thavg,'k','LineWidth',2);  hold on;
    plot(theta(1),thavg(1),'ro','MarkerFaceColor','r','MarkerSize',7);
    plot(theta(end),thavg(end),'ro','MarkerFaceColor','r','MarkerSize',7);
    xlim([-2,362]);  xlabel('Circumferential Position [deg]');  ylabel('Tissue Thickness [pix]');  title('Thickness Variation');
    ax = gca;  ax.XTick = [0,90,180,270,360];  ax.XTickLabels = {'0','90','180','270','360'};  set(gca,'FontSize',13)
    
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',strcat(path,groupnm,'/',fname,'/',mfname,'_thickness-var.',ext));
    close(gcf);
end

% Save .mat file of thickness variation
clear Wavg;    Wavg = cat(2,theta',thavg);

thick.map = W;                 thick.norm = Wb;           thick.avg = Wavg;
thick.bound.inner = ibound;    thick.bound.mid = Wbnd;    thick.bound.outer = obound;
save(strcat(path,groupnm,'/',fname,'/',mfname,'_thickness.mat'),'thick');

close all force


fprintf('...Done! \n')


% _________________________ nested functions ______________________________
% -------------------------------------------------------------------------
    function result = linept(matrice, X1, Y1, X2, Y2)
        % Connect two pixels in a matrice with 1
        %
        % Command line
        % ------------
        % result=linept(matrice, X1, Y1, X2, Y2)
        %   matrice : matrice where I'll write
        %   (X1, Y1), (X2, Y2) : points to connect
        %   result : matrix + the line
        %
        % Georges Cubas 20/11/03
        % georges.c@netcourrier.com
        % Version 1.0
        
        result = matrice;
        
        for xx = max(1,X1):sign(X2 - X1):max(1,X2)
            yy = round(f(xx, X1, Y1, X2, Y2));
            if yy > 0;    result(xx,yy) = 1;    end
        end
        
        for yy = max(1,Y1):sign(Y2 - Y1):max(1,Y2)
            xx = round(f2(yy, X1, Y1, X2, Y2));
            if xx > 0;    result(xx,yy) = 1;     end
        end
        
        function y = f(x, X1, Y1, X2, Y2)
            a = (Y2 - Y1)/(X2 - X1);
            b = Y1 - X1 * a;
            y = a * x + b;
        end
        
        function x = f2(y, X1, Y1, X2, Y2)
            if X1 == X2
                x = X1;
            else
                a = (Y2 - Y1)/(X2 - X1);
                b = Y1 - X1 * a;
                x = (y - b)/a;
            end
        end
    end


    function matrixOut = smooth2(matrixIn, Nr, Nc)
        
        %			MATRIXOUT=SMOOTH2(MATRIXIN,Nr,Nc) smooths the data in MATRIXIN
        %           using a running mean over 2*N+1 successive points, N points on
        %           each side of the current point.  At the ends of the series
        %           skewed or one-sided means are used.
        %
        %           Inputs: matrixIn - original matrix
        %                   Nr - number of points used to smooth rows
        %                   Nc - number of points to smooth columns
        %
        %           Outputs:matrixOut - smoothed version of original matrix
        %
        %           Remark: By default, if Nc is omitted, Nc = Nr.
        %
        %           Written by Kelly Hilands, October 2004
        %           Applied Research Laboratory
        %           Penn State University
        
        %Initial error statements and definitions
        if nargin<2, error('Not enough input arguments!'), end
        
        N(1) = Nr;
        if nargin < 3;    N(2) = N(1);    else    N(2) = Nc;    end
        
        if length(N(1))~=1, error('Nr must be a scalar!'), end
        if length(N(2))~=1, error('Nc must be a scalar!'), end
        
        [row,col] = size(matrixIn);
        eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
        eL = eL./(repmat(sum(eL,2),1,row));
        eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);
        eR = eR./(repmat(sum(eR,1),col,1));
        
        matrixOut = eL*matrixIn*eR;
        
    end


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

end