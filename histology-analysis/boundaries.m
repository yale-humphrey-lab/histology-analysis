function [Ibw, bprops] = boundaries(I,LV)

global stain SF tissue

% Get image size
[x,y,~] = size(I);

% Extract image scale factors
SFr = SF(1);    SFm = SF(2);

% Convert to binary and perform successive open operations to smooth object
if strcmpi(stain,'dPSR') || strcmpi(stain,'POL') || strcmpi(stain,'IF')
    
    % Preprocess and pass disk along edge of binary object to smooth
    IBW = bwareaopen(imcomplement(bwareaopen(imopen(imcomplement(im2bw(I,0.01)),strel('disk',round(15*SFr*SFm))),1000)),1000);
    
    % Force image to have 2-pixel border
    IBW([1:2 end-1:end],:) = 0;     IBW(:,[1:2 end-1:end]) = 0;
    
    % Inverted binary mask
    IBWi = imcomplement(IBW);
else
    % Preprocessing of binary image
    IBW = bwareaopen(imcomplement(bwareaopen(im2bw(I,0.99),5000)),5000);
    
    % Pass disk along edge of binary object to smooth
    IBW = imclose(IBW,strel('disk',round(10*SFr*SFm)));
    
    % Force image to have 2-pixel border
    IBW([1:2 end-1:end],:) = 0;     IBW(:,[1:2 end-1:end]) = 0;
    
    % Inverted binary mask
    IBWi = imcomplement(IBW);
end


% Compute location of centroid of cross-section
cc = bwconncomp(IBWi,8);
S = regionprops(cc,'Area','Centroid');

% Initialize and write pixel areas from structure to array
bwarea = zeros(length(S(:,1)),1);
for II = 1:length(S(:,1));    bwarea(II,1) = S(II,1).Area;    end

% Find correct location for boundary identification
if strcmpi(tissue,'Myocardial Infarction         ');
    if length(bwarea) > 2;
        % Indices of two largest areas (should be left and right ventricle)
        idx1 = find(bwarea == max(bwarea(bwarea < max(bwarea))));                               d1 = pdist(cat(1,S(idx1).Centroid,LV),'euclidean');
        idx2 = find(bwarea == max(bwarea(bwarea < max(bwarea) & bwarea ~= bwarea(idx1))));      d2 = pdist(cat(1,S(idx2).Centroid,LV),'euclidean');
        
        % Index of left ventricle
        if d2 < d1;    idx = idx2;    else    idx = idx1;    end
    else
        % Index of second largest area (should be lumen of vessel)
        idx = find(bwarea == max(bwarea(bwarea < max(bwarea))));
    end
else
    % Index of second largest area (should be lumen of vessel)
    idx = find(bwarea == max(bwarea(bwarea < max(bwarea))));
end

% Define origin as center of sample
cent = cell2mat(struct2cell(S)');       orig = cent(idx,2:3);       orig = floor(orig)+0.5;

% March along image columns to find start of inner and outer contours (b -> w and w -> b)
cols = size(IBW,2);     opt0 = NaN(1,2);
for II = orig(1)-0.5:cols;
    if IBWi(orig(2)-0.5,II) == 1 && IBWi(orig(2)-0.5,II-1) == 0;        opt0(end+1,:) = [II orig(2)-0.5];        end
end

% Extract correct outer point, if multiple were identified
opt0 = opt0(2:end,:);     opt0 = opt0(end,:);

% Trace inner and outer boundaries of raw object
obound = bwtraceboundary(IBWi,[opt0(2), opt0(1)], 'N', 8);


% Build mask of outer boundary to correct location of origin
maskO = false(size(IBW,1),size(IBW,2));    
for II = 1:length(obound(:,1));    maskO(obound(II,1),obound(II,2)) = 1;    end

% Apply mask to binary image to keep only inner region
IBWiO = IBWi;    IBWiO(maskO == 1) = 0;    maskO = imfill(maskO,'holes');    IBWiO(maskO == 0) = 0;

% Compute regional maxima of distance transform of inner region
IBWd = bwdist(~IBWiO);    IBWd = imregionalmax(IBWd,8);
            
% Find closest point on distance transform skeleton
idx = find(IBWd);    xd = NaN(length(idx),1);    yd = NaN(length(idx),1);    d = NaN(length(idx),1);
for II = 1:length(idx);    
    [xd(II),yd(II)] = ind2sub(size(IBW),idx(II));    d(II) = sqrt((orig(2)-xd(II))^2 + (orig(1)-yd(II))^2);    
end

% Re-assign origin to be on skeleton
[~,didx] = min(d);    orig = [yd(didx), xd(didx)];    orig = floor(orig)+0.5;    centpt = orig;


% March along image columns to find start of inner and outer contours (b -> w and w -> b)
cols = size(IBW,2);     ipt = NaN(1,2);     opt = NaN(1,2);
for II = orig(1)-0.5:cols;
    if IBWi(orig(2)-0.5,II) == 1 && IBWi(orig(2)-0.5,II-1) == 0;        opt(end+1,:) = [II orig(2)-0.5];        end
    if IBWi(orig(2)-0.5,II) == 0 && IBWi(orig(2)-0.5,II-1) == 1;        ipt(end+1,:) = [II-1 orig(2)-0.5];      end
end

% Extract correct outer point, if multiple were identified
ipt = ipt(2:end,:);     ipt = ipt(1,:);
opt = opt(2:end,:);     opt = opt(end,:);

% Trace inner and outer boundaries of raw object
ibound = bwtraceboundary(IBWi,[ipt(2), ipt(1)], 'N', 8);
obound = bwtraceboundary(IBWi,[opt(2), opt(1)], 'N', 8);


% 2D smoothing algorithm to smooth boundary contours
smoothval = 50;
smicont = smooth2(ibound,smoothval,0);
smocont = smooth2(obound,smoothval,0);

% Find max and min values of the endpoints of each contour to connect the smoothed points
ymini = round(min([smicont(1,1) smicont(end,1)]));
ymaxi = round(max([smicont(1,1) smicont(end,1)]));
ymino = round(min([smocont(1,1) smocont(end,1)]));
ymaxo = round(max([smocont(1,1) smocont(end,1)]));

% Linear fit to close smoothed contours
Pi = polyfit([smicont(1,1) smicont(end,1)],[smicont(1,2) smicont(end,2)],1);        Fi = polyval(Pi,(ymini:ymaxi));
Po = polyfit([smocont(1,1) smocont(end,1)],[smocont(1,2) smocont(end,2)],1);        Fo = polyval(Po,(ymino:ymaxo));


% Reinitialize storage for image
clear I;    I = zeros(x,y);

% Plot white contours on black background
for II = 1:length(smicont(:,1));        I(round(smicont(II,1)),round(smicont(II,2))) = 1;       end
for II = 1:length(smocont(:,1));        I(round(smocont(II,1)),round(smocont(II,2))) = 1;       end

% Connect smoothed ends of inner/outer contours
I = linept(I,round(smicont(1,1)),round(smicont(1,2)),round(smicont(end,1)),round(smicont(end,2)));
I = linept(I,round(smocont(1,1)),round(smocont(1,2)),round(smocont(end,1)),round(smocont(end,2)));


% Convert to binary image and fill interior of contours using flood-filling algorithm
IBW = im2bw(I,0.99);        IBW = imfill(IBW,[round((ipt(2)+opt(2))/2) round((ipt(1)+opt(1))/2)]);

% Final binary image mask
Ibw = imcomplement(IBW);


% % Re-trace contours using new centroid of filled region
% S = regionprops(Ibw,'Centroid');        centpt = S(2,1).Centroid;       centpt = floor(centpt)+0.5;

% March along image columns to find start of inner and outer contours (b -> w and w -> b)
for II = centpt(1)-0.5:cols;
    if Ibw(centpt(2)-0.5,II) == 1 && Ibw(centpt(2)-0.5,II-1) == 0;      opt = [II centpt(2)-0.5];       end
    if Ibw(centpt(2)-0.5,II) == 0 && Ibw(centpt(2)-0.5,II-1) == 1;      ipt = [II-1 centpt(2)-0.5];     end
end


% Trace inner and outer boundaries of final processed object
bprops.centpt = centpt;
bprops.ibound = bwtraceboundary(Ibw,[ipt(2), ipt(1)], 'N', 8);
bprops.obound = bwtraceboundary(Ibw,[opt(2), opt(1)], 'N', 8);


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
        if nargin<3
            N(2) = N(1);
        else
            N(2) = Nc;
        end
        
        if length(N(1))~=1, error('Nr must be a scalar!'), end
        if length(N(2))~=1, error('Nc must be a scalar!'), end
        
        [row,col] = size(matrixIn);
        eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
        eL = eL./(repmat(sum(eL,2),1,row));
        eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);
        eR = eR./(repmat(sum(eR,1),col,1));
        
        matrixOut = eL*matrixIn*eR;
        
    end

end