function X = HSLfilter(Hpar, Spar, Lpar, H, A, color, bground)

global stain layernm

% Transform HSL image back to RGB colorspace for coupled color comparison
RGB = colorspace('HSL->RGB',H);

% Define different colors in HSL colorspace for plotting constituent images
if strcmpi(color,'w');
    hslcode = [360; 0; 1];
    
elseif strcmpi(color,'r');
    hslcode = [7; 1; 0.5];
    
elseif strcmpi(color,'g');
    hslcode = [120; 1; 0.25];
    
elseif strcmpi(color,'b');
    hslcode = [240; 1; 0.4];
    
elseif strcmpi(color,'o');
    hslcode = [30; 1; 0.5];
    
elseif strcmpi(color,'c');
    hslcode = [180; 1; 0.5];
    
elseif strcmpi(color,'m');
    hslcode = [310; 1; 0.5];
    
elseif strcmpi(color,'y');
    hslcode = [60; 1; 0.5];
    
elseif strcmpi(color,'k');
    hslcode = [360; 0; 0];
    
end

% Define upper and lower bounds on H- S- and L- parameters
eps = 1e-8;     % Small offset to correct for rounding error in colorspace.m

Hlow = Hpar(1) - eps;
Hup  = Hpar(2) + eps;

Slow = Spar(1) - eps;
Sup  = Spar(2) + eps;

Llow = Lpar(1) - eps;
Lup  = Lpar(2) + eps;

% Hard-coded L- and RGB threshold for extracting BLACK pixels
Lblk = 0.30 + eps;

if (strcmpi(stain,'VVG') || strcmpi(stain,'MOV')) && Lup <= Lblk
    
    if strcmpi(layernm,'adventitia')
        RGBblk = 0.5;
    else
        % Select RGB threshold based on G-channel of RGB image
        [cntg,xr] = imhist(RGB(:));
        
        cntgsm = sgolayfilt(cntg(2:end-1),3,41);
        cntgdf = sgolayfilt(diff(cntgsm),3,41);
        
        cross = find(abs(diff(sign(cntgdf))) == 2); % each zero-crossing
        
        % Use color bins for data fitting of normalized histogram data
        if length(cross) >= 2
            if (cross(2) - cross(1)) > 25
                frac = 0.75;
            else
                frac = 0.50;
            end
            
            binnm = round(cross(1) + frac*(cross(2) - cross(1)));
        else
            binnm = 50;
        end
        
        RGBblk = NaN(10,1);
        for k = 1:10;
                        
            options = optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-8,'Display','off');
            x0 = [rand rand rand rand rand];
            xd = 256.*xr(2:binnm);
            
            cntgn = cntg(2:end);    ind = find(cntgn > 1);
            
            yd = cntg(2:binnm)./cntg(ind(1)+1);
            x = lsqcurvefit(@myfun,x0,xd,yd,[],[],options);
            
            % Select threshold based on steady-state/final value of biexponential function
            yf = myfun(x,xd);
            val = yf(end) + 0.001;        % Slight increment to improve thresholding
            
            [~,blkth] = min(abs(yf-val));
            RGBblk(k,1) = (blkth + eps)/255;
            
        end
        
        RGBsort = sort(RGBblk,'descend');
        for i = 1:length(RGBsort)
            RGBct(i,1) = sum(RGBblk == RGBsort(i));
        end
        RGBcsort = flipud(unique(sort(RGBct,'descend')));
        
        val1 = RGBsort(RGBct == RGBcsort(1));
        
        if length(RGBcsort) > 1
            val2 = RGBsort(RGBct == RGBcsort(2));
            RGBblk = 5*mean([val1(1) val2(1)]);
        else
            RGBblk = 5*val1(1);
        end
    end
        
end

[x,y,z] = size(H);

X = ones(x,y,z);

% Iterate through all pixels and check HSL parameters
% -- Note: If isempty(A) ~= 1 then the current consituent is NOT the first constituent analyzed in the image.
%          This helps to ensure that there is not any pixel overlap in the analysis
if strcmpi(bground,'black');
    
    % Split input (H) and output (X) images into individual channels for color thresholding
    H1 = H(:,:,1);      H2 = H(:,:,2);      H3 = H(:,:,3);
    X1 = X(:,:,1);      X2 = X(:,:,2);      X3 = X(:,:,3);
    
    % Generate masks of all pixels within the prescribed HSL thresholds
    if isempty(A) && Hlow < Hup
        mask = (H3 == 1) | (H1 >= Hlow & H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);

    elseif isempty(A) && Hlow > Hup
        mask = (H3 == 1) | (H1 >= Hlow | H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);
        
    elseif ~isempty(A) && Hlow < Hup
        mask = (A(:,:,3) == 1) & (H1 >= Hlow & H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);
       
    elseif ~isempty(A) && Hlow > Hup
        mask = (A(:,:,3) == 1) & (H1 >= Hlow | H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);
    end
    
    % Pseudocolor all positive/true pixels, force negative/false pixels to white (i.e., background)
    X1(mask)  = hslcode(1);    X2(mask)  = hslcode(2);    X3(mask)  = hslcode(3);
    X1(~mask) = 360;           X2(~mask) = 0;             X3(~mask) = 1;
    
    % Recompile thresholded image
    X = cat(3,X1,X2,X3);

    
elseif strcmpi(bground,'white');
    
    % Split input (H) and output (X) images into individual channels for color thresholding
    H1 = H(:,:,1);      H2 = H(:,:,2);      H3 = H(:,:,3);
    X1 = X(:,:,1);      X2 = X(:,:,2);      X3 = X(:,:,3);  
    
    % Generate masks of all pixels within the prescribed HSL thresholds
    if isempty(A) && Hlow < Hup
        mask = (H1 >= Hlow & H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);
        
    elseif isempty(A) && Hlow > Hup
        mask = (H1 >= Hlow | H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);
        
    elseif ~isempty(A) && Hlow < Hup
        mask = (A(:,:,3) == 1) & (H1 >= Hlow & H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);
        
    elseif ~isempty(A) && Hlow > Hup
        mask = (A(:,:,3) == 1) & (H1 >= Hlow | H1 <= Hup) & (H2 >= Slow & H2 <= Sup) & (H3 >= Llow & H3 <= Lup);
    end
        
    % Additional masking is needed when extracting black pixels...
    if Lup <= Lblk
        
        % Generate mask of average RGB values compared to black threshold
        bkmk = (mean(RGB,3) <= RGBblk);
        
        % Both masks must be positive/true to accept pixels
        X1(mask &  bkmk) = hslcode(1);    X2(mask &  bkmk) = hslcode(2);    X3(mask &  bkmk) = hslcode(3);
        X1(mask & ~bkmk) = 360;           X2(mask & ~bkmk) = 0;             X3(mask & ~bkmk) = 1;
        
    else
        % Accept all positive/true pixels in mask
        X1(mask) = hslcode(1);    X2(mask) = hslcode(2);    X3(mask) = hslcode(3);
    end
    
    % Force all negative/false pixels to white (i.e., background)
    X1(~mask) = 360;    X2(~mask) = 0;    X3(~mask) = 1;
    
    % Recompile thresholded image
    X = cat(3,X1,X2,X3);
    
end

% ---------------------------- Nested Functions ---------------------------

% Biexponential Function
    function F = myfun(x,xd)
        F = x(1).*exp(-x(2).*xd) + x(3).*exp(-x(4).*xd) + x(5);
    end

end
