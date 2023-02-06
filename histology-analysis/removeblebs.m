function Efin = removeblebs(E,medadv)

global SF

SFr = SF(1);    SFm = SF(2);

% Invert lightness channel of identified elastin
Ei = imcomplement(im2bw(E(:,:,3)));

% Compute euclidean distance transform
Ed = bwdist(imcomplement(Ei),'euclidean');

% Morphological skeleton
Esk = bwmorph(Ed,'skel',Inf);

% Dilate skeleton by one pixel and recompute skeleton to merge neighboring endpoints
Esk = imdilate(Esk,strel('disk',1));
Esk = bwmorph(Esk,'skel',Inf);

% Identify pixel location of branchpoints in skeletonized image
Eb = bwmorph(Esk,'branchpoints');
[bpti,bptj] = find(Eb);

% Identify pixel location of endpoints in skeletonized image
Ee = bwmorph(Esk,'endpoints');
[epti,eptj] = find(Ee);

% Sizing parameters for analysis window and branch removal threshold
wsz = round((12*SFm + 1)*SFr);     bth = 12*SFm*SFr;    sz = 2*wsz+1;     mpt = wsz+1;

% Initialize branch removal parameters
Esh = Esk;     anyept = 1;    final = 0;    num = 1;

% Remove branches as long as unique endpoints are available
while anyept
    
    h = waitbar(0,'Pruning Branches...','Name',['Iteration ' num2str(num)]);
    
    % Check endpoints...
    for i = 1:length(epti);
        waitbar(i/length(epti))
        
        % Intialize branch and analysis window
        Dmsk = zeros(sz,sz);
        Eshw = Esh(epti(i)-wsz:epti(i)+wsz,eptj(i)-wsz:eptj(i)+wsz);
        
        % Compute all brancpoints within analysis window
        Ebw = bwmorph(Eshw,'branchpoints');
        
        % Geodesic distance transform with current endpoint as seed point
        D = bwdistgeodesic(Eshw,mpt,mpt);
        
        % Distance to closest branchpoint
        bpd = min(D(Ebw > 0));
        
        % If branchpoint exists and is closer than threshold, remove branch
        if ~isempty(bpd) && bpd <= bth;
            Dmsk(D < bpd) = 1;
            
            Esh(epti(i)-wsz:epti(i)+wsz,eptj(i)-wsz:eptj(i)+wsz) = Eshw-Dmsk;
        end
        
    end
    
    if ~final
        
        if num == 1
            % Remove branchpoints from image to reduce branch removal aritfacts
            Esh = Esh - Eb;
            
            % Correct for potential holes from branchpoint removal
            Esh = imdilate(Esh,strel('disk',1));
            Esh = bwmorph(Esh,'skel',Inf);
        end
        
        % Recompute endpoints of cleaned image
        Ee = bwmorph(Esh,'endpoints');
        [eptin,eptjn] = find(Ee);
        
        % Check if same number of endpoints as previous iteration
        ept = [epti eptj];
        eptn = [eptin eptjn];
        
        if length(eptn(:,1)) == length(ept(:,1));
            anyept = 0;    final = 1;
        else
            anyept = 1;    epti = eptin;    eptj = eptjn;
        end
        
    end
    
    num = num + 1;
    
    close(h)
    
end

% Remove branchpoints from image to reduce branch removal artifacts
Esh = Esh - Eb;

% Correct for potential holes from branchpoint removal
Esh = imdilate(Esh,strel('disk',2));
Esh = bwmorph(Esh,'skel',Inf);

% Remove branchpoints that are no longer contacting (or neighboring) skeleton axis
Eshd = bwdist(Esh,'euclidean');
for i = 1:length(bpti);
    if Eshd(bpti(i),bptj(i)) >= 2;
        bpti(i) = NaN;
        bptj(i) = NaN;
    end
end

bpti = bpti(isfinite(bpti));
bptj = bptj(isfinite(bptj));

% Associate euclidean distance with cleaned skeleton axis
[xloc,yloc] = find(Esh);     locsh = [xloc yloc];

Edn = double(Esh);

wsz = 15; mpt = wsz + 1;
h = waitbar(0,'Mapping Distances to Skeleton...');
for i = 1:length(locsh(:,1));
    waitbar(i/length(locsh(:,1)));
    
    Edw = Ed(locsh(i,1)-wsz:locsh(i,1)+wsz,locsh(i,2)-wsz:locsh(i,2)+wsz);
    Eskw = Esk(locsh(i,1)-wsz:locsh(i,1)+wsz,locsh(i,2)-wsz:locsh(i,2)+wsz);
    
    [xloc,yloc] = find(Eskw);
    locsk = [xloc yloc];
    
    dist = pdist2([mpt mpt],locsk);
    dmin = min(dist);
    ind = find(dist == dmin);
    
    for j = 1:length(ind)
        val(j,1) = Edw(locsk(ind(j),1),locsk(ind(j),2));
    end
    
    Edn(locsh(i,1),locsh(i,2)) = mean(val);
    
    clear ind val
end
close(h)

Ed = Edn;

% Sizing parameters for analysis window and branch point removal
r = round((2*SFm + 1)*SFr);   wsz = round((2*r + SFm)*SFr);   sz = 2*wsz+1;
se = zeros(r,r);

% Remove all identified branchpoints and interpolate distance across gap
h = waitbar(0,'Interpolating Distances at Branch Points...');
for k = 1:length(bpti);
    
    waitbar(k/length(bpti))
    
    % Pixel coordinates of current branchpoint
    i = bpti(k);    j = bptj(k);
    
    % Analysis windows in both cleaned and distance map images
    Eshw = Esh(i-wsz:i+wsz,j-wsz:j+wsz);
    Edw = Ed(i-wsz:i+wsz,j-wsz:j+wsz);
    
    % Remove distances surrounding branchpoint
    w = floor(r/2);    mpt = wsz+1;
    Edw(mpt-w:mpt+w,mpt-w:mpt+w) = se;
    
    % Trace all boundaries in current window of cleaned image
    tr = bwboundaries(Eshw,'noholes');
    
    % Create binary mask of extracted boundaries
    Ebou = zeros(sz,sz);
    for n = 1:length(tr);
        for m = 1:length(tr{n,1}(:,1));
            Ebou(tr{n,1}(m,1),tr{n,1}(m,2)) = 1;
        end
    end
    
    % Identify branchpoints in mask
    Eb = bwmorph(Ebou,'branchpoints');
    [bri,brj] = find(Eb);
    
    % Correct for branchpoints by subdividing boundary at identified points
    if ~isempty(bri)
        Eshwtmp = Eshw;
        Eshwtmp(bri,brj) = 0;
        
        Eb = bwmorph(Eshwtmp,'branchpoints');
        [bri,brj] = find(Eb);
        
        % Continue to remove branchpoints until separate boundaries
        while ~isempty(bri)
            Eshwtmp(bri,brj) = 0;
            
            Eb = bwmorph(Eshwtmp,'branchpoints');
            [bri,brj] = find(Eb);
        end
        
        tr = bwboundaries(Eshwtmp,'noholes');
        
        clear Eshwtmp bri brj
    end
    
    %     % Plot analysis window and labeled boundary groups
    %     imshow(Edw)
    %     hold on
    %     color = {'b','r','g','y','m','c','b','r','g','y','m','c'};
    %     for n = 1:length(tr);
    %         plot(tr{n,1}(:,2),tr{n,1}(:,1),color{n},'LineWidth',2)
    %     end
    %     pause(1)
    
    
    % Identify pixel coordinates of missing points/interpolation points
    Esub = logical(Eshw.*Edw);
    [mpti,mptj] = find(Eshw-Esub);
    
    % Interpolate surrounding distance values onto missing points
    for n = 1:length(mpti);
        
        % Identify traced boundary associated with current missing point
        if length(tr) > 1;     % Multiple extracted boundaries
            
            % Search for matching pixel coordinates
            for m = 1:length(tr);
                
                ind = find(tr{m,1}(:,1) == mpti(n) & tr{m,1}(:,2) == mptj(n));
                
                if ~isempty(ind);     % Stop searching if match found
                    trf = tr{m,1};
                    break
                end
                
            end
            
        else                   % Single extracted boundary
            trf = tr{1,1};
        end
        
        % If no matching group was identified...
        if length(tr) > 1 && isempty(ind);
            
            clear dtrmin;    dtrmin = NaN(length(tr),1);
            
            % Search for closest extracted boundary
            for m = 1:length(tr);
                
                clear dtr;    dtr = NaN(length(tr{m,1}),1);
                
                % Point-wise distance between missing point and all boundaries
                dtr = pdist2([mpti(n) mptj(n)],tr{m,1});
                
                dtrmin(m,1) = min(dtr);     % Minimum distance for each boundary
                
            end
            
            % Compute minimum distance (and number of instances) from all boundaries
            dtrm = min(dtrmin);     dtrind = find(dtrmin == dtrm);     numind = length(dtrind);
            
            % Compute length of associated mimimum distance boundaries
            if numind > 1     % Multiple equidistant boundaries
                
                clear len;    len = NaN(numind,1);
                
                % Search for longest equidistant boundary
                for m = 1:numind;
                    len(m,1) = length(tr{dtrind(m),1});
                end
                
                [~,maxi] = max(len);     % Maximum length boundary
                
                trf = tr{dtrind(maxi),1};
                
            else              % Single minimum distance boundary
                trf = tr{dtrind,1};
            end
            
            % Add missing point to identified 'best-associated' boundary
            trf = cat(1,trf,[mpti(n) mptj(n)]);
            
        end
        
        % Remove all non-unique entries in final boundary
        trf = unique(trf,'rows','stable');
        
        % Extract ordered distances from distance map for interpolation
        clear dist;    dist = NaN(length(trf(:,1)),1);
        for m = 1:length(trf(:,1));
            dist(m,1) = Edw(trf(m,1),trf(m,2));
        end
        
        % Get indices of non-zero distance values
        xnz = find(dist);     xr = (1:length(dist))';
        
        % Linearly interpolate distances across all entries in extracted distance vector
        if length(xnz)>=2
            
            clear dinter
            
            % Interpolate zero-valued entries
            for m = 1:length(xnz);
                
                if m == 1                   % Nearest-neighbor interpolation at start
                    dinter(1:xnz(1),1) = dist(xnz(1));
                    
                elseif m == length(xnz);     % Nearest-neighbor interpolation at end
                    dinter(xnz(end):length(dist),1) = dist(xnz(end));
                    
                else                        % Linear interpolation otherwise
                    x = [xnz(m-1) xnz(m)];
                    y = [dist(xnz(m-1)) dist(xnz(m))];
                    
                    lin = polyfit(x,y,1);
                    dinter(xnz(m-1):xnz(m),1) = polyval(lin,xnz(m-1):xnz(m));
                end
                
            end
            
            % Get location of missing point in final boundary
            ind = find(trf(:,1) == mpti(n) & trf(:,2) == mptj(n));
            
            % If not member of boundary, look for closest point in boundary
            if isempty(ind);
                
                clear dtr;    dtr = NaN(length(trf),1);
                
                % Point-wise distance between missing point and boundary
                dtr = pdist2([mpti(n) mptj(n)],trf);
                
                % Minimum distance and indices in current boundary
                dtrm = min(dtr);     ind = find(dtr == dtrm);
                
            end
            
            % Update distance map window with interpolated distance
            dnew = mean(dinter(ind,1));
            Edw(mpti(n),mptj(n)) = dnew;
            
            %             % Plot distance vector before and after interpolation
            %             f1 = figure; plot(xr,dist);    hold on
            %             plot(xr,dinter,'k','LineWidth',2);    pause(1);    close(f1)
            
        end
        
    end
    
    %     % Show updated distance map
    %     f1 = figure; imshow(Edwin2);    pause(1);    close(f1);
    
    % Project updated window back onto full-size distance map image
    Ed(i-wsz:i+wsz,j-wsz:j+wsz) = Edw;
    
end
close(h)

% Get pixel coordinates of non-zero distances in distance map image
[xloc,yloc] = find(Ed);

% Extract, sort, and round distance values for image reconstruction
dpos = [xloc yloc Ed(Ed > 0)];
dpos = double(sortrows(dpos,3));
dpos(:,3) = round(dpos(:,3));

% Initialize radius and identify array of unique radius (distance) values in image
r = 1;    rval = unique(dpos(:,3));

% Allocate memory for reconstructing image
[xsz,ysz] = size(Ed);     Erec = zeros(xsz,ysz);

% Reconstruct image using serial object dilations of increasing size
while r < max(rval);
    
    % Create dummy image for dilations
    Enew = zeros(xsz,ysz);
    
    % Generate mask of pixel coordinates with associated distance r
    pcur = find(dpos(:,3) == r);
    
    if ~isempty(pcur)
        
        for j = 1:length(pcur)
            Enew(dpos(pcur(j),1),dpos(pcur(j),2)) = 1;
        end
        
        % Dilate mask by element of size r to reconstruct distance
        se = strel('disk',r);     % Define disk-shaped structuring element
        Enew = imdilate(Enew,se);
        
        % Add current dilation to final image and convert to binary
        Erec = im2bw(Erec + Enew);
        
    end
    
    % Increase distance/radius increment
    r = r+1;    clear pcur
end

Erec = imdilate(Erec,strel('disk',1));

[xsz,ysz] = size(Erec);

% Initialize all white image in HSL colorspace
Efin = cat(3,360.*ones(xsz,ysz),zeros(xsz,ysz),ones(xsz,ysz));

% Make reconstructed pixels black in HSL colorspace
[xloc,yloc] = find(Erec);
for i = 1:length(xloc)
    Efin(xloc(i),yloc(i),3) = 0;
end

end