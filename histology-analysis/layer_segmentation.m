function [needsingle, allno, wptemp] = layer_segmentation(IDX, numseg, I, chkseg, wp, procopt, format, single, allno)

global SF
global path      ext
global groupnm   outnm   fname   vartype
global mergenm   rednm   stain   tissue   layernm

% Conditions for different procopt values:
% 0 - save non-variational results (group)
% 1 - process image
% 2 - compile/save non-variational results (individual)
% 3 - compile/save variational results (individual)
% 4 - save variational results (group)

% Processing for individual images...'variational' <=> radial/circumferential variation
if procopt == 1 || procopt == 2
    
    % Scaling factors
    SFr = SF(1);    SFm = SF(2);
    
    % Convert RGB image to HSL image
    H = colorspace('RGB->HSL',I);
    
%     % HSL parameters to isolate all NON-BACKGROUND pixels
%     Hlow = 0;      Hup = 360;      % Hue
%     Slow = 0;      Sup = 1;        % Saturation
%     Llow = 0;      Lup = 0.99;     % Lightness
    
    % Extract all foreground pixels    %%% David was here 2022 %%%
    if strcmpi(stain,'IF') || strcmpi(stain,'dPSR')
        % HSL parameters to isolate all NON-BACKGROUND pixels
        Hlow = 0;      Hup = 360;      % Hue
        Slow = 0;      Sup = 1;        % Saturation
        Llow = 0.11;   Lup = 0.99;     % Lightness
        T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','black');
    else
       % HSL parameters to isolate all NON-BACKGROUND pixels
        Hlow = 0;      Hup = 360;      % Hue
        Slow = 0;      Sup = 1;        % Saturation
        Llow = 0;      Lup = 0.99;     % Lightness
         T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
    end
    
    % Compute total number of pixels
    pixtot = length(find(T(:,:,3) < 1));
    
    % Automatic layer segmentation ----------------------------------------
    skip = 0;                     % Don't skip image unless it is necessary
    if strcmpi(chkseg,'no')
        
        % Only allowed for vasculature sections, for now
        if strcmpi(tissue,'Artery / Vasculature');
            
            % Extract mask of relevant regions for segmentation...
            if strcmpi(stain,'VVG') %|| strcmpi(stain,'MOV')
                
                % HSL parameters to isolate BLACK pixels of ELASTIN
                Hlow = 0;   Hup = 360;      % Hue
                Slow = 0;   Sup = 1;        % Saturation
                Llow = 0;   Lup = 0.24;     % Lightness
                
                % Pseudocolor of HSL threshold image
                hcolor = 'k';
                
                % Blob size and circularity thresholds for removal
                bsz = 500;      csz = 0.6;
                
            elseif strcmpi(stain,'MTC');
                
                % HSL parameters to isolate RED pixels of CYTOPLASM
                Hlow = 255;   Hup = 25;       % Hue
                Slow = 0.1;   Sup = 1;        % Saturation
                Llow = 0.0;   Lup = 0.95;     % Lightness
                
                % Pseudocolor of HSL threshold image
                hcolor = 'r';
                
                % Blob size and circularity thresholds for removal
                bsz = 50;      csz = 0.6;
                
            elseif strcmpi(stain,'IF')
                
                % HSL parameters to isolate GREEN pixels of ELASTIN
                Hlow = 60;    Hup = 200;     % Hue
                Slow = 0;     Sup = 1;       % Saturation
                Llow = 0.1;   Lup = 0.9;     % Lightness
                
                % Pseudocolor of HSL threshold image
                hcolor = 'k';
                
                % Blob size and circularity thresholds for removal
                bsz = 150;      csz = 0.6;
            else
                skip = 1;
            end
        else
            skip = 1;
        end
        
        
        % If automatic processing can be performed... ---------------------
        if ~skip
            
            % Color thresholding of HSL image
            H = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],hcolor,'white');
            
            % Blob analysis to identify and remove insignificant pixel groups
            % - Binary morphological operation to smooth mask
            H(:,:,3) = imcomplement(imclose(imcomplement(H(:,:,3)),strel('disk',1)));
            
            % - Generate connectivity map of all pixel groups in image
            cc = bwconncomp(imcomplement(H(:,:,3)),8);
            
            % - Compute morphology/circularity of pixels in each blob and write data from structure to array
            conarea  = regionprops(cc, 'Area');                 conarea  = cell2mat(struct2cell(conarea));
            conperim = regionprops(cc, 'Perimeter');            conperim = cell2mat(struct2cell(conperim));
            concirc  = 4.*pi.*(conarea./(conperim + pi).^2);
            
            % - Pixel coordinate location of blobs in image
            conloc = regionprops(cc,'PixelList');
            
            % - Identify location of small, round blobs to be removed
            blob_loc = find(conarea < round(bsz*SFr*SFm) | concirc > csz);
            
            % - Find each pixel identified to be in a blob
            for II = 1:length(blob_loc);
                
                % Initialize number of pixels in group
                num = length(conloc(blob_loc(II),1).PixelList(:,1));
                
                % Remove blobs from image (i.e., make blob pixels white)
                for JJ = 1:num
                    H(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),1) = 360;       % HSL code
                    H(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),2) = 0;         %   for
                    H(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),3) = 1;         %  WHITE
                end
                
            end
            
            
            % Preprocessing of extracted segmentation mask...
            if strcmpi(stain,'VVG') || strcmpi(stain,'MOV')   % Additional processing of ELASTIN masks
                
                % Additional cleaning of lamellae (needs higher resolution images)
                if SFm*SFr > 1;    H = removeblebs(H,1);    end
                
            elseif strcmpi(stain,'MTC')   % Additional processing of SMOOTH MUSCLE mask
                
                % Extract properties of objects in images/masks
                Tprop = regionprops(im2bw(T(:,:,3)),'Area','Centroid');
                Hprop = regionprops(imcomplement(im2bw(H(:,:,3))),'Area','Centroid','PixelList');
                
                % Find centroid of total pixel mask
                [~,maxind] = max([Tprop.Area]);     cenpt = Tprop(maxind,1).Centroid;
                
                % Store centroids of all other objects
                ind = 1;
                cents = NaN(length(Hprop)-1,2);
                for II = 1:length(Hprop);
                    if II ~= maxind
                        cents(ind,:) = Hprop(II,1).Centroid;    ind = ind + 1;
                    end
                end
                
                % Pairwise distance between individual centroids
                d = pdist2(cenpt,cents,'euclidean');
                
                
                % Compute angle relative to +x-axis
                % -- x- and y- distance from global centroid
                dx = cents(:,1)-cenpt(1);   dy = cenpt(2)-cents(:,2);
                
                % -- Calculate (and adjust) angle
                th = atan2(dy,dx)*180/pi;   th(th < 0) = th(th < 0)+360;
                
                
                % Store and sort centroid by angle, th
                val = [th d' cents];    val = sortrows(val,1);
                th  = val(:,1);           d = val(:,2);
                
                % Initialize parameters to filter centroids by angle and distance
                thstep = 20;    num = 0;    th = sort(th,'ascend');
                
                % Find objects, in an angular sector, who are "far out" from the global centroid
                outind = [];
                for II = 0:thstep:360-thstep;
                    
                    % Theta boundaries of angular secor
                    lb = II;    ub = lb + thstep;
                    
                    % Centroid distances and percentiles
                    pts = d(th >= lb & th <= ub);   pct = prctile(pts,(25:25:75));   D = 1.25*(pct(3) - pct(1));
                    
                    % Identfy and save indices of "far out" centroids
                    indo = find(pts >= pct(2)+D);   if II == 0; outind = cat(1,outind,indo);    else   outind = cat(1,outind,indo+num);    end
                    
                    % Total number of checked objects to updates indices
                    num = num + length(pts);
                    
                end
                
                % Filter identified objects and remove from image, if appropriate
                sz = 5;
                if ~isempty(outind) && (length(outind) > max(outind) + 2*sz)    % must have enough "far out" centroids
                    
                    % Make sure "far out" centroids are far enough out
                    for II = 1:length(outind);
                        
                        % Extract centroids, if there are enough
                        if outind(II) <= sz;
                            pts = val(outind(II):outind(II)+2*sz,3:4);
                        elseif outind(II) > length(val(:,1))-sz;
                            pts = val(outind(II)-2*sz:outind(II),3:4);
                        else
                            pts = val(outind(II)-sz:outind(II)+sz,3:4);
                        end
                        
                        % Compute pairwise distances
                        dist = pdist2(pts(sz+1,:),pts);
                        
                        % Filter based on minimum and mean distance thresholds
                        dme = mean(dist(dist>0));   dmi = min(dist(dist>0));        if dme < 100 && dmi < 50;   outind(II) = NaN;   end
                        
                    end
                    
                    % Keep those that have not been filtered out
                    outind = outind(isfinite(outind));
                    
                    % Remove all remaining centroids from image to improve layer segmentation
                    for II = 1:length(outind)
                        
                        % Extract index of current centroid to remove
                        curcen = val(outind(II),3:4);   ind = find(cents(:,1) == curcen(1) & cents(:,2) == curcen(2));      if ind > maxind;   ind = ind + 1;   end
                        
                        % Remove object by changing lightness in HSL image
                        num = length(Hprop(ind,1).PixelList(:,1));
                        for j = 1:num;
                            H(Hprop(ind,1).PixelList(j,2),Hprop(ind,1).PixelList(j,1),3) = 1;
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % Continue processing image, if applicable ------------------------
        if ~skip
            
            % Extract image dimension and compute 'City-block' distance map
            [ysz,xsz] = size(H(:,:,3));     Hd = bwdist(H(:,:,3),'cityblock');
            
            % Remove non-compact pixels on object borders
            dis = 2;
            for d = dis-1:-1:1;
                for v = 2:xsz-1;
                    for u = 2:ysz-1;
                        
                        % Current pixel in distance map
                        d0 = Hd(u,v);
                        
                        % Only process if at current distance
                        if d0 < d || d0 > d;    continue;    end
                        
                        % Distances in 4-connected neighborhood
                        d1 = Hd(u+1,v);  d2 = Hd(u,v+1);
                        d3 = Hd(u-1,v);  d4 = Hd(u,v-1);
                        
                        % Lower pixel distance if all neighbors are less than current distance
                        if (d1 <= d && d2 <= d && d3 <= d && d4 <= d);    Hd(u,v) = d-1;    end
                        
                    end
                end
            end
            
            % Close preprocessed mask to generate segmentation mask
            Hf = imfill(Hd,'holes');
            
            % Keep track of object area before (arbf) and after (arf) filling
            arbf = length(find(Hd))/numel(Hd);    arf = arbf;
            
            area = bwconncomp(Hf);
            
            maxarea = 1;
            for ii = 1:length(area.PixelIdxList);
                if length(area.PixelIdxList{ii}) > maxarea;
                    maxarea = length(area.PixelIdxList{ii});
                end
            end
            
            areaper = maxarea/(area.ImageSize(1)*area.ImageSize(2));
            
            numarea = area.NumObjects;
            
            % Dilate and fill objects to identify media of vessel cross-section
            r = 1;     dilate = 0;      achk = false;
            while (achk == false) && (numarea > 1) && (areaper > 0.01 && areaper < 0.25) && (r < 100)
                dilate = 1;
                
                arbf = arf;
                
                Hf = imdilate(Hd,strel('disk',r));
                Hf = imfill(Hf,'holes');
                
                area = bwconncomp(Hf);
                numarea = area.NumObjects;
                
                arf = length(find(Hf))/numel(Hf);
                
                if (arf-arbf) > 0.2
                    achk = true;
                end
                
                r = r + 1;
            end
            
            if dilate
                Hf = imerode(Hf,strel('disk',r-1));
            end
            
            zedge = 5;
            r0 = 25*SFm*SFr;
            Hproc1 = zeros(size(Hf));
            Hproc2 = zeros(size(Hf));
            [hysz,hxsz] = size(Hproc1);
            for ii = 1:2;
                
                r = r0*ii;
                
                eval(strcat('Hproc',num2str(ii),' = imdilate(Hf,strel(''disk'',round(r)));'));
                eval(strcat('Hproc',num2str(ii),' = imfill(Hproc',num2str(ii),',''holes'');'));
                eval(strcat('Hproc',num2str(ii),' = imerode(Hproc',num2str(ii),',strel(''disk'',round(r)));'));
                
                objarea = regionprops(im2bw(eval(strcat('Hproc',num2str(ii)))),'Area');
                objarea = cell2mat(struct2cell(objarea));
                
                ind = find(objarea < max(objarea));
                
                objloc = regionprops(im2bw(eval(strcat('Hproc',num2str(ii)))),'PixelList');
                
                if ~isempty(ind);
                    for i = 1:length(ind);
                        num = objarea(ind(i));
                        for j = 1:num;
                            eval(strcat('Hproc',num2str(ii),'(objloc(ind(i),1).PixelList(j,2),objloc(ind(i),1).PixelList(j,1)) = 0;'));
                        end
                    end
                end
                
                eval(strcat('Hproc',num2str(ii),'(1:zedge,:) = 0;'));
                eval(strcat('Hproc',num2str(ii),'(hysz-(zedge-1):hysz,:) = 0;'));
                eval(strcat('Hproc',num2str(ii),'(:,1:zedge) = 0;'));
                eval(strcat('Hproc',num2str(ii),'(:,hxsz-(zedge-1):hxsz) = 0;'));
            end
            
            area = bwconncomp(Hproc1);
            
            if area.NumObjects > 0
                
                if strcmpi(stain,'IF') || strcmpi(stain,'dPSR')
                    areaper = length(area.PixelIdxList{:})/(area.ImageSize(1)*area.ImageSize(2));
                else
                    idx = area.PixelIdxList{:};            areaper = length(idx)/(pixtot);
                end
                
            else
                areaper = 0;
            end
            
        else
            areaper = 0;
        end
        
        if areaper > 0.15
            
            % Pixel coordinates of current contour
            bb = bwboundaries(Hproc1);        bb = bb{1,1};
            sbb = bwboundaries(Hproc2);       sbb = sbb{1,1};
            
            % figure; imshow(I); hold on; plot(bb(:,2),bb(:,1),'LineWidth',2); plot(sbb(:,2),sbb(:,1),'LineWidth',2);
            
            % If no line-indices, assume a x(1) connected with x(2), x(3) with x(4) ...
            seg = [(1:(size(sbb,1)-1))' (2:size(sbb,1))'];
            
            % Calculate tangent vectors
            DT = sbb(seg(:,1),:)-sbb(seg(:,2),:);
            
            % Make influence of tangent vector 1/Distance
            % (Weighted Central Differences. Points which are closer give a more accurate estimate of the normal)
            LL = sqrt(DT(:,1).^2+DT(:,2).^2);
            DT(:,1) = DT(:,1)./max(LL.^2,eps);
            DT(:,2) = DT(:,2)./max(LL.^2,eps);
            
            D1 = zeros(size(sbb)); D1(seg(:,1),:) = DT;
            D2 = zeros(size(sbb)); D2(seg(:,2),:) = DT;
            Ds = D1+D2;
            
            % Normalize the normal
            LL = sqrt(Ds(:,1).^2+Ds(:,2).^2);
            Nsbb(:,1) = -Ds(:,2)./LL;
            Nsbb(:,2) = Ds(:,1)./LL;
            
            for ii = 1:length(sbb(:,1));
                A = sbb(ii,:);  B = sbb(ii,:);
                d0 = pdist2(bb,A,'euclidean','Smallest',1);
                d = d0; idx = 0;
                while d > 0
                    B = round([B(1)+Nsbb(ii,1) B(2)+Nsbb(ii,2)]);
                    [dn,idx] = pdist2(bb,B,'euclidean','Smallest',1);
                    if dn < d0
                        d = dn;
                    else
                        d = 0;
                    end
                end
                D(ii,1) = sqrt((B(2)-A(2))^2+(B(1)-A(1))^2);
                D(ii,2) = idx;
            end
            
            % Find indices of regions to remove
            dthresh = 15*SFm*SFr;
            Didx = find(D(:,1)>dthresh);
            
            midind = round(length(D(:,1))/2);
            Dpad = cat(1,D(midind:end,:),D,D(1:midind,:));
            
            si = NaN(length(Didx),1);
            ei = NaN(length(Didx),1);
            for ii = 1:length(Didx);
                st = Dpad(Didx(ii)+midind,2); off = 1;
                while st > 0;
                    st = Dpad(Didx(ii)+midind-off,2);
                    off = off+1;
                end
                si(ii,1) = Dpad(Didx(ii)+midind-off+2,2)-2;
                
                et = Dpad(Didx(ii)+midind,2); off = 1;
                while et > 0;
                    et = Dpad(Didx(ii)+midind+off,2);
                    off = off+1;
                end
                ei(ii,1) = Dpad(Didx(ii)+midind+off-2,2)+2;
            end
            
            si = unique(si);                          ei = unique(ei);
            si(si <= 0) = 1;                          ei(ei <= 0) = 1;
            
            for ii = 1:length(si);
                bb(si(ii):ei(ii),:) = NaN;
            end
            
            Hb = false(size(Hf));
            for ii = 1:length(bb(:,1));
                if all(isfinite(bb(ii,:)))
                    Hb(bb(ii,1),bb(ii,2)) = true;
                end
            end
            
            for ii = 1:length(si);
                Ic = false(size(Hb));
                
                off = 1;
                if (si(ii)-off) > 1;    A = bb(si(ii)-off,:);   else   A = bb(1,:);    end
                
                while any(isnan(A));
                    off = off + 1;
                    if (si(ii)-off) > 1;    A = bb(si(ii)-off,:);   else   A = bb(end,:);  end
                end
                
                off = 1;
                if (ei(ii)+off) < length(bb);   B = bb(ei(ii)+off,:);   else   B = bb(end,:);    end
                
                while any(isnan(B));
                    off = off + 1;
                    if (ei(ii)+off) < length(bb);   B = bb(ei(ii)+off,:);   else   B = bb(1,:);    end
                end
                
                Ic = linept(Ic,A(1),A(2),B(1),B(2));
                Hb = Hb + Ic;
            end
            
            Hproc  = imfill(Hb,'holes');
            Hproci = logical(imcomplement(Hproc));
            Hproc  = logical(Hproc);
            
            if strcmpi(stain,'IF');
                
                IL1(:,:,1) = Hproc.*I(:,:,1);
                IL1(:,:,2) = Hproc.*I(:,:,2);
                IL1(:,:,3) = Hproc.*I(:,:,3);
                
                IL2(:,:,1) = Hproci.*I(:,:,1);
                IL2(:,:,2) = Hproci.*I(:,:,2);
                IL2(:,:,3) = Hproci.*I(:,:,3);
                
            else
                
                
                I1 = I(:,:,1);  I2 = I(:,:,2);  I3 = I(:,:,3);
                
                IL1_1 = uint8(255.*ones(size(I1))); IL1_2 = IL1_1;  IL1_3 = IL1_1;
                IL2_1 = IL1_1;                      IL2_2 = IL1_1;  IL2_3 = IL1_1;
                
                IL1_1(Hproc == 1) = I1(Hproc == 1);
                IL1_2(Hproc == 1) = I2(Hproc == 1);
                IL1_3(Hproc == 1) = I3(Hproc == 1);
                
                IL1 = cat(3,IL1_1,IL1_2,IL1_3);
                
                IL2_1(Hproci == 1) = I1(Hproci == 1);
                IL2_2(Hproci == 1) = I2(Hproci == 1);
                IL2_3(Hproci == 1) = I3(Hproci == 1);
                
                IL2 = cat(3,IL2_1,IL2_2,IL2_3);
                
            end
            
            needsingle = 0;
            allno = '';
            wptemp = NaN;
            
        else
            
            % Do nothing with current image to allow manual segmentation
            if strcmpi(stain,'IF') || strcmpi(stain,'dPSR')    %%% David was here 2022 %%%
                IL1 = uint8(zeros(size(I)));
            else
                IL1 = uint8(255.*ones(size(I)));
            end
            IL2 = I;
            
            needsingle = 1;
            allno = '';
            wptemp = NaN;
        end
        
    elseif strcmpi(chkseg,'yes')
        
        needsingle = NaN;
        
        % Extract name of image
        fname = outnm{IDX,:};
        
        % Update filename for saving
        if strcmpi(stain,'IF');    mfname = strrep(fname,rednm,mergenm);    else    mfname = fname;    end
        
        if strcmpi(tissue,'Artery / Vasculature') && ~single;
            
            layerlc = cat(2,layernm(numseg),layernm(numseg+1));
            layeruc = cat(2,upper(layernm(numseg)),upper(layernm(numseg+1)));
            
            IL1 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext));
            IL2 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext));
            
        elseif strcmpi(tissue,'Artery / Vasculature') && single;
            
            if numseg == 1 && length(layernm) == 3;
                
                layerlc = cat(2,layernm(numseg+1),layernm(numseg+2));
                layeruc = cat(2,upper(layernm(numseg+1)),upper(layernm(numseg+2)));
                
                IL1 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext));
                IL2 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+2},'.',ext));
                
            elseif numseg == 1 && length(layernm) == 2;
                
                layerlc = cat(2,layernm(numseg),layernm(numseg+1));
                layeruc = cat(2,upper(layernm(numseg)),upper(layernm(numseg+1)));
                
                IL1 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext));
                IL2 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext));
                
            elseif numseg == 2;
                
                layerlc = cat(2,layernm(numseg-1),layernm(numseg));
                layeruc = cat(2,upper(layernm(numseg-1)),upper(layernm(numseg)));
                
                IL1 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg-1},'.',ext));
                IL2 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext));
                
            end
            
        elseif ~strcmpi(tissue,'Artery / Vasculature');
            
            layerlc = cat(2,layernm(end-numseg),layernm(end-numseg+1));
            layeruc = cat(2,upper(layernm(end-numseg)),upper(layernm(end-numseg+1)));
            
            IL1 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{end-numseg},'.',ext));
            IL2 = imread(strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{end-numseg+1},'.',ext));
            
        end
        
        if ~strcmpi(allno,'No to All');
            
            tightsub = @(m,n,p) subtightplot(m, n, p, [0.015 0.015], [0.015 0.015], [0.015 0.015]);
            
            f1 = figure;    tightsub(1,2,1);    imshow(IL1);    title(layeruc{1});  tightsub(1,2,2);    imshow(IL2);    title(layeruc{2});
            set(gcf,'Name',char(strcat(mfname,' -- ',{' '},upper(layerlc{1}(1)),layerlc{1}(2:end),'-',upper(layerlc{2}(1)),layerlc{2}(2:end),' //',{' '},'Previously Processed Images')));   jf = get(gcf,'JavaFrame');   pause(0.1);   set(jf,'Maximized',1);
            proc = questdlg('Adjust segmented images?',char(strcat(upper(layerlc{1}(1)),layerlc{1}(2:end),'-',upper(layerlc{2}(1)),layerlc{2}(2:end),' //',{' '},mfname)),'Yes','No','No to all','No');   pause(0.1);   close all
            
            allno = proc;
            
        else
            proc = allno;
        end
        
        if strcmpi(proc,'Yes');
            
            rmvchk = [1; 0];
            regionuc = {' INSIDE'; ' OUTSIDE'};
            
            imnm = {'IL1','IL2'};
            imth = {'IL2','IL1'};
            
            for II = 1:2;
                
                % Store images to be processed
                Io = eval(imnm{II});        Ith = eval(imth{II});
                
                % Show current image
                f2 = figure;  imshow(Ith)
                jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);
                
                remreg = 'Yes';
                while strcmpi(remreg,'Yes')
                    
                    % Isolate region of interest in thresholded image
                    if strcmpi(tissue,'Artery / Vasculature');
                        set(gcf,'Name',char(strcat('Area', regionuc{II},' of selection will be removed and added back to the',{' '},layeruc{II})))
                    else
                        set(gcf,'Name',char(strcat('Area', regionuc{II},' of selection will be removed and added back to',{' '},layeruc{II})))
                    end
                    
                    % Create mask from traced polygon    %%% OPTIONS %%%
%                     h = impoly(gca);    mask = h.createMask();
                    h = imfreehand(gca);    mask = h.createMask();
                    
                    % Split each image into individual channels
                    Io1 = Io(:,:,1);     Ith1 = Ith(:,:,1);
                    Io2 = Io(:,:,2);     Ith2 = Ith(:,:,2);
                    Io3 = Io(:,:,3);     Ith3 = Ith(:,:,3);
                    
                    % Set background color and transfer pixels between images
                    if strcmpi(stain,'IF') || strcmpi(stain,'dPSR');   %%% David was here 2022 %%%
                        bground = 0;
                    else
                        bground = 255;
                    end
                    
                    Io1(mask == rmvchk(II) & ~all(Ith == bground,3)) = Ith1(mask == rmvchk(II) & ~all(Ith == bground,3));       Ith1(mask == rmvchk(II)) = bground;
                    Io2(mask == rmvchk(II) & ~all(Ith == bground,3)) = Ith2(mask == rmvchk(II) & ~all(Ith == bground,3));       Ith2(mask == rmvchk(II)) = bground;
                    Io3(mask == rmvchk(II) & ~all(Ith == bground,3)) = Ith3(mask == rmvchk(II) & ~all(Ith == bground,3));       Ith3(mask == rmvchk(II)) = bground;
                    
                    % Recompile and store cleaned images
                    Io  = cat(3,Io1,Io2,Io3);       eval(strcat(imnm{II},' = Io;'));
                    Ith = cat(3,Ith1,Ith2,Ith3);    eval(strcat(imth{II},' = Ith;'));
                    
                    % Show cleaned image
                    close(f2);      f2 = figure;    imshow(Ith);
                    jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);
                    
                    % Determine if more processing is needed...
                    if strcmpi(tissue,'Artery / Vasculature');
                        remreg = questdlg(char(strcat('Move more areas back to the',{' '},layerlc{II},'?')),'Continue?','Yes','No','No');
                    else
                        remreg = questdlg(char(strcat('Move more areas back to ',{' '},layerlc{II},'?')),'Continue?','Yes','No','No');
                    end
                    
                end
                
                close all
                
            end
            
        end
        
        % Compute HSL color transformation
        HL1 = colorspace('RGB->HSL',IL1);
        HL2 = colorspace('RGB->HSL',IL2);
        
        %%% David was here 2022 %%%
        if strcmpi(stain,'IF') || strcmpi(stain,'dPSR')
            % HSL parameters to isolate all NON-BLACK pixels
            Hlow = 0;    Hup = 360;      % Hue
            Slow = 0;    Sup = 1;        % Saturation
            Llow = 0.11; Lup = 0.99;     % Lightness
            
            % Isolate all non-white pixels and percentages between the current layers 
            L1 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],HL1,[],'g','black');      pixL1 = length(find(L1(:,:,3) < 1));      pL1 = pixL1/pixtot;
            L2 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],HL2,[],'g','black');      pixL2 = length(find(L2(:,:,3) < 1));      pL2 = pixL2/pixtot;       
        else
            % HSL parameters to isolate all NON-WHITE pixels
            Hlow = 0;    Hup = 360;      % Hue
            Slow = 0;    Sup = 1;        % Saturation
            Llow = 0;    Lup = 0.99;     % Lightness
            
            % Isolate all non-white pixels and percentages between the current layers 
            L1 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],HL1,[],'g','white');      pixL1 = length(find(L1(:,:,3) < 1));      pL1 = pixL1/pixtot;
            L2 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],HL2,[],'g','white');      pixL2 = length(find(L2(:,:,3) < 1));      pL2 = pixL2/pixtot;  
        end
        
        % Store wall percentages
        wptemp = cat(2,pixL1,pixL2,pixtot,pL1,pL2);
        
    end
    
    
    % Extract name of image
    fname = outnm{IDX,:};
    
    % Update filename for saving
    if strcmpi(stain,'IF');    mfname = strrep(fname,rednm,mergenm);    else;    mfname = fname;    end
    
    if strcmpi(chkseg,'no')
        
        if strcmpi(tissue,'Artery / Vasculature') && numseg == 1 && length(layernm) == 3;
            
            if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext,'Quality',100);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+2},'.',ext),ext,'Quality',100);
            else
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+2},'.',ext),ext);
            end
            
        elseif strcmpi(tissue,'Artery / Vasculature') && numseg == 1 && length(layernm) == 2;
            
            if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext,'Quality',100);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext,'Quality',100);
            else
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext);
            end
            
        elseif strcmpi(tissue,'Artery / Vasculature') && numseg == 2;
            
            if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg-1},'.',ext),ext,'Quality',100);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext,'Quality',100);
            else
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg-1},'.',ext),ext);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext);
            end
            
        elseif ~strcmpi(tissue,'Artery / Vasculature')
            
            if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext,'Quality',100);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext,'Quality',100);
            else
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext);
            end
            
        end
        
    elseif strcmpi(chkseg,'yes') && procopt == 1;
        
        if strcmpi(tissue,'Artery / Vasculature') && ~single;
            
            if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext,'Quality',100);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext,'Quality',100);
            else
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext);
            end
            
        elseif strcmpi(tissue,'Artery / Vasculature') && single;
            
            if numseg == 1 && length(layernm) == 3;
                
                if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                    imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext,'Quality',100);
                    imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+2},'.',ext),ext,'Quality',100);
                else
                    imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext);
                    imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+2},'.',ext),ext);
                end
                
            elseif numseg == 1 && length(layernm) == 2;
                
                if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                    imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext,'Quality',100);
                    imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext,'Quality',100);
                else
                    imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext);
                    imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg+1},'.',ext),ext);
                end
                
            elseif numseg == 2;
                
                if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                    imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg-1},'.',ext),ext,'Quality',100);
                    imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext,'Quality',100);
                else
                    imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg-1},'.',ext),ext);
                    imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{numseg},'.',ext),ext);
                end
                
            end
            
        elseif ~strcmpi(tissue,'Artery / Vasculature');
            
            if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{end-numseg},'.',ext),ext,'Quality',100);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{end-numseg+1},'.',ext),ext,'Quality',100);
            else
                imwrite(IL1,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{end-numseg},'.',ext),ext);
                imwrite(IL2,strcat(path,groupnm,'\',fname,'\',mfname,'_',layernm{end-numseg+1},'.',ext),ext);
            end
            
        end
    end
    
    fprintf('...Done! \n')
    
    
    
elseif procopt == 3;
    
    if strcmpi(chkseg,'var');
        
        needsingle = NaN;
        allno = NaN;
        
        % Extract name of image
        fname = outnm{IDX,:};
        
        % Update filename for saving
        if strcmpi(stain,'IF');    mfname = strrep(fname,rednm,mergenm);    else    mfname = fname;    end
        
        % Load layer segmentation images for local wall percentages
        if length(layernm) == 2;
            fpstr1 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{1});    IL1 = imread(strcat(fpstr1,'.',ext));    HL1 = colorspace('RGB->HSL',IL1);    HL1_3 = HL1(:,:,3);
            fpstr2 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{2});    IL2 = imread(strcat(fpstr2,'.',ext));    HL2 = colorspace('RGB->HSL',IL2);    HL2_3 = HL2(:,:,3);
        elseif length(layernm) == 3;
            fpstr1 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{1});    IL1 = imread(strcat(fpstr1,'.',ext));    HL1 = colorspace('RGB->HSL',IL1);    HL1_3 = HL1(:,:,3);
            fpstr2 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{2});    IL2 = imread(strcat(fpstr2,'.',ext));    HL2 = colorspace('RGB->HSL',IL2);    HL2_3 = HL2(:,:,3);
            fpstr3 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{3});    IL3 = imread(strcat(fpstr3,'.',ext));    HL3 = colorspace('RGB->HSL',IL3);    HL3_3 = HL3(:,:,3);
        elseif length(layernm) == 4;
            fpstr1 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{1});    IL1 = imread(strcat(fpstr1,'.',ext));    HL1 = colorspace('RGB->HSL',IL1);    HL1_3 = HL1(:,:,3);
            fpstr2 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{2});    IL2 = imread(strcat(fpstr2,'.',ext));    HL2 = colorspace('RGB->HSL',IL2);    HL2_3 = HL2(:,:,3);
            fpstr3 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{3});    IL3 = imread(strcat(fpstr3,'.',ext));    HL3 = colorspace('RGB->HSL',IL3);    HL3_3 = HL3(:,:,3);
            fpstr4 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{4});    IL4 = imread(strcat(fpstr4,'.',ext));    HL4 = colorspace('RGB->HSL',IL4);    HL4_3 = HL4(:,:,3);
        elseif length(layernm) == 5;
            fpstr1 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{1});    IL1 = imread(strcat(fpstr1,'.',ext));    HL1 = colorspace('RGB->HSL',IL1);    HL1_3 = HL1(:,:,3);
            fpstr2 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{2});    IL2 = imread(strcat(fpstr2,'.',ext));    HL2 = colorspace('RGB->HSL',IL2);    HL2_3 = HL2(:,:,3);
            fpstr3 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{3});    IL3 = imread(strcat(fpstr3,'.',ext));    HL3 = colorspace('RGB->HSL',IL3);    HL3_3 = HL3(:,:,3);
            fpstr4 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{4});    IL4 = imread(strcat(fpstr4,'.',ext));    HL4 = colorspace('RGB->HSL',IL4);    HL4_3 = HL4(:,:,3);
            fpstr5 = strcat(path,groupnm,'\',outnm{IDX,:},'\',outnm{IDX,:},'_',layernm{5});    IL5 = imread(strcat(fpstr5,'.',ext));    HL5 = colorspace('RGB->HSL',IL5);    HL5_3 = HL5(:,:,3);
        end
              
        % Load partition data structure
        localpart = importdata(strcat(path,groupnm,'\',fname,'\',mfname,'_partition_',char(vartype),'.mat'));

        % Convert RGB image to HSL image
        H = colorspace('RGB->HSL',I);    H3 = H(:,:,3);    [xsz,ysz,~] = size(H3);
        
        h = waitbar(0,'Calculating Local Wall Percentages...','Name',outnm{IDX,:});
        
        % Compute local wall percentages
        npartq = size(localpart,1);    npartr = size(localpart,2);    wptemp = cell(npartq,npartr);
        for II = 1:npartq;
            for JJ = 1:npartr;
                
                waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                
                % Generate mask from current partition
                Imsk = false(xsz,ysz);    
                for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                Imsk = imfill(Imsk,'holes');
                
                % Extract all non-background pixels
                pixtot = length(find(H3(Imsk == 1) < 1));
                
                % Compute total number of pixels and area fraction
                wpd = NaN(length(layernm),3);
                if length(layernm) == 2;
                    pixL1 = length(find(HL1_3(Imsk == 1) < 1));    pL1 = pixL1/pixtot;    wpd(1,:) = cat(2,pixL1,pixtot,pL1);
                    pixL2 = length(find(HL2_3(Imsk == 1) < 1));    pL2 = pixL2/pixtot;    wpd(2,:) = cat(2,pixL2,pixtot,pL2);
                elseif length(layernm) == 3;
                    pixL1 = length(find(HL1_3(Imsk == 1) < 1));    pL1 = pixL1/pixtot;    wpd(1,:) = cat(2,pixL1,pixtot,pL1);
                    pixL2 = length(find(HL2_3(Imsk == 1) < 1));    pL2 = pixL2/pixtot;    wpd(2,:) = cat(2,pixL2,pixtot,pL2);
                    pixL3 = length(find(HL3_3(Imsk == 1) < 1));    pL3 = pixL3/pixtot;    wpd(3,:) = cat(2,pixL3,pixtot,pL3);
                elseif length(layernm) == 4;
                    pixL1 = length(find(HL1_3(Imsk == 1) < 1));    pL1 = pixL1/pixtot;    wpd(1,:) = cat(2,pixL1,pixtot,pL1);
                    pixL2 = length(find(HL2_3(Imsk == 1) < 1));    pL2 = pixL2/pixtot;    wpd(2,:) = cat(2,pixL2,pixtot,pL2);
                    pixL3 = length(find(HL3_3(Imsk == 1) < 1));    pL3 = pixL3/pixtot;    wpd(3,:) = cat(2,pixL3,pixtot,pL3);
                    pixL4 = length(find(HL4_3(Imsk == 1) < 1));    pL4 = pixL4/pixtot;    wpd(4,:) = cat(2,pixL4,pixtot,pL4);
                elseif length(layernm) == 5;
                    pixL1 = length(find(HL1_3(Imsk == 1) < 1));    pL1 = pixL1/pixtot;    wpd(1,:) = cat(2,pixL1,pixtot,pL1);
                    pixL2 = length(find(HL2_3(Imsk == 1) < 1));    pL2 = pixL2/pixtot;    wpd(2,:) = cat(2,pixL2,pixtot,pL2);
                    pixL3 = length(find(HL3_3(Imsk == 1) < 1));    pL3 = pixL3/pixtot;    wpd(3,:) = cat(2,pixL3,pixtot,pL3);
                    pixL4 = length(find(HL4_3(Imsk == 1) < 1));    pL4 = pixL4/pixtot;    wpd(4,:) = cat(2,pixL4,pixtot,pL4);
                    pixL5 = length(find(HL5_3(Imsk == 1) < 1));    pL5 = pixL5/pixtot;    wpd(5,:) = cat(2,pixL5,pixtot,pL5);
                end
                                
                % Organize wall percentage arrays
                wptemp{II,JJ} = cat(2,wpd(:,1)',wpd(1,2),wpd(:,3)');
                
            end
        end
        
        close(h)
                                
    end
    
    fprintf('...Done! \n')
    
end

% Saving of wall percentage data
if ~isempty(wp)
    
    needsingle = NaN;
    allno = NaN;
    wptemp = NaN;
    
    % Extract name of image
    fname = outnm{IDX,:};
    
    % Update filename for saving
    if strcmpi(stain,'IF');    mfname = strrep(fname,rednm,mergenm);    else    mfname = fname;    end
    
    % Number of images included in wall percentage input
    im = size(wp,1);
    
    % Build array of current file headings
    const_per = cell(length(layernm),1);
    const_pix = cell(length(layernm),1);
    
    for kk = 1:length(layernm)
        const_per{kk,1} = strcat(upper(layernm{kk}(1)),layernm{kk}(2:end),' Percent');
        const_pix{kk,1} = strcat(upper(layernm{kk}(1)),layernm{kk}(2:end),' Pix');
    end
    
    const = cat(1,const_pix,{'Total Pix'},const_per);
    
    for II = 1:length(format);
        
        if strcmpi(format{II},'.txt');
            
            if procopt == 3 || procopt == 4;
                txt1 = '''%s\t%s\t%s\t';
                txt2 = '''Image Number'',''Circ. Part.'',''Rad. Part''';
                txt3 = '''%2.0f\t%2.0f\t%2.0f\t';
            else
                txt1 = '''%s\t';
                txt2 = '''Image Number''';
                txt3 = '''%2.0f\t';
            end
            
            for kk = 1:length(const);
                
                txt1 = strcat(txt1,'%s\t');
                txt2 = strcat(txt2,',''',const{kk},'''');
                
                if kk <= ceil(length(const)/2);
                    txt3 = strcat(txt3,'%8.0f\t');
                else
                    txt3 = strcat(txt3,'%7.4f\t');
                end
            end
            
            if procopt == 0;
                
                txt1 = strcat(txt1,'%s\t\n''');
                txt2 = strcat(txt2,',''Image Name''');
                txt3 = strcat(txt3,'%s\t\n''');
                
                fid = fopen(strcat(path,groupnm,'\All_',groupnm,'_wall percentage.txt'),'w');
                
                eval(strcat('fprintf(fid,',txt1,',',txt2,');'))
                
                for JJ = 1:im
                    eval(strcat('fprintf(fid,',txt3,',JJ,wp(JJ,:),outnm(JJ,:));'));
                end
                
            elseif procopt == 2;
                
                txt1 = strcat(txt1,'\n''');
                txt3 = strcat(txt3,'\n''');
                
                fid = fopen(strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage.txt'),'w');
                
                eval(strcat('fprintf(fid,',txt1,',',txt2,');'))
                eval(strcat('fprintf(fid,',txt3,',im,wp);'));
                
            elseif procopt == 3;
                
                txt1 = strcat(txt1,'\n''');
                txt3 = strcat(txt3,'\n''');
                
                fid = fopen(strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage_',char(vartype),'.txt'),'w');
                
                eval(strcat('fprintf(fid,',txt1,',',txt2,');'))
                
                npartq = size(wp{im,1},1);    npartr = size(wp{im,1},2);
                
                for JJ = 1:npartq;
                    for KK = 1:npartr;
                        eval(strcat('fprintf(fid,',txt3,',im,JJ,KK,wp{im,1}{JJ,KK});'));
                    end
                end
                
            elseif procopt == 4
                
                txt1 = strcat(txt1,'%s\t\n''');
                txt2 = strcat(txt2,',''Image Name''');
                txt3 = strcat(txt3,'%s\t\n''');
                
                fid = fopen(strcat(path,groupnm,'\All_',groupnm,'_wall percentage_',char(vartype),'.txt'),'w');

                eval(strcat('fprintf(fid,',txt1,',',txt2,');'))
                
                for JJ = 1:im
                    
                    npartq = size(wp{JJ,1},1);    npartr = size(wp{JJ,1},2);
                    
                    for KK = 1:npartq;
                        for PP = 1:npartr;
                            eval(strcat('fprintf(fid,',txt3,',JJ,KK,PP,wp{JJ,1}{KK,PP},outnm(JJ,:));'));
                        end
                    end
                    
                end
                
            end
            
            fclose(fid);
            
        elseif strcmpi(format{II},'.xls');
            
            if procopt == 0;
                cwidth = 16*ones(1,length(const)+2);    cwidth(ceil(length(const)/2)+1:end-1) = 20;    cwidth(end) = 25;
            elseif procopt == 2;
                cwidth = 16*ones(1,length(const)+1);    cwidth(ceil(length(const)/2)+1:end)   = 20;
            elseif procopt == 3;
                cwidth = 16*ones(1,length(const)+3);    cwidth(ceil(length(const)/2)+4:end)   = 20;
            elseif procopt == 4;
                cwidth = 16*ones(1,length(const)+3);    cwidth(ceil(length(const)/2)+4:end)   = 20;    cwidth(end) = 25;
            end
            
            xlsformat(im, const, cwidth, wp, mfname, procopt);
            
        elseif strcmpi(format{II},'.mat');
            
            if procopt == 0;
                save(strcat(path,groupnm,'\All_',groupnm,'_wall percentage.mat'),'wp');
            elseif procopt == 2;
                save(strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage.mat'),'wp');
            elseif procopt == 3;
                save(strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage_',char(vartype),'.mat'),'wp');
            elseif procopt == 4;
                save(strcat(path,groupnm,'\All_',groupnm,'_wall percentage_',char(vartype),'.mat'),'wp');
            end
            
        end
        
    end
    
end



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
        
        for x = max(1,X1):sign(X2 - X1):max(1,X2)
            y = round(f(x, X1, Y1, X2, Y2));
            if y > 0;    result(x,y) = 1;    end
        end
        
        for y = max(1,Y1):sign(Y2 - Y1):max(1,Y2)
            x = round(f2(y, X1, Y1, X2, Y2));
            if x > 0;    result(x,y) = 1;     end
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


    function h = subtightplot(m, n, p, gap, marg_h, marg_w, varargin)
        % Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the gap between
        % neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the gap between
        % subplots can reach 40% of figure area, which is pretty lavish.
        %
        % Input arguments (defaults exist):
        %   gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
        %            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
        %            relatively large axis.
        %   marg_h  margins in height in normalized units (0...1)
        %            or [lower uppper] for different lower and upper margins
        %   marg_w  margins in width in normalized units (0...1)
        %            or [left right] for different left and right margins
        %
        % Output arguments: same as subplot- none, or axes handle according to function call.
        %
        % Issues & Comments: Note that if additional elements are used in order to be passed to subplot, gap parameter must
        %       be defined. For default gap value use empty element- [].
        %
        % Usage example: h=subtightplot((2,3,1:2,[0.5,0.2])
        
        if (nargin<4) || isempty(gap),          gap = 0.01;      end
        if (nargin<5) || isempty(marg_h),    marg_h = 0.05;      end
        if (nargin<5) || isempty(marg_w),    marg_w = marg_h;    end
        if isscalar(gap),                    gap(2) = gap;       end
        if isscalar(marg_h),              marg_h(2) = marg_h;    end
        if isscalar(marg_w),              marg_w(2) = marg_w;    end
        
        gap_vert   = gap(1);       gap_horz   = gap(2);
        marg_lower = marg_h(1);    marg_upper = marg_h(2);
        marg_left  = marg_w(1);    marg_right = marg_w(2);
        
        %note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
        [subplot_col,subplot_row] = ind2sub([n,m],p);
        
        % note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
        subplot_cols = 1 + max(subplot_col) - min(subplot_col); % number of column elements in merged subplot
        subplot_rows = 1 + max(subplot_row) - min(subplot_row); % number of row elements in merged subplot
        
        % single subplot dimensions:
        height = (1 - (marg_lower + marg_upper) - (m-1)*gap_vert)/m;
        width  = (1 - (marg_left  + marg_right) - (n-1)*gap_horz)/n;
        
        % merged subplot dimensions:
        merged_height = subplot_rows*(height + gap_vert) - gap_vert;
        merged_width  = subplot_cols*(width  + gap_horz) - gap_horz;
        
        % merged subplot position:
        merged_bottom = (m - max(subplot_row))*(height + gap_vert) + marg_lower;
        merged_left   = (min(subplot_col) - 1)*(width + gap_horz) + marg_left;
        pos_vec       = [merged_left merged_bottom merged_width merged_height];
        
        h = subplot('Position',pos_vec,varargin{:});
        
        if (nargout < 1),  clear h;  end
        
    end


    function xlsformat(im, const, cwidth, wp, mfname, procopt)
        
        if ~isempty(wp)
            
            % Define column headers
            labels = ceil(length(const(:,1))/2);    
            header = cell(1,2*labels);              header(1) = {'Image Number'};    header(2:end) = const;
            
            % Compile data and set file paths based on type of analysis 
            if procopt == 0;
                
                header(1,end+1) = {'Image Name'};    fpath = strcat(path,groupnm,'\All_',groupnm,'_wall percentage.xls');
                
                for NN = 1:im;    
                    data(NN,:) = cat(2,num2cell([NN wp(NN,:)]),outnm(NN,:));   
                end
                                
            elseif procopt == 2;
                
                data = num2cell([1 wp(1,:)]);    fpath = strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage.xls');
                
            elseif procopt == 3;
                
                header = cat(2,header(1),{'Circ. Part.'},{'Rad. Part'},header(2:end));
                
                data = [];    npq = size(wp{im,1},1);    npr = size(wp{im,1},2);
                
                for JJJ = 1:npq;
                    for KKK = 1:npr
                        data = cat(1,data,num2cell([1 JJJ KKK wp{im,1}{JJJ,KKK}]));
                    end
                end
                
                fpath = strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage_',char(vartype),'.xls');
                
            elseif procopt == 4;
                
                header = cat(2,header(1),{'Circ. Part.'},{'Rad. Part'},header(2:end));    header(1,end+1) = {'Image Name'};
                
                data = [];    npq = NaN(im,1);    npr = NaN(im,1);
                for NN = 1:im;
                    
                    npq(NN,1) = size(wp{NN,1},1);    npr(NN,1) = size(wp{NN,1},2);
                    
                    for JJJ = 1:npq(NN,1);
                        for KKK = 1:npr(NN,1)
                            data = cat(1,data,cat(2,num2cell([NN JJJ KKK wp{NN,1}{JJJ,KKK}]),outnm(NN,:)));
                        end
                    end
                    
                end
                
                fpath = strcat(path,groupnm,'\All_',groupnm,'_wall percentage_',char(vartype),'.xls');
                
            end
            
            % Compile into single array
            data_write = [header; data];
                        
            if procopt == 4;
                pszo = 0;
                for III = 1:im;
                    pszn = npq(III,1)*npr(III,1);
                    data_write(3+pszo:pszo+pszn+1,end) = cell(pszn-1,1);
                    pszo = pszo + pszn;
                end
            end
            
            % Write data to file
            if exist(fpath,'file') ~= 0
                delete(fpath);    xlswrite(fpath,data_write,'Sheet1');
            else
                xlswrite(fpath,data_write,'Sheet1');
            end
            
            % Initiate connection with Excel
            hExcel = actxserver('Excel.Application');
            
            % Open workbook and define sheet names
            if procopt == 0;
                hWorkbook = hExcel.Workbooks.Open(strcat(path,groupnm,'\All_',groupnm,'_wall percentage.xls'));
            elseif procopt == 2;
                hWorkbook = hExcel.Workbooks.Open(strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage.xls'));
            elseif procopt == 3;
                hWorkbook = hExcel.Workbooks.Open(strcat(path,groupnm,'\',fname,'\',mfname,'_wall percentage_',char(vartype),'.xls'));
            elseif procopt == 4;
                hWorkbook = hExcel.Workbooks.Open(strcat(path,groupnm,'\All_',groupnm,'_wall percentage_',char(vartype),'.xls'));
            end
            
            hWorksheet = hWorkbook.Sheets.Item('Sheet1');    hWorksheet.Name = 'Wall Percentage';
            
            % Set column widths and justification
            for NN = 1:length(header);     hWorksheet.Columns.Item(NN).columnWidth = cwidth(NN);      end
            for NN = 1:length(header);    hWorksheet.Columns.Item(NN).HorizontalAlignment = -4108;    end
            
            % Set cell formatting based on type of analysis
            for NN = 1:length(header);

                cells = get(hWorksheet.cells,'item',1,NN);
                
                hB = cells.borders;    hBL = get(hB,'item',4);    set(hBL,'linestyle',3)
                
                if procopt == 3 || procopt == 4;
                    
                    if NN == 1 || NN == 3;
                        for MM = 1:npq*npr+1;
                            cells = get(hWorksheet.cells,'item',MM,NN);
                            hB = cells.borders;    hBR = get(hB,'item',2);    set(hBR,'linestyle',3)
                        end
                    end

                else
                    
                    if NN == 1;
                        for MM = 1:im+1;
                            cells = get(hWorksheet.cells,'item',MM,NN);
                            hB = cells.borders;    hBR = get(hB,'item',2);    set(hBR,'linestyle',3)
                        end
                    end
                    
                end
                
                if procopt == 0;
                    colform = length(header)-labels;      cind = length(header)-labels+1;
                elseif procopt == 2;
                    colform = length(header)-labels+1;    cind = length(header)-labels+2;
                elseif procopt == 3;
                    colform = length(header)-labels+1;    cind = length(header)-labels+2;
                elseif procopt == 4;
                    colform = length(header)-labels;      cind = length(header)-labels+1;
                end
                
                if procopt == 0 || procopt == 2;
                    
                    if NN > colform;
                        
                        if NN == cind || (procopt == 0 && NN == length(header));
                            for MM = 1:im+1;
                                cells = get(hWorksheet.cells,'item',MM,NN);
                                hB = cells.borders;    hBL = get(hB,'item',1);    set(hBL,'linestyle',3)
                            end
                        end
                        
                        for MM = 1:im;
                            cells = get(hWorksheet.cells,'item',MM+1,NN);    set(cells,'NumberFormat','0.0000')
                        end
                    end
                    
                elseif procopt == 3
                    
                    if NN > colform;
                        
                        if NN == cind
                            for MM = 1:npq*npr+1;
                                cells = get(hWorksheet.cells,'item',MM,NN);
                                hB = cells.borders;    hBL = get(hB,'item',1);    set(hBL,'linestyle',3)
                            end
                        end
                        
                        for MM = 1:npq*npr;
                            cells = get(hWorksheet.cells,'item',MM+1,NN);    set(cells,'NumberFormat','0.0000')
                        end
                    end
                    
                elseif procopt == 4
                    
                    if NN > colform;
                        
                        if NN == cind || NN == length(header)
                            
                            pszo = 0;
                            for III = 1:im;
                                pszn = npq(III,1)*npr(III,1);
                                for MM = 1:pszn+1;
                                    cells = get(hWorksheet.cells,'item',MM,NN);
                                    hB = cells.borders;    hBL = get(hB,'item',1);    set(hBL,'linestyle',3)
                                end
                                pszo = pszo + pszn;
                            end
                            
                        end
                        
                        for MM = 1:npq*npr;
                            cells = get(hWorksheet.cells,'item',MM+1,NN);    set(cells,'NumberFormat','0.0000')
                        end
                        
                    end
                end
                                
            end
            
            % Save and close Excel file
            hWorkbook.Save;    hWorkbook.Close;    hExcel.Quit;
            
        end
        
    end



end