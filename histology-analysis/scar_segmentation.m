function scar_segmentation(IDX,fileform)

% ===== Define global variables ===========================================
global path groupnm outnm
global stain ext 
global imsave 

% ===== Scar morphological analysis ============================================
% Extract name and open processed image
fname = outnm{IDX,:};
I = imread(strcat(path,groupnm,'/',fname,'/',fname,'.',ext));

% Open results of partitioning analysis
af = importdata(strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_AF-HSL_circ.mat'));    af = cell2mat(af.af);    % Area fractions
TH = importdata(strcat(path,groupnm,'/',fname,'/',fname,'_thickness.mat'));                th = TH.avg;             % Thickness
LP = importdata(strcat(path,groupnm,'/',fname,'/',fname,'_partition_circ.mat'));                                    % Local partitions

% Register thickness and area fraction locations
midx = NaN(size(LP.localpart,1),1);
for II = 1:size(LP.localpart,1)
    
    % Centroid of current partition
    xbar = mean(LP.localpart{II,1}(:,1));    ybar = mean(LP.localpart{II,1}(:,2));
    
    % Find point in thickness array closest to centroid
    D = pdist2([xbar,ybar],TH.bound.mid,'euclidean');    [~,midx(II,1)] = min(D);
end

% Build array for thickness, area fractions, and circumferential position (theta)
th_avg = th(:,2);    th_mid = TH.bound.mid;
th_pos = th(:,1);    af_pos = th_pos(midx,1);

if strcmpi(stain,'bPSR')
    af_col = af(:,4);    af_tis = af(:,5);
elseif strcmpi(stain,'MTC')
    af_col = af(:,5);    af_tis = af(:,4);
end

% Duplicate periodic data to better identify scar boundaries
doubleth = true;
if doubleth
    th_avg = cat(1,th_avg,th_avg);        th_mid = cat(1,th_mid,th_mid);
    th_pos = cat(1,th_pos,th_pos+360);    af_pos = cat(1,af_pos,af_pos+360);
    af_tis = cat(1,af_tis,af_tis);        af_col = cat(1,af_col,af_col);
end

% Register arrays by sorting based on theta position
afsort = cat(2,af_pos,af_tis,af_col);    afsort = sortrows(afsort,1);

% Interpolate area fractions at each theta position from thickness (increases resolution)
th_tis = interp1(afsort(:,1),afsort(:,2),th_pos,'linear','extrap');
th_col = interp1(afsort(:,1),afsort(:,3),th_pos,'linear','extrap');

% Define locations in arrays for analysis (doubleth takes middle 360deg)
if doubleth
    sz = size(th_col,1);    sidx = round(0.25*sz);    eidx = sidx + (round(0.5*sz)-1);
else
    sidx = 1;   eidx = size(th_col,1);  %#ok
end

% Define location of scar based on intersection (or minimum distance)
% between local collagen and tissue area fractions
[~,~,idx,~] = intersections(th_pos(sidx:eidx),th_tis(sidx:eidx),th_pos(sidx:eidx),th_col(sidx:eidx),true);

idx = round(idx);
if length(idx) >= 2;
    idx1 = idx(1);    idx2 = idx(end);
else
    daf = abs(th_tis(sidx:eidx)-th_col(sidx:eidx));
    dafsort = cat(2,daf,linspace(1,length(daf),length(daf))');    dafsort = sortrows(dafsort,1);
    idx1 = dafsort(1,2);    idx2 = dafsort(2,2);
end

% Correct indices of scar location
chk = round(0.1*(idx2-idx1));
if (th_col(sidx+idx1+chk) < th_tis(sidx+idx1+chk));      idx2 = idx2 - (eidx-sidx);                 end
if (sidx+idx2-1) < 0;                   idx2 = idx2 + (eidx-sidx);    idx1 = idx1 + (eidx-sidx);    end
if (idx1 > idx2);                               tidx = idx1;    idx1 = idx2;    idx2 = tidx;        end

% Extract scar morphology ----------------
% - Circumferential position (theta)
intersect.scar.theta = th_pos(idx1+sidx-1:idx2+sidx-1);

% - Local thickness variation
intersect.scar.thick.local = th_avg(idx1+sidx-1:idx2+sidx-1);
intersect.scar.thick.avg   = mean(intersect.scar.thick.local);
intersect.scar.thick.std   = std(intersect.scar.thick.local);

% - Local area fraction
intersect.scar.af.local = cat(2,th_col(idx1+sidx-1:idx2+sidx-1),th_tis(idx1+sidx-1:idx2+sidx-1));
intersect.scar.af.avg   = mean(intersect.scar.af.local);
intersect.scar.af.std   = std(intersect.scar.af.local);

% - Length and percent circumference
intersect.scar.length.degree = intersect.scar.theta(end) - intersect.scar.theta(1);
intersect.scar.length.percent = 100*(length(intersect.scar.theta)/(length(th_pos)/2));

% Extract non-scar morphology ----------------
if (idx1+eidx-1 < eidx) || (idx2+sidx-1 > eidx);
    idxnon = (idx2+sidx-1:idx1+eidx-1)';
else
    idxnon = cat(1,(idx2+sidx-1:eidx)',(sidx:idx1+sidx-1)');
end

% - Circumferential position (theta)
intersect.nonscar.theta = th_pos(idxnon);

% - Local thickness variation
intersect.nonscar.thick.local = th_avg(idxnon);
intersect.nonscar.thick.avg   = mean(intersect.nonscar.thick.local);
intersect.nonscar.thick.std   = std(intersect.nonscar.thick.local);

% - Local area fraction
intersect.nonscar.af.local  = cat(2,th_col(idxnon),th_tis(idxnon));
intersect.nonscar.af.avg    = mean(intersect.nonscar.af.local);
intersect.nonscar.af.std    = std(intersect.nonscar.af.local);

% - Length and percent circumference
intersect.nonscar.length.degree = 360-intersect.scar.length.degree;
intersect.nonscar.length.percent = 100-intersect.scar.length.percent; %#ok

save(strcat(path,groupnm,'/',fname,'/',fname,'_scar-properties.mat'),'intersect')

% ----- Plotting ----------------------------------------------------------
scrsz = get(0,'ScreenSize');

% Thickness/Area Fraction curves
figure('Position',[0.15*scrsz(3) 0.25*scrsz(4) 0.7*scrsz(3) 0.5*scrsz(4)]);
p1 = plot(th_pos,th_tis,'r','LineWidth',1.5);
hold on;
p2 = plot(th_pos,th_col,'b','LineWidth',1.5);
p3 = plot(th_pos,th_avg/max(th_avg),'k','LineWidth',1.5);
p4 = plot(th_pos(idx1+sidx-1),th_col(idx1+sidx-1),'go','MarkerFaceColor','g','LineWidth',1.5,'MarkerSize',8);
plot(th_pos(idx2+sidx-1),th_col(idx2+sidx-1),'go','MarkerFaceColor','g','LineWidth',1.5,'MarkerSize',8);

xlabel('Circumferential Position (deg)');    ylabel('Histology Metric');
legend([p1 p2 p3 p4],'Tissue','Collagen','Thickness','Intersection','Location','NorthWest');   legend boxoff
if doubleth;    xlim([0,720]);    else    xlim([0,360]);    end  %#ok
ylim([0,1.15]);    set(gca,'FontSize',13);

% Save figure
set(gcf,'PaperPositionMode','auto');
print('-dtiff','-r300',strcat(path,groupnm,'/',fname,'/',fname,'_local-properties.tif'));
close(gcf);

% Identified scar locations
figure;    imshow(I);    hold on;
plot(th_mid(:,2),th_mid(:,1),'k','LineWidth',3)
plot(th_mid(idx1+sidx-1:idx2+sidx-1,2),th_mid(idx1+sidx-1:idx2+sidx-1,1),'y','LineWidth',5)
plot(th_mid(idxnon,2),th_mid(idxnon,1),'c','LineWidth',5)
plot(th_mid(idx1+sidx-1,2),th_mid(idx1+sidx-1,1),'ro','MarkerFaceColor','r','MarkerSize',10)
plot(th_mid(idx2+sidx-1,2),th_mid(idx2+sidx-1,1),'ro','MarkerFaceColor','r','MarkerSize',10)

% Save image
set(gcf,'color','w');    F = getframe(gcf);
imwrite(F.cdata,strcat(path,groupnm,'/',fname,'/',fname,'_scar-location.tif'),'tif');


% % Partition centroids for theta offsets
% figure; imshow(I); hold on
% for II = 1:size(LP.localpart,1);
%     for KK = 1:size(LP.localpart,2);
%         plot(LP.localpart{II,KK}(:,2),LP.localpart{II,KK}(:,1),'LineWidth',3,'Color',[0.0,0.0,0.0]);
%     end
% end
% plot(th_mid(:,2),th_mid(:,1),'b','LineWidth',3)
%
% for II = 1:size(midx,1);
%     plot(th_mid(midx(II,1),2),th_mid(midx(II,1),1),'mo','MarkerFaceColor','m')
% end


% First automatic, then manual correction if needed
proc = questdlg('Manually adjust scar location?','Adjust','Yes','No','No');
while strcmpi(proc,'Yes');
    
    % Show image and midline overlay
    close(gcf);    figure;    imshow(I);    hold on;
    plot(th_mid(:,2),th_mid(:,1),'k','LineWidth',3)
    set(gcf,'Name','Select points at scar boundary in CLOCKWISE order...');
    
    % Select points
    [X,Y] = ginput(1);    X = round(X);    Y = round(Y);
    D = pdist2([Y,X],th_mid(sidx:eidx,:),'euclidean');    [~,idx1m] = min(D);
    plot(th_mid(sidx+idx1m-1,2),th_mid(sidx+idx1m-1,1),'mo','MarkerFaceColor','m','MarkerSize',10)
    
    [X,Y] = ginput(1);    X = round(X);    Y = round(Y);
    D = pdist2([Y,X],th_mid(sidx:eidx,:),'euclidean');    [~,idx2m] = min(D);
    plot(th_mid(sidx+idx2m-1,2),th_mid(sidx+idx2m-1,1),'mo','MarkerFaceColor','m','MarkerSize',10)
    pause(0.5);    close(gcf)
    
    % Correct indices of scar location
    if (sum(th_col(sidx+idx1m-1:sidx+idx2m-1) > th_tis(sidx+idx1m-1:sidx+idx2m-1))/length(sidx+idx1m-1:sidx+idx2m-1)) < 0.9;      idx2m = idx2m - (eidx-sidx);                   end
    if (sidx+idx2m-1) < 0;                                                                                       idx2m = idx2m + (eidx-sidx);    idx1m = idx1m + (eidx-sidx);    end
    if (idx1m > idx2m);                                                                                                  tidx = idx1m;    idx1m = idx2m;    idx2m = tidx;        end

    % Extract scar morphology ----------------
    % - Circumferential position (theta)
    manual.scar.theta = th_pos(sidx+idx1m-1:sidx+idx2m-1);
    
    % - Local thickness variation
    manual.scar.thick.local = th_avg(sidx+idx1m-1:sidx+idx2m-1);
    manual.scar.thick.avg = mean(manual.scar.thick.local);
    manual.scar.thick.std = std(manual.scar.thick.local);
    
    % - Local area fraction
    manual.scar.af.local = cat(2,th_tis(sidx+idx1m-1:sidx+idx2m-1),th_col(sidx+idx1m-1:sidx+idx2m-1));
    manual.scar.af.avg   = mean(manual.scar.af.local);
    manual.scar.af.std   = mean(manual.scar.af.local);
    
    % - Length and percent circumference
    manual.scar.length.degree = manual.scar.theta(end) - manual.scar.theta(1);
    manual.scar.length.percent = 100*(length(manual.scar.theta)/(length(th_pos)/2));
    
    % Extract non-scar morphology ----------------
    if (idx1m+eidx-1 < eidx) || (idx2m+sidx-1 > eidx);
        idxnon = (idx2m+sidx-1:idx1m+eidx-1)';
    else
        idxnon = cat(1,(idx2m+sidx-1:eidx)',(sidx:idx1m+sidx-1)');
    end
    
    % - Circumferential position (theta)
    manual.nonscar.theta     = th_pos(idxnon);
    
    % - Local thickness variation
    manual.nonscar.thick.local = th_avg(idxnon);
    manual.nonscar.thick.avg   = mean(manual.nonscar.thick.local);
    manual.nonscar.thick.std   = std(manual.nonscar.thick.local);
    
    % - Local area fraction
    manual.nonscar.af.local  = cat(2,th_tis(idxnon),th_col(idxnon));
    manual.nonscar.af.avg = mean(manual.nonscar.af.local);
    manual.nonscar.af.std = std(manual.nonscar.af.local);
    
    % - Length and percent circumference
    manual.nonscar.length.degree = 360-manual.scar.length.degree;
    manual.nonscar.length.percent = 100-manual.scar.length.percent;
    
    save(strcat(path,groupnm,'/',fname,'/',fname,'_scar-properties.mat'),'intersect','manual')
    
    % ----- Re-Plotting ----------------------------------------------------------
    scrsz = get(0,'ScreenSize');
    
    % Thickness/Area Fraction curves
    figure('Position',[0.15*scrsz(3) 0.25*scrsz(4) 0.7*scrsz(3) 0.5*scrsz(4)]);
    p1 = plot(th_pos,th_tis,'r','LineWidth',1.5);
    hold on;
    p2 = plot(th_pos,th_col,'b','LineWidth',1.5);
    p3 = plot(th_pos,th_avg/max(th_avg),'k','LineWidth',1.5);
    p4 = plot(th_pos(idx1+sidx-1),th_col(idx1+sidx-1),'go','MarkerFaceColor','g','LineWidth',1.5,'MarkerSize',8);
    plot(th_pos(idx2+sidx-1),th_col(idx2+sidx-1),'go','MarkerFaceColor','g','LineWidth',1.5,'MarkerSize',8);
    p5 = plot(th_pos(idx1m+sidx-1),th_col(idx1m+sidx-1),'mo','MarkerFaceColor','m','LineWidth',1.5,'MarkerSize',8);
    plot(th_pos(idx2m+sidx-1),th_col(idx2m+sidx-1),'mo','MarkerFaceColor','m','LineWidth',1.5,'MarkerSize',8);
    
    xlabel('Circumferential Position (deg)');    ylabel('Histology Metric');
    legend([p1 p2 p3 p4 p5],'Myocardium','Collagen','Thickness','Intersection','Manual','Location','NorthWest');   legend boxoff
    if doubleth;    xlim([0,720]);  else    xlim([0,360]);    end  %#ok
    ylim([0,1.2]);    set(gca,'FontSize',13);
    
    % Save figure
    set(gcf,'PaperPositionMode','auto');
    print('-dtiff','-r300',strcat(path,groupnm,'/',fname,'/',fname,'_local-properties.tif'));
    close(gcf);
    
    % Identified scar locations
    figure;    imshow(I);    hold on;
    plot(th_mid(:,2),th_mid(:,1),'k','LineWidth',3)
    plot(th_mid(sidx+idx1m-1:sidx+idx2m-1,2),th_mid(sidx+idx1m-1:sidx+idx2m-1,1),'y','LineWidth',5)
    plot(th_mid(idxnon,2),th_mid(idxnon,1),'c','LineWidth',5)
    plot(th_mid(idx1m+sidx-1,2),th_mid(idx1m+sidx-1,1),'mo','MarkerFaceColor','m','MarkerSize',10)
    plot(th_mid(idx2m+sidx-1,2),th_mid(idx2m+sidx-1,1),'mo','MarkerFaceColor','m','MarkerSize',10)
    
    % Save image
    set(gcf,'color','w');    F = getframe(gcf);
    imwrite(F.cdata,strcat(path,groupnm,'/',fname,'/',fname,'_scar-location.tif'),'tif');
    
    proc = questdlg('Manually adjust scar location?','Adjust','Yes','No','No');
end
close(gcf)

if exist('idx1m','var');    idx1 = idx1m;    end
if exist('idx2m','var');    idx2 = idx2m;    end

% Split scar and nonscar at identified points
% - Find partitions with scar endpoints
chk1 = false(size(LP.localpart,1));    chk2 = false(size(LP.localpart,1));
for II = 1:size(LP.localpart,1);
    chk1(II,1) = inpoly(cat(2,th_mid(idx1+sidx-1,1),th_mid(idx1+sidx-1,2)),LP.localpart{II});
    chk2(II,1) = inpoly(cat(2,th_mid(idx2+sidx-1,1),th_mid(idx2+sidx-1,2)),LP.localpart{II});
end

% - Update indices, if necessary
if idx1 < 0;    idx1 = idx1 + (eidx-sidx+1);    end
if idx2 < 0;    idx2 = idx2 + (eidx-sidx+1);    end

% - Partition index
pidx1 = find(chk1);    pidx2 = find(chk2);

% - Extract points on inner and outer boundaries of partition
in1  = inpoly(TH.bound.inner,LP.localpart{pidx1});         in2  = inpoly(TH.bound.inner,LP.localpart{pidx2});
out1 = inpoly(TH.bound.outer,LP.localpart{pidx1});         out2 = inpoly(TH.bound.outer,LP.localpart{pidx2});
mid1 = inpoly(th_mid(sidx:eidx,:),LP.localpart{pidx1});    mid2 = inpoly(th_mid(sidx:eidx,:),LP.localpart{pidx2});

% - Reorder boundary points
[~,iP1] = max(diff(in1));     [~,iN1] = min(diff(in1));    if iP1 > iN1;   inR1  = cat(2,iP1:length(in1) ,1:iN1)';   else;   inR1  = (iP1:iN1)';   end
[~,iP2] = max(diff(in2));     [~,iN2] = min(diff(in2));    if iP2 > iN2;   inR2  = cat(2,iP2:length(in2) ,1:iN2)';   else;   inR2  = (iP2:iN2)';   end
[~,oP1] = max(diff(out1));    [~,oN1] = min(diff(out1));   if oP1 > oN1;   outR1 = cat(2,oP1:length(out1),1:oN1)';   else;   outR1 = (oP1:oN1)';   end
[~,oP2] = max(diff(out2));    [~,oN2] = min(diff(out2));   if oP2 > oN2;   outR2 = cat(2,oP2:length(out2),1:oN2)';   else;   outR2 = (oP2:oN2)';   end
[~,mP1] = max(diff(mid1));    [~,mN1] = min(diff(mid1));   if mP1 > mN1;   midR1 = cat(2,mP1:length(mid1),1:mN1)';   else;   midR1 = (mP1:mN1)';   end
[~,mP2] = max(diff(mid2));    [~,mN2] = min(diff(mid2));   if mP2 > mN2;   midR2 = cat(2,mP2:length(mid2),1:mN2)';   else;   midR2 = (mP2:mN2)';   end

% - Extract boundaries
ibound1 = TH.bound.inner(inR1,:);     ibound2 = TH.bound.inner(inR2,:);
obound1 = TH.bound.outer(outR1,:);    obound1 = flipud(obound1);
obound2 = TH.bound.outer(outR2,:);    obound2 = flipud(obound2);

% - Order partition boundaries
ind1 = find(midR1 == idx1);    
if isempty(ind1) && all(idx1 > midR1);    ind1 = find(midR1 == (idx1 - (eidx-sidx)));   elseif isempty(ind1) && all(idx1 < midR1);   ind1 = find(midR1 == (idx1 + (eidx-sidx)));    end    
pos1 = ind1/length(midR1);

ind2 = find(midR2 == idx2);    
if isempty(ind2) && all(idx2 > midR2);    ind2 = find(midR2 == (idx2 - (eidx-sidx)));   elseif isempty(ind2) && all(idx2 < midR2);   ind2 = find(midR2 == (idx2 + (eidx-sidx)));    end    
pos2 = ind2/length(midR2);

% - Extract points on inner and outer boundaries
iind1 = round(pos1*length(inR1));     [~,Iidx1] = min(pdist2(ibound1(iind1,:),TH.bound.inner));
iind2 = round(pos2*length(inR2));     [~,Iidx2] = min(pdist2(ibound2(iind2,:),TH.bound.inner));
oind1 = round(pos1*length(outR1));    [~,Oidx1] = min(pdist2(obound1(oind1,:),TH.bound.outer));
oind2 = round(pos2*length(outR2));    [~,Oidx2] = min(pdist2(obound2(oind2,:),TH.bound.outer));

% - Define inner boundaries of scar/nonscar regions
if Iidx1 > Iidx2;
    IR1 = cat(2,Iidx1:length(TH.bound.inner(:,1)),1:Iidx2)';   IR2 = (Iidx2:Iidx1)';
elseif Iidx2 > Iidx1;
    IR1 = (Iidx1:Iidx2)';    IR2 = cat(2,Iidx2:length(TH.bound.inner(:,1)),1:Iidx1)';
end

% - Define outer boundaries of scar/nonscar regions
outer = flipud(TH.bound.outer);
if Oidx1 > Oidx2;
    OR2 = cat(2,Oidx1:length(outer(:,1)),1:Oidx2)';   OR1 = (Oidx2:Oidx1)';
elseif Oidx2 > Oidx1;
    OR2 = (Oidx1:Oidx2)';    OR1 = cat(2,Oidx2:length(outer(:,1)),1:Oidx1)';
end

% - Generate masks based on boundaries
[Nx0,Ny0,~] = size(I);
[X,Y] = meshgrid(linspace(1,Ny0,Ny0),linspace(1,Nx0,Nx0));
Imsk1 = inpoly(cat(2,X(:),Y(:)),cat(2,cat(1,TH.bound.inner(IR1,2),TH.bound.outer(OR1,2)),cat(1,TH.bound.inner(IR1,1),TH.bound.outer(OR1,1))));    Imsk1 = reshape(Imsk1,size(X));
Imsk2 = inpoly(cat(2,X(:),Y(:)),cat(2,cat(1,TH.bound.inner(IR2,2),TH.bound.outer(OR2,2)),cat(1,TH.bound.inner(IR2,1),TH.bound.outer(OR2,1))));    Imsk2 = reshape(Imsk2,size(X));

% Load pseudocolored masks for global area fraction calculation
COL = imread(strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_collagen.tif'));
TOT = imread(strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_total.tif'));
if strcmpi(stain,'MTC');
    OTH = imread(strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_cytoplasm.tif'));
elseif strcmpi(stain,'bPSR');
    OTH = imread(strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_tissue.tif'));
end
COM = imread(strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_combined.tif'));


% Mask images of scar/nonscar regions
% - collagen
COL1s = COL(:,:,1);    COL1n = COL1s;    COL1s(Imsk1 == false) = 255;    COL1n(Imsk2 == false) = 255;
COL2s = COL(:,:,2);    COL2n = COL2s;    COL2s(Imsk1 == false) = 255;    COL2n(Imsk2 == false) = 255;
COL3s = COL(:,:,3);    COL3n = COL3s;    COL3s(Imsk1 == false) = 255;    COL3n(Imsk2 == false) = 255;
COL_scar = cat(3,COL1s,COL2s,COL3s);     COL_non = cat(3,COL1n,COL2n,COL3n);

% - other tissue
OTH1s = OTH(:,:,1);    OTH1n = OTH1s;    OTH1s(Imsk1 == false) = 255;    OTH1n(Imsk2 == false) = 255;
OTH2s = OTH(:,:,2);    OTH2n = OTH2s;    OTH2s(Imsk1 == false) = 255;    OTH2n(Imsk2 == false) = 255;
OTH3s = OTH(:,:,3);    OTH3n = OTH3s;    OTH3s(Imsk1 == false) = 255;    OTH3n(Imsk2 == false) = 255;
OTH_scar = cat(3,OTH1s,OTH2s,OTH3s);     OTH_non = cat(3,OTH1n,OTH2n,OTH3n);

% - total
TOT1s = TOT(:,:,1);    TOT1n = TOT1s;    TOT1s(Imsk1 == false) = 255;    TOT1n(Imsk2 == false) = 255;
TOT2s = TOT(:,:,2);    TOT2n = TOT2s;    TOT2s(Imsk1 == false) = 255;    TOT2n(Imsk2 == false) = 255;
TOT3s = TOT(:,:,3);    TOT3n = TOT3s;    TOT3s(Imsk1 == false) = 255;    TOT3n(Imsk2 == false) = 255;
TOT_scar = cat(3,TOT1s,TOT2s,TOT3s);     TOT_non = cat(3,TOT1n,TOT2n,TOT3n);

% - combined
COM1s = COM(:,:,1);    COM1n = COM1s;    COM1s(Imsk1 == false) = 255;    COM1n(Imsk2 == false) = 255;
COM2s = COM(:,:,2);    COM2n = COM2s;    COM2s(Imsk1 == false) = 255;    COM2n(Imsk2 == false) = 255;
COM3s = COM(:,:,3);    COM3n = COM3s;    COM3s(Imsk1 == false) = 255;    COM3n(Imsk2 == false) = 255;
COM_scar = cat(3,COM1s,COM2s,COM3s);     COM_non = cat(3,COM1n,COM2n,COM3n);

imstack_scar = cat(4,COL_scar,OTH_scar,TOT_scar,COM_scar);
imstack_non  = cat(4, COL_non, OTH_non, TOT_non, COM_non);

if strcmpi(stain,'MTC');
    cnm = {'Collagen'; 'Cytoplasm'; 'Total';'Combined'};
elseif strcmpi(stain,'bPSR');
    cnm = {'Collagen'; 'Tissue'; 'Total';'Combined'};
end

% Compute global area fractions
COL_H = colorspace('RGB->HSL',COL_scar);    pix_col_scar = length(find(COL_H(:,:,3) < 1));
OTH_H = colorspace('RGB->HSL',OTH_scar);    pix_oth_scar = length(find(OTH_H(:,:,3) < 1));
TOT_H = colorspace('RGB->HSL',TOT_scar);    pix_tot_scar = length(find(TOT_H(:,:,3) < 1));

COL_H = colorspace('RGB->HSL',COL_non);     pix_col_non = length(find(COL_H(:,:,3) < 1));
OTH_H = colorspace('RGB->HSL',OTH_non);     pix_oth_non = length(find(OTH_H(:,:,3) < 1));
TOT_H = colorspace('RGB->HSL',TOT_non);     pix_tot_non = length(find(TOT_H(:,:,3) < 1));

af_col_scar = pix_col_scar/pix_tot_scar;    af_col_non = pix_col_non/pix_tot_non;
af_oth_scar = pix_oth_scar/pix_tot_scar;    af_oth_non = pix_oth_non/pix_tot_non;

af_scar = cat(2,pix_col_scar,pix_oth_scar,pix_tot_scar,af_col_scar,af_oth_scar);
af_non = cat(2,pix_col_non,pix_oth_non,pix_tot_non,af_col_non,af_oth_non);

% Make sure scar and non-scar are defined properly
if af_col_non > af_col_scar;
    af_temp = af_scar;              af_scar = af_non;              af_non = af_temp;
    imstack_temp = imstack_scar;    imstack_scar = imstack_non;    imstack_non = imstack_temp;
end
    
% Save images of scar segmentation
if imsave
    for II = 1:length(cnm);
        imwrite(imstack_scar(:,:,:,II),strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_',strcat(lower(cnm{II,1}(1)),cnm{II,1}(2:end)),'-scar.tif'),'tif');
        imwrite(imstack_non(:,:,:,II),strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'_',strcat(lower(cnm{II,1}(1)),cnm{II,1}(2:end)),'-nonscar.tif'),'tif');
    end
end

% Save area fraction data
for II = 1:length(fileform);
    
    if strcmpi(fileform{II},'.xls');
        
        afcon_per = cell(length(cnm)-2,1);
        afcon_pix = cell(length(cnm)-2,1);
        
        for KK = 1:length(cnm)-2
            afcon_per{KK,1} = strcat(cnm{KK},' AF');
            afcon_pix{KK,1} = strcat(cnm{KK},' Pix');
        end
        
        afcon = cat(1,{'Image Number'},afcon_pix,{'Total Pix'},afcon_per);
        cwidth = 20*ones(1,length(afcon));    cwidth(1) = 15;
        fpath = strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'-interstitial-fib.xls');
        
        xlsformat(fpath,afcon,cwidth,af_scar,af_non);
        
    elseif strcmpi(fileform{II},'.mat');
        
        save(strcat(path,groupnm,'/',fname,'/',stain,'_',fname,'-interstitial-fib.mat'),'af_scar','af_non');

    end
end



% _________________________ nested functions ______________________________
% -------------------------------------------------------------------------
    function xlsformat(fpath, header_af, cwidth_af, af_s, af_n)   % Write formatted text file
        
        % Arrange area fraction / HSL data
        afdata_s = num2cell([1 af_s]);    afdata_n = num2cell([1 af_n]);   
            
        % Compile into single array
        afdata_write_s  = [header_af'; afdata_s];      Na_s = size(afdata_write_s,1);
        afdata_write_n  = [header_af'; afdata_n];      Na_n = size(afdata_write_n,1);
        
        % Write data to file
        if exist(fpath,'file') ~= 0
            delete(fpath);    xlswrite(fpath,afdata_write_s,'Sheet1');    xlswrite(fpath,afdata_write_n,'Sheet2');
        else
            xlswrite(fpath,afdata_write_s,'Sheet1');    xlswrite(fpath,afdata_write_n,'Sheet2');  
        end
        
        % Initiate connection with Excel
        hExcel = actxserver('Excel.Application');
        
        % Open workbook and define sheet names
        hWorkbook = hExcel.Workbooks.Open(fpath);
        hWorksheetAF_s = hWorkbook.Sheets.Item('Sheet1');    hWorksheetAF_s.Name  = 'AF Scar';
        hWorksheetAF_n = hWorkbook.Sheets.Item('Sheet2');    hWorksheetAF_n.Name  = 'AF Nonscar';
        
        % Set column widths
        for ii = 1:length(cwidth_af);   hWorksheetAF_s.Columns.Item(ii).columnWidth = cwidth_af(ii);   end
        for ii = 1:length(cwidth_af);   hWorksheetAF_n.Columns.Item(ii).columnWidth = cwidth_af(ii);   end
        
        % Set column justification
        for ii = 1:length(cwidth_af);   hWorksheetAF_s.Columns.Item(ii).HorizontalAlignment = -4108;   end
        for ii = 1:length(cwidth_af);   hWorksheetAF_n.Columns.Item(ii).HorizontalAlignment = -4108;   end
        
        % Set cell formatting based on type of analysis - AF scar
        for ii = 1:length(header_af);
            
            cells = get(hWorksheetAF_s.cells,'item',1,ii);
            hB = cells.borders;    hBL = get(hB,'item',4);    set(hBL,'linestyle',3);
            
            % First column
            if ii == 1;
                for jj = 1:Na_s;
                    cells = get(hWorksheetAF_s.cells,'item',jj,ii);
                    hB = cells.borders;    hBR = get(hB,'item',2);      set(hBR,'linestyle',3);
                end
            end
            
            % Formatting column index
            fidx = floor(1+((((length(header_af)-1)-1)/2)+1));
            
            % Pix-to-AF column and Last column
            chk = (ii == fidx+1);
            
            if chk
                for jj = 1:Na_s;
                    cells = get(hWorksheetAF_s.cells,'item',jj,ii);
                    hB = cells.borders;    hBL = get(hB,'item',1);      set(hBL,'linestyle',3);
                end
            end
            
            % AF formatting
            chk = (ii > fidx);
            
            if chk;
                for jj = 2:Na_s;
                    cells = get(hWorksheetAF_s.cells,'item',jj,ii);    set(cells,'NumberFormat','0.0000');
                end
            end
            
        end
            
        % Set cell formatting based on type of analysis - AF nonscar
        for ii = 1:length(header_af);
            
            cells = get(hWorksheetAF_n.cells,'item',1,ii);
            hB = cells.borders;    hBL = get(hB,'item',4);    set(hBL,'linestyle',3);
            
            % First column
            if ii == 1;
                for jj = 1:Na_n;
                    cells = get(hWorksheetAF_n.cells,'item',jj,ii);
                    hB = cells.borders;    hBR = get(hB,'item',2);      set(hBR,'linestyle',3);
                end
            end
            
            % Formatting column index
            fidx = floor(1+((((length(header_af)-1)-1)/2)+1));
            
            % Pix-to-AF column and Last column
            chk = (ii == fidx+1);
            
            if chk
                for jj = 1:Na_n;
                    cells = get(hWorksheetAF_n.cells,'item',jj,ii);
                    hB = cells.borders;    hBL = get(hB,'item',1);      set(hBL,'linestyle',3);
                end
            end
            
            % AF formatting
            chk = (ii > fidx);
            
            if chk;
                for jj = 2:Na_n;
                    cells = get(hWorksheetAF_n.cells,'item',jj,ii);    set(cells,'NumberFormat','0.0000');
                end
            end
            
        end
        
        % Save and close Excel file
        hWorkbook.Save;    hWorkbook.Close;    hExcel.Quit;
        
    end


    function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
        %INTERSECTIONS Intersections of curves.
        %   Computes the (x,y) locations where two curves intersect.  The curves
        %   can be broken with NaNs or have vertical segments.
        %
        % Example:
        %   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
        %
        % where X1 and Y1 are equal-length vectors of at least two points and
        % represent curve 1.  Similarly, X2 and Y2 represent curve 2.
        % X0 and Y0 are column vectors containing the points at which the two
        % curves intersect.
        %
        % ROBUST (optional) set to 1 or true means to use a slight variation of the
        % algorithm that might return duplicates of some intersection points, and
        % then remove those duplicates.  The default is true, but since the
        % algorithm is slightly slower you can set it to false if you know that
        % your curves don't intersect at any segment boundaries.  Also, the robust
        % version properly handles parallel and overlapping segments.
        %
        % The algorithm can return two additional vectors that indicate which
        % segment pairs contain intersections and where they are:
        %
        %   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
        %
        % For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
        % (how far along this segment the intersection is).  For example, if I(k) =
        % 45.25 then the intersection lies a quarter of the way between the line
        % segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
        % the vector J and the segments in (X2,Y2).
        %
        % You can also get intersections of a curve with itself.  Simply pass in
        % only one curve, i.e.,
        %
        %   [X0,Y0] = intersections(X1,Y1,ROBUST);
        %
        % where, as before, ROBUST is optional.
        % Version: 2.0, 25 May 2017
        % Author:  Douglas M. Schwarz
        % Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
        % Real_email = regexprep(Email,{'=','*'},{'@','.'})
        % Theory of operation:
        %
        % Given two line segments, L1 and L2,
        %
        %   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
        %   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
        %
        % we can write four equations with four unknowns and then solve them.  The
        % four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
        % L1 and L2, t1 is the distance from the starting point of L1 to the
        % intersection relative to the length of L1 and t2 is the distance from the
        % starting point of L2 to the intersection relative to the length of L2.
        %
        % So, the four equations are
        %
        %    (x1(2) - x1(1))*t1 = x0 - x1(1)
        %    (x2(2) - x2(1))*t2 = x0 - x2(1)
        %    (y1(2) - y1(1))*t1 = y0 - y1(1)
        %    (y2(2) - y2(1))*t2 = y0 - y2(1)
        %
        % Rearranging and writing in matrix form,
        %
        %  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
        %        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
        %   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
        %        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
        %
        % Let's call that A*T = B.  We can solve for T with T = A\B.
        %
        % Once we have our solution we just have to look at t1 and t2 to determine
        % whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
        % line segments cross and we can include (x0,y0) in the output.
        %
        % In principle, we have to perform this computation on every pair of line
        % segments in the input data.  This can be quite a large number of pairs so
        % we will reduce it by doing a simple preliminary check to eliminate line
        % segment pairs that could not possibly cross.  The check is to look at the
        % smallest enclosing rectangles (with sides parallel to the axes) for each
        % line segment pair and see if they overlap.  If they do then we have to
        % compute t1 and t2 (via the A\B computation) to see if the line segments
        % cross, but if they don't then the line segments cannot cross.  In a
        % typical application, this technique will eliminate most of the potential
        % line segment pairs.
        % Input checks.
        if verLessThan('matlab','7.13')
            error(nargchk(2,5,nargin)) %#ok<NCHKN>
        else
            narginchk(2,5)
        end
        % Adjustments based on number of arguments.
        switch nargin
            case 2
                robust = true;
                x2 = x1;
                y2 = y1;
                self_intersect = true;
            case 3
                robust = x2;
                x2 = x1;
                y2 = y1;
                self_intersect = true;
            case 4
                robust = true;
                self_intersect = false;
            case 5
                self_intersect = false;
        end
        % x1 and y1 must be vectors with same number of points (at least 2).
        if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
                length(x1) ~= length(y1)
            error('X1 and Y1 must be equal-length vectors of at least 2 points.')
        end
        % x2 and y2 must be vectors with same number of points (at least 2).
        if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
                length(x2) ~= length(y2)
            error('X2 and Y2 must be equal-length vectors of at least 2 points.')
        end
        % Force all inputs to be column vectors.
        x1 = x1(:);
        y1 = y1(:);
        x2 = x2(:);
        y2 = y2(:);
        % Compute number of line segments in each curve and some differences we'll
        % need later.
        n1 = length(x1) - 1;
        n2 = length(x2) - 1;
        xy1 = [x1 y1];
        xy2 = [x2 y2];
        dxy1 = diff(xy1);
        dxy2 = diff(xy2);
        % Determine the combinations of i and j where the rectangle enclosing the
        % i'th line segment of curve 1 overlaps with the rectangle enclosing the
        % j'th line segment of curve 2.
        % Original method that works in old MATLAB versions, but is slower than
        % using binary singleton expansion (explicit or implicit).
        % [i,j] = find( ...
        % 	repmat(mvmin(x1),1,n2) <= repmat(mvmax(x2).',n1,1) & ...
        % 	repmat(mvmax(x1),1,n2) >= repmat(mvmin(x2).',n1,1) & ...
        % 	repmat(mvmin(y1),1,n2) <= repmat(mvmax(y2).',n1,1) & ...
        % 	repmat(mvmax(y1),1,n2) >= repmat(mvmin(y2).',n1,1));
        % Select an algorithm based on MATLAB version and number of line
        % segments in each curve.  We want to avoid forming large matrices for
        % large numbers of line segments.  If the matrices are not too large,
        % choose the best method available for the MATLAB version.
        if n1 > 1000 || n2 > 1000 || verLessThan('matlab','7.4')
            % Determine which curve has the most line segments.
            if n1 >= n2
                % Curve 1 has more segments, loop over segments of curve 2.
                ijc = cell(1,n2);
                min_x1 = mvmin(x1);
                max_x1 = mvmax(x1);
                min_y1 = mvmin(y1);
                max_y1 = mvmax(y1);
                for k = 1:n2
                    k1 = k + 1;
                    ijc{k} = find( ...
                        min_x1 <= max(x2(k),x2(k1)) & max_x1 >= min(x2(k),x2(k1)) & ...
                        min_y1 <= max(y2(k),y2(k1)) & max_y1 >= min(y2(k),y2(k1)));
                    ijc{k}(:,2) = k;
                end
                ij = vertcat(ijc{:});
                i = ij(:,1);
                j = ij(:,2);
            else
                % Curve 2 has more segments, loop over segments of curve 1.
                ijc = cell(1,n1);
                min_x2 = mvmin(x2);
                max_x2 = mvmax(x2);
                min_y2 = mvmin(y2);
                max_y2 = mvmax(y2);
                for k = 1:n1
                    k1 = k + 1;
                    ijc{k}(:,2) = find( ...
                        min_x2 <= max(x1(k),x1(k1)) & max_x2 >= min(x1(k),x1(k1)) & ...
                        min_y2 <= max(y1(k),y1(k1)) & max_y2 >= min(y1(k),y1(k1)));
                    ijc{k}(:,1) = k;
                end
                ij = vertcat(ijc{:});
                i = ij(:,1);
                j = ij(:,2);
            end
            
        elseif verLessThan('matlab','9.1')
            % Use bsxfun.
            [i,j] = find( ...
                bsxfun(@le,mvmin(x1),mvmax(x2).') & ...
                bsxfun(@ge,mvmax(x1),mvmin(x2).') & ...
                bsxfun(@le,mvmin(y1),mvmax(y2).') & ...
                bsxfun(@ge,mvmax(y1),mvmin(y2).'));
            
        else
            % Use implicit expansion.
            [i,j] = find( ...
                mvmin(x1) <= mvmax(x2).' & mvmax(x1) >= mvmin(x2).' & ...
                mvmin(y1) <= mvmax(y2).' & mvmax(y1) >= mvmin(y2).');
            
        end
        % Find segments pairs which have at least one vertex = NaN and remove them.
        % This line is a fast way of finding such segment pairs.  We take
        % advantage of the fact that NaNs propagate through calculations, in
        % particular subtraction (in the calculation of dxy1 and dxy2, which we
        % need anyway) and addition.
        % At the same time we can remove redundant combinations of i and j in the
        % case of finding intersections of a line with itself.
        if self_intersect
            remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
        else
            remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
        end
        i(remove) = [];
        j(remove) = [];
        % Initialize matrices.  We'll put the T's and B's in matrices and use them
        % one column at a time.  AA is a 3-D extension of A where we'll use one
        % plane at a time.
        n = length(i);
        T = zeros(4,n);
        AA = zeros(4,4,n);
        AA([1 2],3,:) = -1;
        AA([3 4],4,:) = -1;
        AA([1 3],1,:) = dxy1(i,:).';
        AA([2 4],2,:) = dxy2(j,:).';
        B = -[x1(i) x2(j) y1(i) y2(j)].';
        % Loop through possibilities.  Trap singularity warning and then use
        % lastwarn to see if that plane of AA is near singular.  Process any such
        % segment pairs to determine if they are colinear (overlap) or merely
        % parallel.  That test consists of checking to see if one of the endpoints
        % of the curve 2 segment lies on the curve 1 segment.  This is done by
        % checking the cross product
        %
        %   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
        %
        % If this is close to zero then the segments overlap.
        % If the robust option is false then we assume no two segment pairs are
        % parallel and just go ahead and do the computation.  If A is ever singular
        % a warning will appear.  This is faster and obviously you should use it
        % only when you know you will never have overlapping or parallel segment
        % pairs.
        if robust
            overlap = false(n,1);
            warning_state = warning('off','MATLAB:singularMatrix');
            % Use try-catch to guarantee original warning state is restored.
            try
                lastwarn('')
                for k = 1:n
                    T(:,k) = AA(:,:,k)\B(:,k);
                    [unused,last_warn] = lastwarn; %#ok<ASGLU>
                    lastwarn('')
                    if strcmp(last_warn,'MATLAB:singularMatrix')
                        % Force in_range(k) to be false.
                        T(1,k) = NaN;
                        % Determine if these segments overlap or are just parallel.
                        overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
                    end
                end
                warning(warning_state)
            catch err
                warning(warning_state)
                rethrow(err)
            end
            % Find where t1 and t2 are between 0 and 1 and return the corresponding
            % x0 and y0 values.
            in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
            % For overlapping segment pairs the algorithm will return an
            % intersection point that is at the center of the overlapping region.
            if any(overlap)
                ia = i(overlap);
                ja = j(overlap);
                % set x0 and y0 to middle of overlapping region.
                T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
                    min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
                T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
                    min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
                selected = in_range | overlap;
            else
                selected = in_range;
            end
            xy0 = T(3:4,selected).';
            
            % Remove duplicate intersection points.
            [xy0,index] = unique(xy0,'rows');
            x0 = xy0(:,1);
            y0 = xy0(:,2);
            
            % Compute how far along each line segment the intersections are.
            if nargout > 2
                sel_index = find(selected);
                sel = sel_index(index);
                iout = i(sel) + T(1,sel).';
                jout = j(sel) + T(2,sel).';
            end
        else % non-robust option
            for k = 1:n
                [L,U] = lu(AA(:,:,k));
                T(:,k) = U\(L\B(:,k));
            end
            
            % Find where t1 and t2 are between 0 and 1 and return the corresponding
            % x0 and y0 values.
            in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
            x0 = T(3,in_range).';
            y0 = T(4,in_range).';
            
            % Compute how far along each line segment the intersections are.
            if nargout > 2
                iout = i(in_range) + T(1,in_range).';
                jout = j(in_range) + T(2,in_range).';
            end
        end
        % Plot the results (useful for debugging).
        % plot(x1,y1,x2,y2,x0,y0,'ok');
    end

    function y = mvmin(x)
        % Faster implementation of movmin(x,k) when k = 1.
        y = min(x(1:end-1),x(2:end));
    end

    function y = mvmax(x)
        % Faster implementation of movmax(x,k) when k = 1.
        y = max(x(1:end-1),x(2:end));
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