function [imstack, afrac, hslavg] = color_segmentation(IDX,LDX,localpart,fileform,indiv)

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ADD: CASES FOR DIFFERENT TISSUES
%      IMMUNOFLUORESCENCE ANALYSIS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% ===== Define global variables ===========================================
global path groupnm outnm tissue
global stain ext SF layer layernm rembleb
global imsave vartype enhance

% ===== Area fraction analysis ============================================
% Extract name of image
fname = outnm{IDX,:};

% For single images...
if indiv
    
    % Load current image for processing
    if layer && ~isempty(LDX)
        I = imread(strcat(path,groupnm,'\',fname,'\',fname,'_',layernm{LDX,:},'.',ext));
    else
        I = imread(strcat(path,groupnm,'\',fname,'\',fname,'.',ext));
    end
    
    % Convert RGB image to binary and compute area/location of all pixel groups
    if strcmpi(stain,'PSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
        Ibw = im2bw(I,0.01);    bwprop = regionprops(Ibw,'Area','PixelList');
    else
        Ibw = im2bw(I,0.99);    bwprop = regionprops(imcomplement(Ibw),'Area','PixelList');
    end
    
    % Initialize and write pixel areas from structure to array
    bwarea = zeros(length(bwprop(:,1)),1);
    for II = 1:length(bwprop(:,1));    bwarea(II,1) = bwprop(II,1).Area;    end
    
    % Remove small (below threshold) groups (i.e., make pixels white/black) to clean image after global thresholding
    SFr = SF(1);    SFm = SF(2);
    for II = 1:length(bwarea);
        if bwarea(II,1) < round(12*SFr*SFm);
            
            num = length(bwprop(II,1).PixelList(:,1));
            
            if strcmpi(stain,'PSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
                for JJ = 1:num;    I(bwprop(II,1).PixelList(JJ,2),bwprop(II,1).PixelList(JJ,1),:) = 0;      end
            else
                for JJ = 1:num;    I(bwprop(II,1).PixelList(JJ,2),bwprop(II,1).PixelList(JJ,1),:) = 255;    end
            end
            
        end
    end
    
    % ===== Area fraction calculations for current image/stain ============
    % ---------------------------------------------------------------------
    
    % Convert RGB image to HSL image
    if strcmpi(stain,'IHC');    I = 255 - I;    end
    
    if strcmpi(enhance,'Yes');
        Io = imread(strcat(path,groupnm,'\',fname,'\',fname,'.',ext));    
        H = colorspace('RGB->HSL',imadjust(I,stretchlim(Io)));
    elseif  strcmpi(enhance,'No');
        H = colorspace('RGB->HSL',I);
    end
    
    % ===== GLOBAL color extraction =======================================
    if isempty(localpart)                         % no partitions specified
        
        % ----- Birefringent Polymer (POL) --------------------------------
        if strcmpi(stain,'POL');
            
            % HSL parameters to isolate NON-BLACK pixels of POLYMER
            Hlow = 0;      Hup = 360;      % Hue
            Slow = 0;      Sup = 1;        % Saturation
            Llow = 0.01;   Lup = 1;        % Lightness
            
            P = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'b','black');
            
            [Havg_pol, Savg_pol, Lavg_pol] = HSLavg(P,H);
            
            % Compute number of identified constituent pixels
            pix_pol = length(find(P(:,:,3) < 1));      % Polymer
            
            % Convert HSL images back to RGB for plotting/saving images
            Pall = colorspace('HSL->RGB',P);
            
            % Plot routines for output of each group
            imstack = Pall;                                   % Images
            afrac   = pix_pol;                                % Area fractions
            hslavg  = cat(1,Havg_pol,Savg_pol,Lavg_pol);      % HSL values
            
            
            % ----- Brightfield Immunohistochemistry (IHC) --------------------
        elseif strcmpi(stain,'IHC');
            
            % HSL parameters to isolate LIGHT BLUE pixels of inverted IHC
            Hlow = 160;    Hup = 210;      % Hue
            Slow = 0.1;    Sup = 1;        % Saturation
            
            Llow = 0.47;    Lup = 1.0;     % Lightness
%             Llow = 0.27;    Lup = 1.0;     % Lightness
%             Llow = 0.01;    Lup = 1.0;     % Lightness
            
            C = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'k','black');
            
            [Havg_ihc, Savg_ihc, Lavg_ihc] = HSLavg(C,H);
            
            % HSL parameters to isolate all NON-BLACK pixels
            Hlow = 0;    Hup = 360;      % Hue
            Slow = 0;    Sup = 1;        % Saturation
            Llow = 0.01; Lup = 1.0;      % Lightness
%             Llow = 0;   Lup = 1.0;      % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','black');
            
            % Compute number of identified constituent pixels
            pix_ihc = length(find(C(:,:,3) < 1));      % IHC antibody
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_ihc = pix_ihc./pix_tot;      % IHC antibody
            
            % Convert HSL images back to RGB for plotting/saving images
            Call = colorspace('HSL->RGB',C);
            Tall = colorspace('HSL->RGB',T);
            
            % Plot/output routines for each constituent
            imstack = cat(4,Call,Tall);                       % Images
            afrac   = cat(1,pix_ihc,pix_tot,af_ihc)';         % Area fractions
            hslavg  = cat(1,Havg_ihc,Savg_ihc,Lavg_ihc);      % HSL values
            
            
            % ----- Hematoxylin & Eosin (H&E) -----------------------------
        elseif strcmpi(stain,'H&E');
            
            % HSL parameters to isolate FG pixels of NUCLEI
            Hlow = 0;      Hup = 42;      % Hue
            Slow = 0.05;   Sup = 1.0;     % Saturation
            Llow = 0.0;    Lup = 0.78;    % Lightness
            
            F = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'b','white');
            
            % HSL parameters to isolate FG pixels of NUCLEI
            Hlow = 212;    Hup = 42;      % Hue
            Slow = 0.05;   Sup = 1.0;     % Saturation
            Llow = 0.0;    Lup = 0.78;    % Lightness
            
            B = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,F,'o','white');    %#ok
            
            fprintf('H&E ANALYSIS IS NOT YET FUNCTIONAL!!!!!\n');
            imstack = NaN;    afrac = NaN;    hslavg = NaN;
            return
            
            % --------- Verhoeff Van Gieson (VVG) -------------------------
        elseif strcmpi(stain,'VVG');
            
            iso = 2;    % Number of constituents to isolate (Elastin and 'Tissue')
            for con = 1:iso;
                if con == 1
                    
                    % HSL parameters to isolate BLACK pixels of ELASTIN
                    Hlow = 0;   Hup = 360;      % Hue
                    Slow = 0;   Sup = 1;        % Saturation
                    
                    if strcmpi(enhance,'Yes')
                        Llow = 0;   Lup = 0.16;     % Lightness
                    else
                        Llow = 0;   Lup = 0.24;     % Lightness
                    end
                    
                    E = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'k','white');
                    
                    [Havg_eln, Savg_eln, Lavg_eln] = HSLavg(E,H);
                    
                    if strcmpi(tissue,'Artery / Vasculature');
                        
                        % Blob analysis to identify and remove cell nuclei and other small pixel groups
                        %     Note: Cell nuclei appear in VVG stains and look like Elastin when stained
                        
                        E(:,:,3) = imcomplement(imclose(imcomplement(E(:,:,3)),strel('disk',1)));
                        cc = bwconncomp(imcomplement(E(:,:,3)),8);   % Mask of all pixel groups in image
                        
                        % Compute and store area, perimeter, and coordinate location of blobs in image
                        conarea  = regionprops(cc,'Area');          conarea = cell2mat(struct2cell(conarea));
                        conperim = regionprops(cc,'Perimeter');    conperim = cell2mat(struct2cell(conperim));
                        conloc   = regionprops(cc,'PixelList');     concirc = 4.*pi.*(conarea./(conperim + pi).^2);
                        
                        % Find small, round blobs as these are not Elastin
                        if ~isempty(LDX) && ~strcmpi(layernm{LDX,:},'adventitia')
                            blob_loc = find(conarea < round(50*SFr*SFm) | concirc > 0.6);
                        else
                            blob_loc = find(concirc > 0.6);
                        end
                        
                        % Find each pixel identified to be in a blob
                        for II = 1:length(blob_loc);
                            
                            num = length(conloc(blob_loc(II),1).PixelList(:,1));
                            
                            % Remove blobs from image (i.e., make blob pixels white)
                            for JJ = 1:num
                                E(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),1) = 360;
                                E(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),2) = 0;
                                E(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),3) = 1;
                            end
                            
                        end
                        
                        % High-resolution bleb removal
                        if ~isempty(LDX) && ~strcmpi(layernm{LDX,:},'adventitia') && SFm*SFr > 1
                            
                            if strcmpi(rembleb,'Yes') || strcmpi(rembleb,'No')
                                qoption = {'Yes','No','Yes to All','No to All'};
                                qans = genquestdlg(qoption, 'Blebs',qoption,'Perform Bleb Removal?');
                                rembleb = qoption{qans};
                            end
                            
                            if strcmpi(rembleb,'Yes') || strcmpi(rembleb,'Yes to All');    E = removeblebs(E,0);    end
                            
                        end
                    end
                    
                elseif con == 2
                    
                    % HSL parameters to isolate PINK/RED pixels of 'TISSUE'
                    Hlow = 300;    Hup = 17;       % Hue
                    Slow = 0.1;    Sup = 1;        % Saturation
                    Llow = 0.1;    Lup = 0.93;     % Lightness
                    
                    O = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,E,'r','white');
                    
                    [Havg_tis, Savg_tis, Lavg_tis] = HSLavg(O,H);
                    
                end
                
            end
            
            % HSL parameters to isolate all NON-WHITE pixels
            Hlow = 0;   Hup = 360;      % Hue
            Slow = 0;   Sup = 1;        % Saturation
            Llow = 0;   Lup = 0.99;     % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
            
            % Compute number of identified constituent pixels
            pix_eln = length(find(E(:,:,3) < 1));      % Elastin
            pix_tis = length(find(O(:,:,3) < 1));      % Other Tissue
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_eln = pix_eln./pix_tot;      % Elastin
            af_tis = pix_tis./pix_tot;      % Other Tissue
            
            % Convert HSL images back to RGB for plotting/saving images
            Eall = colorspace('HSL->RGB',E);
            Oall = colorspace('HSL->RGB',O);
            Tall = colorspace('HSL->RGB',T);
            
            % Plot/output routines for each constituent
            imstack = cat(4,Eall,Oall,Tall);                                                                  % Images
            afrac   = cat(1,pix_eln,pix_tis,pix_tot,af_eln,af_tis)';                                          % Area fractions
            hslavg  = cat(1,cat(2,Havg_eln,Havg_tis),cat(2,Savg_eln,Savg_tis),cat(2,Lavg_eln,Lavg_tis));      % HSL values
            
            
            % --------- Oil Red O (ORO) -----------------------------------
        elseif strcmpi(stain,'ORO') 
                                           
            % HSL parameters to isolate RED pixels of LIPID
            Hlow = 320;   Hup = 40;       % Hue
            Slow = 0.4;   Sup = 1;        % Saturation
            Llow = 0.3;   Lup = 0.75;     % Lightness
            
            L = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'r','white');
            
            [Havg_lip, Savg_lip, Lavg_lip] = HSLavg(L,H);
            
            % HSL parameters to isolate all NON-WHITE pixels
            Hlow = 0;   Hup = 360;      % Hue
            Slow = 0;   Sup = 1;        % Saturation
            Llow = 0;   Lup = 0.99;     % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
            
            % Compute number of identified constituent pixels
            pix_lip = length(find(L(:,:,3) < 1));      % Lipid
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_lip = pix_lip./pix_tot;      % Lipid
            
            % Convert HSL images back to RGB for plotting/saving images
            Lall = colorspace('HSL->RGB',L);
            Tall = colorspace('HSL->RGB',T);
            
            % Plot/output routines for each constituent
            imstack = cat(4,Lall,Tall);                       % Images
            afrac   = cat(1,pix_lip,pix_tot,af_lip)';         % Area fractions
            hslavg  = cat(1,Havg_lip,Savg_lip,Lavg_lip);      % HSL values
            
            
            % ------ Calcification: Von Kossa (VnK) & Alizarin (AzR) ------
        elseif strcmpi(stain,'VnK') || strcmpi(stain,'AzR');
            
            iso = 2;    % Number of constituents to isolate (Calcification and 'Tissue')
            for con = 1:iso;
                if con == 1
                    
                    % HSL parameters to isolate BLACK pixels of CALCIFICATION
                    Hlow = 0;   Hup = 360;      % Hue
                    Slow = 0;   Sup = 1;        % Saturation

                    if strcmpi(stain,'VnK');
                        Llow = 0;   Lup = 0.45;     % Lightness
                    elseif strcmpi(stain,'AzR');
                        Llow = 0;   Lup = 0.3;      % Lightness
                    end
                    
                    C = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'k','white');
                    
                    [Havg_cal, Savg_cal, Lavg_cal] = HSLavg(C,H);
                                        
                elseif con == 2
                    
                    % HSL parameters to isolate PINK/RED pixels of 'TISSUE'
                    Hlow = 300;    Hup = 17;       % Hue
                    Slow = 0.1;    Sup = 1;        % Saturation
                    Llow = 0.1;    Lup = 0.93;     % Lightness
                    
                    O = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,C,'r','white');
                    
                    [Havg_tis, Savg_tis, Lavg_tis] = HSLavg(O,H);
                    
                end
                
            end
            
            % HSL parameters to isolate all NON-WHITE pixels
            Hlow = 0;   Hup = 360;      % Hue
            Slow = 0;   Sup = 1;        % Saturation
            Llow = 0;   Lup = 0.99;     % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
            
            % Compute number of identified constituent pixels
            pix_cal = length(find(C(:,:,3) < 1));      % Elastin
            pix_tis = length(find(O(:,:,3) < 1));      % Other Tissue
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_cal = pix_cal./pix_tot;      % Elastin
            af_tis = pix_tis./pix_tot;      % Other Tissue
            
            % Convert HSL images back to RGB for plotting/saving images
            Call = colorspace('HSL->RGB',C);
            Oall = colorspace('HSL->RGB',O);
            Tall = colorspace('HSL->RGB',T);
            
            % Plot/output routines for each constituent
            imstack = cat(4,Call,Oall,Tall);                                                                  % Images
            afrac   = cat(1,pix_cal,pix_tis,pix_tot,af_cal,af_tis)';                                          % Area fractions
            hslavg  = cat(1,cat(2,Havg_cal,Havg_tis),cat(2,Savg_cal,Savg_tis),cat(2,Lavg_cal,Lavg_tis));      % HSL values
            
            
            % ----- Masson's Trichrome (MTC) ------------------------------
        elseif strcmpi(stain,'MTC');
            
            iso = 2;    % Number of constituents to isolate (Cytoplasm and Collagen)
            for con = 1:iso;
                if con == 1
                    
                    % HSL parameters to isolate RED pixels of CYTOPLASM
                    Hlow = 250;   Hup = 25;       % Hue
                    Slow = 0.1;   Sup = 1;        % Saturation
                    
                    if strcmpi(tissue,'Myocardial Infarction         ');
                        Llow = 0.01;   Lup = 0.93;     % Lightness
                    else
                        Llow = 0.1;   Lup = 0.93;     % Lightness
                    end
                    
                    S = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'r','white');
                    
                    [Havg_cyt, Savg_cyt, Lavg_cyt] = HSLavg(S,H);
                    
                    if strcmpi(tissue,'Artery / Vasculature');
                        
                        % Blob analysis to identify and remove (adventitial)cell nuclei and other small pixel groups
                        cc = bwconncomp(imcomplement(S(:,:,3)),8);   % Mask of all pixel groups in image
                        
                        % Compute and store area, perimeter, and coordinate location of blobs in image
                        conarea  = regionprops(cc,'Area');          conarea = cell2mat(struct2cell(conarea));
                        conperim = regionprops(cc,'Perimeter');    conperim = cell2mat(struct2cell(conperim));
                        conloc   = regionprops(cc,'PixelList');     concirc = 4.*pi.*(conarea./(conperim + pi).^2);
                        
                        % Find small, round blobs as these are not SMCs
                        blob_loc = find(conarea < round(20*SFr*SFm)  | concirc > 0.6);
                        
                        % Find each pixel identified to be in a blob
                        for II = 1:length(blob_loc);
                            
                            num = length(conloc(blob_loc(II),1).PixelList(:,1));
                            
                            % Remove blobs from image (i.e., make blob pixels white)
                            for JJ = 1:num
                                S(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),1) = 360;
                                S(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),2) = 0;
                                S(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),3) = 1;
                            end
                            
                        end
                    end
                    
                elseif con == 2
                    
                    % HSL parameters to isolate BLUE pixels of COLLAGEN
                    Hlow = 150;    Hup = 250;     % Hue
                    Slow = 0.1;    Sup = 1.0;     % Saturation
                    Llow = 0.1;    Lup = 0.93;    % Lightness
                    
                    C = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,S,'b','white');
                    
                    [Havg_col, Savg_col, Lavg_col] = HSLavg(C,H);
                    
                end
            end
            
            % HSL parameters to isolate all NON-WHITE pixels
            Hlow = 0;   Hup = 360;      % Hue
            Slow = 0;   Sup = 1;        % Saturation
            Llow = 0;   Lup = 0.99;     % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
            
            % Compute number of identified constituent pixels
            pix_cyt = length(find(S(:,:,3) < 1));      % Cytoplasm
            pix_col = length(find(C(:,:,3) < 1));      % Collagen
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_cyt = pix_cyt./pix_tot;      % Cytoplasm
            af_col = pix_col./pix_tot;      % Collagen
            
            % Convert HSL images back to RGB for plotting/saving images
            Sall = colorspace('HSL->RGB',S);      % Cytoplasm
            Call = colorspace('HSL->RGB',C);      % Collagen
            Tall = colorspace('HSL->RGB',T);      % Total
            
            % Plot/output routines for each constituent
            imstack = cat(4,Sall,Call,Tall);                                                                % Images
            afrac   = cat(1,pix_cyt,pix_col,pix_tot,af_cyt,af_col)';                                        % Area fractions
            hslavg  = cat(1,cat(2,Havg_cyt,Havg_col),cat(2,Savg_cyt,Savg_col),cat(2,Lavg_cyt,Lavg_col));    % HSL values
            
            
            % ----- Picrosirius Red: Brightfield (bPSR) -------------------
        elseif strcmpi(stain,'bPSR');
            
            iso = 2;    % Number of constituents to isolate (Collagen and 'Tissue')
            for con = 1:iso;
                if con == 1
                    
                    % HSL parameters to isolate RED pixels of COLLAGEN
                    Hlow = 250;   Hup = 25;       % Hue
                    Slow = 0.1;   Sup = 1;        % Saturation
                    Llow = 0.01;   Lup = 0.93;    % Lightness
                    
                    C = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'r','white');
                    
                    [Havg_col, Savg_col, Lavg_col] = HSLavg(C,H);
                                        
                elseif con == 2
                    
                    % HSL parameters to isolate YELLOW pixels of TISSUE
                    Hlow = 25;     Hup = 75;     % Hue
                    Slow = 0.1;    Sup = 1.0;     % Saturation
                    Llow = 0.1;    Lup = 0.93;    % Lightness
                    
                    O = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,C,'y','white');
                    
                    [Havg_tis, Savg_tis, Lavg_tis] = HSLavg(O,H);
                    
                end
            end
            
            % HSL parameters to isolate all NON-WHITE pixels
            Hlow = 0;   Hup = 360;      % Hue
            Slow = 0;   Sup = 1;        % Saturation
            Llow = 0;   Lup = 0.99;     % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
            
            % Compute number of identified constituent pixels
            pix_col = length(find(C(:,:,3) < 1));      % Collagen
            pix_tis = length(find(O(:,:,3) < 1));      % 'Tissue'
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_col = pix_col./pix_tot;      % Collagen
            af_tis = pix_tis./pix_tot;      % 'Tissue'
            
            % Convert HSL images back to RGB for plotting/saving images
            Call = colorspace('HSL->RGB',C);      % Collagen
            Oall = colorspace('HSL->RGB',O);      % 'Tissue'
            Tall = colorspace('HSL->RGB',T);      % Total
            
            % Plot/output routines for each constituent
            imstack = cat(4,Call,Oall,Tall);                                                                % Images
            afrac   = cat(1,pix_col,pix_tis,pix_tot,af_col,af_tis)';                                        % Area fractions
            hslavg  = cat(1,cat(2,Havg_col,Havg_tis),cat(2,Savg_col,Savg_tis),cat(2,Lavg_col,Lavg_tis));    % HSL values
            
            
            % ----- Picrosirius Red: Darkfield (dPSR) ---------------------
        elseif strcmpi(stain,'dPSR')
            
            iso = 4;    % Number of constituents to isolate (Red, Orange, Yellow, and Green fibers)
            for con = 1:iso
                if con == 1
                    
                    % HSL parameters to isolate RED pixels of COLLAGEN FIBERS
                    Hlow = 324;     Hup = 12;       % Hue
                    Slow = 0.1;     Sup = 1;        % Saturation
                    Llow = 0.10;    Lup = 0.93;     % Lightness    DSL 2022
%                     Llow = 0.06;    Lup = 0.93;     % Lightness
                    
                    R = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'r','black');
                    
                    [Havg_red, Savg_red, Lavg_red] = HSLavg(R,H);
                    
                elseif con == 2
                    
                    % HSL parameters to isolate ORANGE pixels of COLLAGEN FIBERS
                    Hlow = 13;     Hup = 52;      % Hue
                    Slow = 0.1;    Sup = 1;       % Saturation
                    Llow = 0.10;    Lup = 0.93;     % Lightness    DSL 2022
%                     Llow = 0.06;    Lup = 0.93;     % Lightness
                    
                    O = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,R,'o','black');
                    
                    [Havg_orn, Savg_orn, Lavg_orn] = HSLavg(O,H);
                    
                elseif con == 3
                    
                    % HSL parameters to isolate YELLOW pixels of COLLAGEN FIBERS
                    Hlow = 53;     Hup = 72;       % Hue
                    Slow = 0.1;    Sup = 1;        % Saturation
                    Llow = 0.10;    Lup = 0.93;     % Lightness    DSL 2022
%                     Llow = 0.06;    Lup = 0.93;     % Lightness
                    
                    Y = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(R+O-1),'y','black');
                    
                    [Havg_yel, Savg_yel, Lavg_yel] = HSLavg(Y,H);
                    
                elseif con == 4
                    
                    % HSL parameters to isolate GREEN pixels of COLLAGEN FIBERS
                    Hlow = 73;     Hup = 180;      % Hue
                    Slow = 0.1;    Sup = 1;        % Saturation
                    Llow = 0.10;    Lup = 0.93;     % Lightness    DSL 2022
%                     Llow = 0.06;    Lup = 0.93;     % Lightness
                    
                    G = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(R+O+Y-2),'g','black');
                    
                    [Havg_grn, Savg_grn, Lavg_grn] = HSLavg(G,H);
                    
                end
            end
            
            % HSL parameters to isolate all NON-BLACK pixels
            Hlow = 0;      Hup = 360;   % Hue
            Slow = 0;      Sup = 1;     % Saturation
            Llow = 0.10;   Lup = 1;     % Lightness    DSL 2022
            % Llow = 0.01;   Lup = 1;     % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'k','black');
            
            % Compute number of identified constituent pixels
            pix_red = length(find(R(:,:,3) < 1));      % Red
            pix_orn = length(find(O(:,:,3) < 1));      % Orange
            pix_yel = length(find(Y(:,:,3) < 1));      % Yellow
            pix_grn = length(find(G(:,:,3) < 1));      % Green
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_red = pix_red./pix_tot;      % Red
            af_orn = pix_orn./pix_tot;      % Orange
            af_yel = pix_yel./pix_tot;      % Yellow
            af_grn = pix_grn./pix_tot;      % Green
            
            % Convert HSL images back to RGB for plotting/saving images
            Rall = colorspace('HSL->RGB',R);
            Oall = colorspace('HSL->RGB',O);
            Yall = colorspace('HSL->RGB',Y);
            Gall = colorspace('HSL->RGB',G);
            Tall = colorspace('HSL->RGB',T);
            
            % Plot/output routines for each constituent
            imstack = cat(4,Rall,Oall,Yall,Gall,Tall);                                                                                                              % Images
            afrac   = cat(1,pix_red,pix_orn,pix_yel,pix_grn,pix_tot,af_red,af_orn,af_yel,af_grn)';                                                                  % Area fractions
            hslavg  = cat(1,cat(2,Havg_red,Havg_orn,Havg_yel,Havg_grn),cat(2,Savg_red,Savg_orn,Savg_yel,Savg_grn),cat(2,Lavg_red,Lavg_orn,Lavg_yel,Lavg_grn));      % HSL values
            
            
            % ----- Movat's Pentachrome (MOV) -----------------------------
        elseif strcmpi(stain,'MOV');
            
            iso = 5;    % Number of constituents to isolate (Ground Substance (GAG), Fibrin, Elastin, Cytoplasm, Collagen)
            for con = 1:iso
                if con == 2    %%%% 1
                    
                    % HSL parameters to isolate BLUE-GREEN pixels of GROUND SUBSTANCE (GAG)
%                     Hlow = 173;    Hup = 240;      % Hue
%                     Slow = 0.09;   Sup = 1;        % Saturation
%                     Llow = 0.1;    Lup = 0.8;      % Lightness
                    
                    % HSL parameters to isolate BLUE-GREEN pixels of GROUND SUBSTANCE (GAG) DSL 2021
%                     Hlow = 140;    Hup = 220;      % Hue
%                     Slow = 0.10;   Sup = 0.35;     % Saturation
%                     Llow = 0.51;   Lup = 0.76;     % Lightness
                    
                    % HSL parameters to isolate BLUE-GREEN pixels of GROUND SUBSTANCE (GAG) DW Aging
%                     Hlow = 145;    Hup = 226;      % Hue
%                     Slow = 0.12;   Sup = 0.87;     % Saturation
%                     Llow = 0.22;   Lup = 0.70;     % Lightness
                    
                    % HSL parameters to isolate BLUE-GREEN pixels of GROUND SUBSTANCE (GAG) C1041G+BAPN 2022
%                     Hlow = 140;    Hup = 245;      % Hue, old low=159
%                     Slow = 0.12;   Sup = 0.87;     % Saturation
%                     Llow = 0.31;   Lup = 0.75;     % Lightness
                    
                    % HSL parameters to isolate BLUE-GREEN pixels of GROUND SUBSTANCE (GAG) WT vs mGR
%                     Hlow = 130;    Hup = 226;      % WT, a5/2, C1041G Hue, hPA
%                     Slow = 0.09;   Sup = 0.99;     % WT, a5/2 Saturation, hPA
%                     Hlow = 130;    Hup = 226;      % C1041G ATA Hue
%                     Slow = 0.05;   Sup = 0.99;     % C1041G ATA Saturation
                    Hlow = 138;    Hup = 226;      % mgR ATA Hue
                    Slow = 0.06;   Sup = 0.99;     % mgR ATA Saturation
%                     Hlow = 135;    Hup = 226;      % mgR DTA Hue
%                     Slow = 0.05;   Sup = 0.99;     % mgR DTA Saturation
                    Llow = 0.22;   Lup = 0.81;     % Lightness
                    
%                     G = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'b','white');
                    G = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,E,'b','white');
                    
                    [Havg_gag, Savg_gag, Lavg_gag] = HSLavg(G,H);
                    
                elseif con == 3    %%%% 4
                    
                    % HSL parameters to isolate BRIGHT RED pixels of FIBRIN
%                     Hlow = 262;    Hup = 11;       % Hue
%                     Slow = 0.2;    Sup = 1.0;      % Saturation
%                     Llow = 0.1;    Lup = 0.93;     % Lightness
                    
                    % HSL parameters to isolate BRIGHT RED pixels of FIBRIN DW Aging
%                     Hlow = 340;    Hup = 11;       % Hue
%                     Slow = 0.02;   Sup = 0.99;     % Saturation
%                     Llow = 0.27;   Lup = 0.93;     % Lightness
                    
                    % HSL parameters to isolate BRIGHT RED pixels of FIBRIN C1041G+BAPN 2022
%                     Hlow = 280;     Hup = 45;       % Hue
%                     Slow = 0.05;    Sup = 0.83;     % Saturation
%                     Llow = 0.35;    Lup = 0.68;     % Lightness
                    
                    % HSL parameters to isolate BRIGHT RED pixels of FIBRIN WT vs mGR, hPA
                    Hlow = 282;    Hup = 11;       % Hue
                    Slow = 0.37;   Sup = 0.70;     % Saturation
                    Llow = 0.32;   Lup = 0.93;     % Lightness
                    
%                     F = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,(C+E+G-2),'m','white');
%                     F = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,(M+E+G-2),'m','white');
                    F = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,(E+G-1),'m','white');
                    
                    [Havg_fib, Savg_fib, Lavg_fib] = HSLavg(F,H);
                    
                    % Blob analysis to identify and remove cell nuclei and other small pixel groups
                    %     Note: Cell nuclei appear in VVG stains and look like Elastin when stained
                    
                    F(:,:,3) = imcomplement(imclose(imcomplement(F(:,:,3)),strel('disk',1)));
                    
                    cc = bwconncomp(imcomplement(F(:,:,3)),8);   % Mask of all pixel groups in image
                    
                    % Compute morphology of pixels in each blob and write data from structure to array
                    conarea  = regionprops(cc,'Area');          conarea = cell2mat(struct2cell(conarea));
                    conperim = regionprops(cc,'Perimeter');    conperim = cell2mat(struct2cell(conperim));
                    conloc   = regionprops(cc,'PixelList');     concirc = 4.*pi.*(conarea./(conperim + pi).^2);
                    
                    % Find small, round blobs as these are not Fibrin
                    blob_loc = find((conarea > round(20*SFr*SFm) & conarea < round(170*SFr*SFm)) | concirc > 0.6);
                    
                    % Find each pixel identified to be in a blob
                    for II = 1:length(blob_loc);
                        
                        num = length(conloc(blob_loc(II),1).PixelList(:,1));
                        
                        % Remove blobs from image (i.e., make blob pixels white)
                        for JJ = 1:num
                            F(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),1) = 360;
                            F(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),2) = 0;
                            F(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),3) = 1;
                        end
                        
                    end
                    
                elseif con == 1    %%%% 2
                    
                    % HSL parameters to isolate BLACK pixels of ELASTIN
                    Hlow = 0;   Hup = 360;      % Hue
                    Slow = 0;   Sup = 1;        % Saturation
                    
                    if strcmpi(enhance,'Yes')
                        Llow = 0;   Lup = 0.16;     % Lightness
                    else
                        Llow = 0;   Lup = 0.24;     % Lightness DSL
%                         Llow = 0;   Lup = 0.28;     % Lightness DW Aging
                    end
                    
%                     E = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,G,'k','white');
                    E = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'k','white');
                    
                    [Havg_eln, Savg_eln, Lavg_eln] = HSLavg(E,H);
                    
                    if strcmpi(tissue,'Artery / Vasculature');
                        
                        % Blob analysis to identify and remove cell nuclei and other small pixel groups
                        %     Note: Cell nuclei appear in VVG stains and look like Elastin when stained
                        
                        E(:,:,3) = imcomplement(imclose(imcomplement(E(:,:,3)),strel('disk',1)));
                        cc = bwconncomp(imcomplement(E(:,:,3)),8);   % Mask of all pixel groups in image
                        
                        % Compute and store area, perimeter, and coordinate location of blobs in image
                        conarea  = regionprops(cc,'Area');          conarea = cell2mat(struct2cell(conarea));
                        conperim = regionprops(cc,'Perimeter');    conperim = cell2mat(struct2cell(conperim));
                        conloc   = regionprops(cc,'PixelList');     concirc = 4.*pi.*(conarea./(conperim + pi).^2);
                        
                        % Find small, round blobs as these are not Elastin
                        if ~isempty(LDX) && ~strcmpi(layernm{LDX,:},'adventitia')
                            blob_loc = find(conarea < round(50*SFr*SFm) | concirc > 0.6);
                        else
                            blob_loc = find(concirc > 0.6);
                        end
                        
                        % Find each pixel identified to be in a blob
                        for II = 1:length(blob_loc);
                            
                            num = length(conloc(blob_loc(II),1).PixelList(:,1));
                            
                            % Remove blobs from image (i.e., make blob pixels white)
                            for JJ = 1:num
                                E(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),1) = 360;
                                E(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),2) = 0;
                                E(conloc(blob_loc(II),1).PixelList(JJ,2),conloc(blob_loc(II),1).PixelList(JJ,1),3) = 1;
                            end
                            
                        end
                        
                        % High-resolution bleb removal
                        if ~isempty(LDX) && ~strcmpi(layernm{LDX,:},'adventitia') && SFm*SFr > 1
                            
                            if strcmpi(rembleb,'Yes') || strcmpi(rembleb,'No')
                                qoption = {'Yes','No','Yes to All','No to All'};
                                qans = genquestdlg(qoption, 'Blebs',qoption,'Perform Bleb Removal?');
                                rembleb = qoption{qans};
                            end
                            
                            if strcmpi(rembleb,'Yes') || strcmpi(rembleb,'Yes to All');    E = removeblebs(E,0);    end
                            
                        end
                    end
                    
                elseif con == 4    %%%% 5
                    
                    % HSL parameters to isolate RED pixels of CYTOPLASM
%                     Hlow = 270;     Hup = 357;      % Hue
%                     Slow = 0.12;    Sup = 0.85;     % Saturation
%                     Llow = 0.1;     Lup = 0.93;     % Lightness
                    
                    % HSL parameters to isolate RED pixels of CYTOPLASM DW Aging
%                     Hlow = 226;     Hup =  335;      % Hue
%                     Slow = 0.02;    Sup = 1.00;     % Saturation
%                     Llow = 0.22;    Lup = 0.70;     % Lightness
                    
                    % HSL parameters to isolate RED pixels of CYTOPLASM C1041G+BAPN 2022
%                     Hlow = 230;     Hup = 355;      % Hue
%                     Slow = 0.02;    Sup = 0.87;     % Saturation
%                     Llow = 0.24;    Lup = 0.70;     % Lightness
                    
                    % HSL parameters to isolate RED pixels of CYTOPLASM WT vs mgR 2022
                    % Hlow =  195;    Hup =   35;     % WT Hue
                    % Slow = 0.03;    Sup = 0.87;     % WT Saturation
                    Hlow =  195;    Hup =   40;     % a5/2 C1041G Hue
                    Slow = 0.01;    Sup = 0.87;     % a5/2 C1041G Saturation
%                     Hlow =  195;    Hup =   40;     % mgR ATA Hue
%                     Slow = 0.05;    Sup = 0.87;     % mgR ATA Saturation
%                     Hlow =  195;    Hup =   40;     % mgR DTA Hue
%                     Slow = 0.04;    Sup = 0.87;     % mgR DTA Saturation
                    Llow = 0.19;    Lup = 0.92;     % Lightness
%                     Hlow =  195;    Hup =  355;     % hPA Hue
%                     Slow = 0.05;    Sup = 0.87;     % hPA Saturation
%                     Llow = 0.19;    Lup = 0.92;     % hPA Lightness					
                    
%                     M = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(C+E+G+F-3),'r','white');
%                     M = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(E+G-1),'r','white');
                    M = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(E+G+F-2),'r','white');
                    
                    [Havg_cyt, Savg_cyt, Lavg_cyt] = HSLavg(M,H);
                    
                elseif con == 5    %%%% 3
                    
                    % HSL parameters to isolate all NON-WHITE pixels
                    Hlow = 0;   Hup = 360;      % Hue
                    Slow = 0;   Sup = 1;        % Saturation
                    Llow = 0;   Lup = 0.99;     % Lightness
                    
                    T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
                    
                    % HSL parameters to isolate all FIBRIN pixels
%                     Hlow = 262;    Hup = 11;       % Hue
%                     Slow = 0.2;    Sup = 1.0;      % Saturation
%                     Llow = 0.1;    Lup = 0.93;     % Lightness
                    
                    % HSL parameters to isolate BRIGHT RED pixels of FIBRIN DW Aging
%                     Hlow = 340;    Hup = 11;       % Hue
%                     Slow = 0.02;   Sup = 0.99;     % Saturation
%                     Llow = 0.27;   Lup = 0.93;     % Lightness
                    
                    % HSL parameters to isolate BRIGHT RED pixels of FIBRIN C1041G+BAPN 2022
%                     Hlow = 280;     Hup = 45;       % Hue
%                     Slow = 0.05;    Sup = 0.83;     % Saturation
%                     Llow = 0.35;    Lup = 0.68;     % Lightness
                    
                    % HSL parameters to isolate BRIGHT RED pixels of FIBRIN WT vs mGR, hPA
                    Hlow = 282;    Hup = 11;       % Hue
                    Slow = 0.37;   Sup = 0.70;     % Saturation
                    Llow = 0.22;   Lup = 0.93;     % Lightness
                    
                    F = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,(E+G-1),'m','white');
%                     F = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,(E+G-1),'m','white');
                    
                    % Compute number of identified constituent pixels
                    pix_fib = length(find(F(:,:,3) < 1));      % Fibrin
                    pix_tot = length(find(T(:,:,3) < 1));      % Total
            
                    % HSL parameters to isolate YELLOW-BROWN pixels of COLLAGEN (part 1)
                    Hlow = 0;      Hup = 360;     % Hue
                    Slow = 0.0;    Sup = (1+(pix_fib/pix_tot))*0.2;     % Saturation
                    Llow = 0.0;    Lup = 0.99;     % Lightness
                    
%                     C1 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(E+G-1),'y','white');
%                     C1 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(M+E+G-2),'y','white');
                    C1 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(M+E+G+F-3),'y','white');

                    % HSL parameters to isolate YELLOW-BROWN pixels of COLLAGEN (part 3)
%                     Hlow = 30;     Hup = 88;       % Hue
%                     Slow = 0.0;    Sup = 1;        % Saturation
%                     Llow = 0.1;    Lup = 0.93;     % Lightness
                    
                    % HSL parameters to isolate YELLOW-BROWN pixels of COLLAGEN (part 3) C1041G+BAPN & DW Aging 2022
%                     Hlow = 300;    Hup = 90;       % Hue
%                     Slow = 0.1;    Sup = 0.27;     % Saturation
%                     Llow = 0.32;   Lup = 0.93;     % Lightness
                    
                    % HSL parameters to isolate YELLOW-BROWN pixels of COLLAGEN (part 3) WT vs mgR 2022, hPA
                    Hlow =   27;    Hup = 120;      % Hue
                    Slow = 0.04;    Sup = 0.82;     % Saturation
                    Llow = 0.32;    Lup = 0.75;     % Lightness
                    
%                     C3 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(C1+E+G-2),'y','white');
                    C3 = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,abs(C1+M+E+G-3),'y','white');
                    
                    % Compare color isolations to identify YELLOW-BROWN pixels of COLLAGEN based on saturation level
                    Cmsk = C1(:,:,3) + C3(:,:,3);
                    
                    C = zeros(size(C1));
                    
                    C1 = C(:,:,1) + 360;    C1(Cmsk ~= 2) = 60;
                    C2 = C(:,:,1);          C2(Cmsk ~= 2) = 1;
                    C3 = C(:,:,1) + 1;      C3(Cmsk ~= 2) = 0.5;
                    
                    C = cat(3,C1,C2,C3);
                    [Havg_col, Savg_col, Lavg_col] = HSLavg(C,H);
                                      
                end
                
            end
            
            % HSL parameters to isolate all NON-WHITE pixels
            Hlow = 0;   Hup = 360;      % Hue
            Slow = 0;   Sup = 1;        % Saturation
            Llow = 0;   Lup = 0.99;     % Lightness
            
            T = HSLfilter([Hlow Hup],[Slow Sup],[Llow Lup],H,[],'g','white');
            
            % Compute number of identified constituent pixels
            pix_gag = length(find(G(:,:,3) < 1));      % Ground Substance (GAG)
            pix_fib = length(find(F(:,:,3) < 1));      % Fibrin
            pix_eln = length(find(E(:,:,3) < 1));      % Elastin
            pix_cyt = length(find(M(:,:,3) < 1));      % Cytoplasm
            pix_col = length(find(C(:,:,3) < 1));      % Collagen
            pix_tot = length(find(T(:,:,3) < 1));      % Total
            
            % Compute area fractions for each constituent
            af_gag = pix_gag./pix_tot;      % Ground Substance (GAG)
            af_fib = pix_fib./pix_tot;      % Fibrin
            af_eln = pix_eln./pix_tot;      % Elastin
            af_cyt = pix_cyt./pix_tot;      % Cytoplasm
            af_col = pix_col./pix_tot;      % Collagen
            
            % Convert HSL images back to RGB for plotting/saving images
            Gall = colorspace('HSL->RGB',G);      % Ground Substance (GAG)
            Fall = colorspace('HSL->RGB',F);      % Fibrin
            Eall = colorspace('HSL->RGB',E);      % Elastin
            Mall = colorspace('HSL->RGB',M);      % Cytoplasm
            Call = colorspace('HSL->RGB',C);      % Collagen
            Tall = colorspace('HSL->RGB',T);      % Total
            
            % Plot/output routines for each constituent
            imstack = cat(4,Gall,Fall,Eall,Mall,Call,Tall);                                                                                                                                   % Images
            afrac = cat(1,pix_gag,pix_fib,pix_eln,pix_cyt,pix_col,pix_tot,af_gag,af_fib,af_eln,af_cyt,af_col)';                                                                               % Area fractions
            hslavg = cat(1,cat(2,Havg_gag,Havg_fib,Havg_eln,Havg_cyt,Havg_col),cat(2,Savg_gag,Savg_fib,Savg_eln,Savg_cyt,Savg_col),cat(2,Lavg_gag,Lavg_fib,Lavg_eln,Lavg_cyt,Lavg_col));      % HSL values
            
        end
        
        % ===== LOCAL color extraction ===========================================
    else                                  % partitions successfully defined
        
        % ----- Brightfield Immunohistochemistry (IHC) --------------------
        if strcmpi(stain,'IHC');
            
            % Load images
            if layer && ~isempty(LDX)
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\IHC_',fname,'_ihc-',layernm{LDX,:},'.tif'));      HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\IHC_',fname,'_total-',layernm{LDX,:},'.tif'));    HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            else
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\IHC_',fname,'_ihc.tif'));      HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\IHC_',fname,'_total.tif'));    HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            end
            
            % Count total pixels
            pix_tot = length(find(HT3 < 1));
            
            % Initialize local area fraction / hsl storage
            npartq = size(localpart,1);     npartr = size(localpart,2);      [xsz,ysz,~] = size(HT(:,:,3));
            afrac = cell(npartq,npartr);    hslavg = cell(npartq,npartr);
            
            h = waitbar(0,'Calculating Local Area Fractions...','Name',char(strcat(stain,{' '},'/',{' '},outnm{IDX,:})));
            
            % Calculate local area fractions
            for II = 1:npartq;
                for JJ = 1:npartr;
                    
                    waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                    
                    % Generate mask from current partition
                    Imsk = false(xsz,ysz);
                    for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                    Imsk = imfill(Imsk,'holes');
                    
                    % Extract all non-background pixels
                    pix_loc = length(find(HT3(Imsk == 1) < 1));
                    pix_ihc = length(find(HC3(Imsk == 1) < 1));    af_ihc = pix_ihc/pix_loc;    af_ihc_tot = pix_ihc/pix_tot;
                    
                    % Extract HSL of positive pixels
                    Htemp = H1(Imsk == 1 & HC3 < 1);             Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_ihc = mean(Htemp) + 360;      % IHC
                    Savg_ihc = mean(H2(Imsk == 1 & HC3 < 1));    Lavg_ihc = mean(H3(Imsk == 1 & HC3 < 1));
                    
                    afrac{II,JJ}  = cat(2,pix_ihc,pix_loc,af_ihc,af_ihc_tot);
                    hslavg{II,JJ} = cat(1,Havg_ihc,Savg_ihc,Lavg_ihc);
                    
                end
            end
            close(h)
            
            % Compile images for saving
            imstack = cat(4,C,T);
            
            
            % ----- Hematoxylin & Eosin (H&E) -----------------------------
        elseif strcmpi(stain,'H&E');
            
            fprintf('H&E ANALYSIS IS NOT YET FUNCTIONAL!!!!!\n');
            imstack = NaN;    afrac = NaN;    hslavg = NaN;
            return
            
            % ----- Verhoeff Van Gieson (VVG) -----------------------------
        elseif strcmpi(stain,'VVG');
            
            % Load images
            if layer && ~isempty(LDX)
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                E = imread(strcat(path,groupnm,'\',fname,'\VVG_',fname,'_elastin-',layernm{LDX,:},'.tif'));      HE = colorspace('RGB->HSL',E);     HE3 = HE(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\VVG_',fname,'_tissue-',layernm{LDX,:},'.tif'));       HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\VVG_',fname,'_total-',layernm{LDX,:},'.tif'));        HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            else
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                E = imread(strcat(path,groupnm,'\',fname,'\VVG_',fname,'_elastin.tif'));      HE = colorspace('RGB->HSL',E);     HE3 = HE(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\VVG_',fname,'_tissue.tif'));       HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\VVG_',fname,'_total.tif'));        HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            end
            
            % Count total pixels
            pix_tot = length(find(HT3 < 1));
            
            % Initialize local area fraction / hsl storage
            npartq = size(localpart,1);     npartr = size(localpart,2);      [xsz,ysz,~] = size(HT(:,:,3));
            afrac = cell(npartq,npartr);    hslavg = cell(npartq,npartr);
            
            h = waitbar(0,'Calculating Local Area Fractions...','Name',char(strcat(stain,{' '},'/',{' '},outnm{IDX,:})));
            
            % Calculate local area fractions
            for II = 1:npartq;
                for JJ = 1:npartr;
                    
                    waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                    
                    % Generate mask from current partition
                    Imsk = false(xsz,ysz);
                    for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                    Imsk = imfill(Imsk,'holes');
                    
                    % Extract all non-background pixels
                    pix_loc = length(find(HT3(Imsk == 1) < 1));
                    pix_eln = length(find(HE3(Imsk == 1) < 1));    af_eln = pix_eln/pix_loc;    af_eln_tot = pix_eln/pix_tot;
                    pix_tis = length(find(HO3(Imsk == 1) < 1));    af_tis = pix_tis/pix_loc;    af_tis_tot = pix_tis/pix_tot;
                    
                    % Extract HSL of positive pixels
                    Htemp = H1(Imsk == 1 & HE3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_eln = mean(Htemp) + 360;      % Elastin
                    Htemp = H1(Imsk == 1 & HO3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_tis = mean(Htemp) + 360;      % Tissue
                    
                    Savg_eln = mean(H2(Imsk == 1 & HE3 < 1));    Lavg_eln = mean(H3(Imsk == 1 & HE3 < 1));      % Elastin
                    Savg_tis = mean(H2(Imsk == 1 & HO3 < 1));    Lavg_tis = mean(H3(Imsk == 1 & HO3 < 1));      % Tissue
                    
                    afrac{II,JJ}  = cat(2,pix_eln,pix_tis,pix_loc,af_eln,af_tis,af_eln_tot,af_tis_tot);
                    hslavg{II,JJ} = cat(1,cat(2,Havg_eln,Havg_tis),cat(2,Savg_eln,Savg_tis),cat(2,Lavg_eln,Lavg_tis));
                    
                end
            end
            close(h)
            
            % Compile images for saving
            imstack = cat(4,E,O,T);
            
            
            % ----- Calcification -----------------------------
        elseif strcmpi(stain,'VnK') || strcmpi(stain,'AzR');
            
            % Load images
            if layer && ~isempty(LDX)
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_calcification-',layernm{LDX,:},'.tif'));      HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_tissue-',layernm{LDX,:},'.tif'));             HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_total-',layernm{LDX,:},'.tif'));              HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            else
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_calcification.tif'));      HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_tissue.tif'));             HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_total.tif'));              HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            end
            
            % Count total pixels
            pix_tot = length(find(HT3 < 1));
            
            % Initialize local area fraction / hsl storage
            npartq = size(localpart,1);     npartr = size(localpart,2);      [xsz,ysz,~] = size(HT(:,:,3));
            afrac = cell(npartq,npartr);    hslavg = cell(npartq,npartr);
            
            h = waitbar(0,'Calculating Local Area Fractions...','Name',char(strcat(stain,{' '},'/',{' '},outnm{IDX,:})));
            
            % Calculate local area fractions
            for II = 1:npartq;
                for JJ = 1:npartr;
                    
                    waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                    
                    % Generate mask from current partition
                    Imsk = false(xsz,ysz);
                    for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                    Imsk = imfill(Imsk,'holes');
                    
                    % Extract all non-background pixels
                    pix_loc = length(find(HT3(Imsk == 1) < 1));
                    pix_cal = length(find(HC3(Imsk == 1) < 1));    af_cal = pix_cal/pix_loc;    af_cal_tot = pix_cal/pix_tot;
                    pix_tis = length(find(HO3(Imsk == 1) < 1));    af_tis = pix_tis/pix_loc;    af_tis_tot = pix_tis/pix_tot;
                    
                    % Extract HSL of positive pixels
                    Htemp = H1(Imsk == 1 & HC3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_cal = mean(Htemp) + 360;      % Calcification
                    Htemp = H1(Imsk == 1 & HO3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_tis = mean(Htemp) + 360;      % 'Tissue'
                    
                    Savg_cal = mean(H2(Imsk == 1 & HC3 < 1));    Lavg_cal = mean(H3(Imsk == 1 & HC3 < 1));      % Calcification
                    Savg_tis = mean(H2(Imsk == 1 & HO3 < 1));    Lavg_tis = mean(H3(Imsk == 1 & HO3 < 1));      % 'Tissue'
                    
                    afrac{II,JJ}  = cat(2,pix_cal,pix_tis,pix_loc,af_cal,af_tis,af_cal_tot,af_tis_tot);
                    hslavg{II,JJ} = cat(1,cat(2,Havg_cal,Havg_tis),cat(2,Savg_cal,Savg_tis),cat(2,Lavg_cal,Lavg_tis));
                    
                end
            end
            close(h)
            
            % Compile images for saving
            imstack = cat(4,C,O,T);
            
            
            % ----- Masson's Trichrome (MTC) ------------------------------
        elseif strcmpi(stain,'MTC');
            
            % Load images
            if layer && ~isempty(LDX)
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                S = imread(strcat(path,groupnm,'\',fname,'\MTC_',fname,'_cytoplasm-',layernm{LDX,:},'.tif'));    HS = colorspace('RGB->HSL',S);     HS3 = HS(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\MTC_',fname,'_collagen-',layernm{LDX,:},'.tif'));     HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\MTC_',fname,'_total-',layernm{LDX,:},'.tif'));        HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            else
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                S = imread(strcat(path,groupnm,'\',fname,'\MTC_',fname,'_cytoplasm.tif'));    HS = colorspace('RGB->HSL',S);     HS3 = HS(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\MTC_',fname,'_collagen.tif'));     HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\MTC_',fname,'_total.tif'));        HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            end
            
            % Count total pixels
            pix_tot = length(find(HT3 < 1));
            
            % Initialize local area fraction / hsl storage
            npartq = size(localpart,1);     npartr = size(localpart,2);      [xsz,ysz,~] = size(HT(:,:,3));
            afrac = cell(npartq,npartr);    hslavg = cell(npartq,npartr);
            
            h = waitbar(0,'Calculating Local Area Fractions...','Name',char(strcat(stain,{' '},'/',{' '},outnm{IDX,:})));
            
            % Calculate local area fractions
            for II = 1:npartq;
                for JJ = 1:npartr;
                    
                    waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                    
                    % Generate mask from current partition
                    Imsk = false(xsz,ysz);
                    for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                    Imsk = imfill(Imsk,'holes');
                    
                    % Extract all non-background pixels
                    pix_loc = length(find(HT3(Imsk == 1) < 1));
                    pix_cyt = length(find(HS3(Imsk == 1) < 1));    af_cyt = pix_cyt/pix_loc;    af_cyt_tot = pix_cyt/pix_tot;
                    pix_col = length(find(HC3(Imsk == 1) < 1));    af_col = pix_col/pix_loc;    af_col_tot = pix_col/pix_tot;
                    
                    % Extract HSL of positive pixels
                    Htemp = H1(Imsk == 1 & HS3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_cyt = mean(Htemp) + 360;      % Cytoplasm
                    Htemp = H1(Imsk == 1 & HC3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_col = mean(Htemp) + 360;      % Collagen
                    
                    Savg_cyt = mean(H2(Imsk == 1 & HS3 < 1));    Lavg_cyt = mean(H3(Imsk == 1 & HS3 < 1));      % Cytoplasm
                    Savg_col = mean(H2(Imsk == 1 & HC3 < 1));    Lavg_col = mean(H3(Imsk == 1 & HC3 < 1));      % Collagen
                    
                    afrac{II,JJ}  = cat(2,pix_cyt,pix_col,pix_loc,af_cyt,af_col,af_cyt_tot,af_col_tot);
                    hslavg{II,JJ} = cat(1,cat(2,Havg_cyt,Havg_col),cat(2,Savg_cyt,Savg_col),cat(2,Lavg_cyt,Lavg_col));
                    
                end
            end
            close(h)
            
            % Compile images for saving
            imstack = cat(4,S,C,T);
            
            
            % ----- Picrosirius Red: Brightfield (bPSR) ---------------------
        elseif strcmpi(stain,'bPSR');
            
            % Load images
            if layer && ~isempty(LDX)
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\bPSR_',fname,'_collagen-',layernm{LDX,:},'.tif'));       HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\bPSR_',fname,'_tissue-',layernm{LDX,:},'.tif'));         HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\bPSR_',fname,'_total-',layernm{LDX,:},'.tif'));          HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            else
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\bPSR_',fname,'_collagen.tif'));       HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\bPSR_',fname,'_tissue.tif'));         HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\bPSR_',fname,'_total.tif'));          HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            end
            
            % Count total pixels
            pix_tot = length(find(HT3 < 1));
            
            % Initialize local area fraction / hsl storage
            npartq = size(localpart,1);     npartr = size(localpart,2);      [xsz,ysz,~] = size(HT(:,:,3));
            afrac = cell(npartq,npartr);    hslavg = cell(npartq,npartr);
            
            h = waitbar(0,'Calculating Local Area Fractions...','Name',char(strcat(stain,{' '},'/',{' '},outnm{IDX,:})));
            
            % Calculate local area fractions
            for II = 1:npartq;
                for JJ = 1:npartr;
                    
                    waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                    
                    % Generate mask from current partition
                    Imsk = false(xsz,ysz);
                    for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                    Imsk = imfill(Imsk,'holes');
                    
                    % Extract all non-background pixels
                    pix_loc = length(find(HT3(Imsk == 1) < 1));
                    pix_col = length(find(HC3(Imsk == 1) < 1));    af_col = pix_col/pix_loc;    af_col_tot = pix_col/pix_tot;
                    pix_tis = length(find(HO3(Imsk == 1) < 1));    af_tis = pix_tis/pix_loc;    af_tis_tot = pix_tis/pix_tot;
                    
                    % Extract HSL of positive pixels
                    Htemp = H1(Imsk == 1 & HC3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_col = mean(Htemp) + 360;      % Red
                    Htemp = H1(Imsk == 1 & HO3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_tis = mean(Htemp) + 360;      % Orange
                    
                    Savg_col = mean(H2(Imsk == 1 & HC3 < 1));    Lavg_col = mean(H3(Imsk == 1 & HC3 < 1));      % Red
                    Savg_tis = mean(H2(Imsk == 1 & HO3 < 1));    Lavg_tis = mean(H3(Imsk == 1 & HO3 < 1));      % Orange
                    
                    afrac{II,JJ}  = cat(2,pix_col,pix_tis,pix_loc,af_col,af_tis,af_col_tot,af_tis_tot);
                    hslavg{II,JJ} = cat(1,cat(2,Havg_col,Havg_tis),cat(2,Savg_col,Savg_tis),cat(2,Lavg_col,Lavg_tis));
                    
                end
            end
            close(h)
            
            % Compile images for saving
            imstack = cat(4,C,O,T);
            
            
            % ----- Picrosirius Red: Darkfield (dPSR) ---------------------
        elseif strcmpi(stain,'dPSR');
            
            % Load images
            if layer && ~isempty(LDX)
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                R = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_red-',layernm{LDX,:},'.tif'));       HR = colorspace('RGB->HSL',R);     HR3 = HR(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_orange-',layernm{LDX,:},'.tif'));    HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                Y = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_yellow-',layernm{LDX,:},'.tif'));    HY = colorspace('RGB->HSL',Y);     HY3 = HY(:,:,3);
                G = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_green-',layernm{LDX,:},'.tif'));     HG = colorspace('RGB->HSL',G);     HG3 = HG(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_total-',layernm{LDX,:},'.tif'));     HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            else
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                R = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_red.tif'));       HR = colorspace('RGB->HSL',R);     HR3 = HR(:,:,3);
                O = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_orange.tif'));    HO = colorspace('RGB->HSL',O);     HO3 = HO(:,:,3);
                Y = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_yellow.tif'));    HY = colorspace('RGB->HSL',Y);     HY3 = HY(:,:,3);
                G = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_green.tif'));     HG = colorspace('RGB->HSL',G);     HG3 = HG(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\dPSR_',fname,'_total.tif'));     HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            end
            
            % Count total pixels
            pix_tot = length(find(HT3 < 1));
            
            % Initialize local area fraction / hsl storage
            npartq = size(localpart,1);     npartr = size(localpart,2);      [xsz,ysz,~] = size(HT(:,:,3));
            afrac = cell(npartq,npartr);    hslavg = cell(npartq,npartr);
            
            h = waitbar(0,'Calculating Local Area Fractions...','Name',char(strcat(stain,{' '},'/',{' '},outnm{IDX,:})));
            
            % Calculate local area fractions
            for II = 1:npartq;
                for JJ = 1:npartr;
                    
                    waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                    
                    % Generate mask from current partition
                    Imsk = false(xsz,ysz);
                    for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                    Imsk = imfill(Imsk,'holes');
                    
                    % Extract all non-background pixels
                    pix_loc = length(find(HT3(Imsk == 1) < 1));
                    pix_red = length(find(HR3(Imsk == 1) < 1));    af_red = pix_red/pix_loc;    af_red_tot = pix_red/pix_tot;
                    pix_orn = length(find(HO3(Imsk == 1) < 1));    af_orn = pix_orn/pix_loc;    af_orn_tot = pix_orn/pix_tot;
                    pix_yel = length(find(HY3(Imsk == 1) < 1));    af_yel = pix_yel/pix_loc;    af_yel_tot = pix_yel/pix_tot;
                    pix_grn = length(find(HG3(Imsk == 1) < 1));    af_grn = pix_grn/pix_loc;    af_grn_tot = pix_grn/pix_tot;
                    
                    % Extract HSL of positive pixels
                    Htemp = H1(Imsk == 1 & HR3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_red = mean(Htemp) + 360;      % Red
                    Htemp = H1(Imsk == 1 & HO3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_orn = mean(Htemp) + 360;      % Orange
                    Htemp = H1(Imsk == 1 & HY3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_yel = mean(Htemp) + 360;      % Yellow
                    Htemp = H1(Imsk == 1 & HG3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_grn = mean(Htemp) + 360;      % Green
                    
                    Savg_red = mean(H2(Imsk == 1 & HR3 < 1));    Lavg_red = mean(H3(Imsk == 1 & HR3 < 1));      % Red
                    Savg_orn = mean(H2(Imsk == 1 & HO3 < 1));    Lavg_orn = mean(H3(Imsk == 1 & HO3 < 1));      % Orange
                    Savg_yel = mean(H2(Imsk == 1 & HY3 < 1));    Lavg_yel = mean(H3(Imsk == 1 & HY3 < 1));      % Yellow
                    Savg_grn = mean(H2(Imsk == 1 & HG3 < 1));    Lavg_grn = mean(H3(Imsk == 1 & HG3 < 1));      % Green
                    
                    afrac{II,JJ}  = cat(2,pix_red,pix_orn,pix_yel,pix_grn,pix_loc,af_red,af_orn,af_yel,af_grn,af_red_tot,af_orn_tot,af_yel_tot,af_grn_tot);
                    hslavg{II,JJ} = cat(1,cat(2,Havg_red,Havg_orn,Havg_yel,Havg_grn),cat(2,Savg_red,Savg_orn,Savg_yel,Savg_grn),cat(2,Lavg_red,Lavg_orn,Lavg_yel,Lavg_grn));
                    
                end
            end
            close(h)
            
            % Compile images for saving
            imstack = cat(4,R,O,Y,G,T);
            
            
            % ----- Movat's Pentachrome (MOV) -----------------------------
        elseif strcmpi(stain,'MOV');
            
            % Load images
            if layer && ~isempty(LDX)
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                G = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_gag-',layernm{LDX,:},'.tif'));           HG = colorspace('RGB->HSL',G);     HG3 = HG(:,:,3);
                F = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_fibrin-',layernm{LDX,:},'.tif'));        HF = colorspace('RGB->HSL',F);     HF3 = HF(:,:,3);
                E = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_elastin-',layernm{LDX,:},'.tif'));       HE = colorspace('RGB->HSL',E);     HE3 = HE(:,:,3);
                M = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_cytoplasm-',layernm{LDX,:},'.tif'));     HM = colorspace('RGB->HSL',M);     HM3 = HM(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_collagen-',layernm{LDX,:},'.tif'));      HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_total-',layernm{LDX,:},'.tif'));         HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            else
                H1 = H(:,:,1);   H2 = H(:,:,2);   H3 = H(:,:,3);
                G = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_gag.tif'));           HG = colorspace('RGB->HSL',G);     HG3 = HG(:,:,3);
                F = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_fibrin.tif'));        HF = colorspace('RGB->HSL',F);     HF3 = HF(:,:,3);
                E = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_elastin.tif'));       HE = colorspace('RGB->HSL',E);     HE3 = HE(:,:,3);
                M = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_cytoplasm.tif'));     HM = colorspace('RGB->HSL',M);     HM3 = HM(:,:,3);
                C = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_collagen.tif'));      HC = colorspace('RGB->HSL',C);     HC3 = HC(:,:,3);
                T = imread(strcat(path,groupnm,'\',fname,'\MOV_',fname,'_total.tif'));         HT = colorspace('RGB->HSL',T);     HT3 = HT(:,:,3);
            end
            
            % Count total pixels
            pix_tot = length(find(HT3 < 1));
            
            % Initialize local area fraction / hsl storage
            npartq = size(localpart,1);     npartr = size(localpart,2);      [xsz,ysz,~] = size(HT(:,:,3));
            afrac = cell(npartq,npartr);    hslavg = cell(npartq,npartr);
            
            h = waitbar(0,'Calculating Local Area Fractions...','Name',char(strcat(stain,{' '},'/',{' '},outnm{IDX,:})));
            
            % Calculate local area fractions
            for II = 1:npartq;
                for JJ = 1:npartr;
                    
                    waitbar(((II-1)*npartr + JJ)/(npartq*npartr))
                    
                    % Generate mask from current partition
                    Imsk = false(xsz,ysz);
                    for KK = 1:length(localpart{II,JJ}(:,1));  Imsk(localpart{II,JJ}(KK,1),localpart{II,JJ}(KK,2)) = true;  end
                    Imsk = imfill(Imsk,'holes');
                    
                    % Extract all non-background pixels
                    pix_loc = length(find(HT3(Imsk == 1) < 1));
                    pix_gag = length(find(HG3(Imsk == 1) < 1));    af_gag = pix_gag/pix_loc;    af_gag_tot = pix_gag/pix_tot;
                    pix_fib = length(find(HF3(Imsk == 1) < 1));    af_fib = pix_fib/pix_loc;    af_fib_tot = pix_fib/pix_tot;
                    pix_eln = length(find(HE3(Imsk == 1) < 1));    af_eln = pix_eln/pix_loc;    af_eln_tot = pix_eln/pix_tot;
                    pix_cyt = length(find(HM3(Imsk == 1) < 1));    af_cyt = pix_cyt/pix_loc;    af_cyt_tot = pix_cyt/pix_tot;
                    pix_col = length(find(HC3(Imsk == 1) < 1));    af_col = pix_col/pix_loc;    af_col_tot = pix_col/pix_tot;
                    
                    % Extract HSL of positive pixels
                    Htemp = H1(Imsk == 1 & HG3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_gag = mean(Htemp) + 360;      % Ground substance (GAG)
                    Htemp = H1(Imsk == 1 & HF3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_fib = mean(Htemp) + 360;      % Fibrin
                    Htemp = H1(Imsk == 1 & HE3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_eln = mean(Htemp) + 360;      % Elastin
                    Htemp = H1(Imsk == 1 & HM3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_cyt = mean(Htemp) + 360;      % Cytoplasm
                    Htemp = H1(Imsk == 1 & HC3 < 1);    Htemp(Htemp > 180) = Htemp(Htemp > 180) - 360;    Havg_col = mean(Htemp) + 360;      % Collagen
                    
                    Savg_gag = mean(H2(Imsk == 1 & HG3 < 1));    Lavg_gag = mean(H3(Imsk == 1 & HG3 < 1));      % Ground substance (GAG)
                    Savg_fib = mean(H2(Imsk == 1 & HF3 < 1));    Lavg_fib = mean(H3(Imsk == 1 & HF3 < 1));      % Fibrin
                    Savg_eln = mean(H2(Imsk == 1 & HE3 < 1));    Lavg_eln = mean(H3(Imsk == 1 & HE3 < 1));      % Elastin
                    Savg_cyt = mean(H2(Imsk == 1 & HM3 < 1));    Lavg_cyt = mean(H3(Imsk == 1 & HM3 < 1));      % Cytoplasm
                    Savg_col = mean(H2(Imsk == 1 & HC3 < 1));    Lavg_col = mean(H3(Imsk == 1 & HC3 < 1));      % Collagen
                    
                    afrac{II,JJ}  = cat(2,pix_gag,pix_fib,pix_eln,pix_cyt,pix_col,pix_tot,af_gag,af_fib,af_eln,af_cyt,af_col,af_gag_tot,af_fib_tot,af_eln_tot,af_cyt_tot,af_col_tot);
                    hslavg{II,JJ} = cat(1,cat(2,Havg_gag,Havg_fib,Havg_eln,Havg_cyt,Havg_col),cat(2,Savg_gag,Savg_fib,Savg_eln,Savg_cyt,Savg_col),cat(2,Lavg_gag,Lavg_fib,Lavg_eln,Lavg_cyt,Lavg_col));
                    
                end
            end
            close(h)
            
            % Compile images for saving
            imstack = cat(4,G,F,E,M,C,T);
            
        end
        
    end
    
end


% ===== Save area fraction/HSL data =======================================

% Define stain-dependent file headings
if strcmpi(stain,'POL')
    cnm = {'Polymer'};
    
elseif strcmpi(stain,'IHC')
    cnm = {'IHC'};
    
elseif strcmpi(stain,'ORO')
    cnm = {'Lipid'};

elseif strcmpi(stain,'VVG')
    cnm = {'Elastin'; 'Tissue'};

elseif strcmpi(stain,'MTC')
    cnm = {'Cytoplasm'; 'Collagen'};
    
elseif strcmpi(stain,'bPSR')
    cnm = {'Collagen'; 'Tissue'};

elseif strcmpi(stain,'dPSR')
    cnm = {'Red'; 'Orange'; 'Yellow'; 'Green'};
    
elseif strcmpi(stain,'VnK') || strcmpi(stain,'AzR')
    cnm = {'Calcification'; 'Tissue'};

elseif strcmpi(stain,'MOV')
    cnm = {'Gag'; 'Fibrin'; 'Elastin'; 'Cytoplasm'; 'Collagen'};
    
elseif strcmpi(stain,'IF')
    %     'IFAb'; 'Nuclear IFAb'; 'Matrix IFAb'
    %     'Positive Nuclei','Positive Nuclei Fraction','Overlap Threshold'
end

% Build array of current output file headings
if strcmpi(stain,'POL');
    
    afcon  = cell(length(cnm),1);
    hslcon = cell(3*length(cnm),1);
    
    for KK = 1:length(cnm)
        afcon{KK,1} = strcat(cnm{KK},' Pix');
        
        hslcon{3*KK-2+0,1} = strcat(cnm{KK},' (Havg)');
        hslcon{3*KK-2+1,1} = strcat(cnm{KK},' (Savg)');
        hslcon{3*KK-2+2,1} = strcat(cnm{KK},' (Lavg)');
    end
    
elseif strcmpi(stain,'IF');
    
    % %     ADD CODE HERE!!
    
else
    afcon_per = cell(length(cnm),1);
    afcon_pix = cell(length(cnm),1);
    hslcon    = cell(3*length(cnm),1);
    
    for KK = 1:length(cnm)
        afcon_per{KK,1} = strcat(cnm{KK},' AF');
        afcon_pix{KK,1} = strcat(cnm{KK},' Pix');
        
        hslcon{3*KK-2+0,1} = strcat(cnm{KK},' (Havg)');
        hslcon{3*KK-2+1,1} = strcat(cnm{KK},' (Savg)');
        hslcon{3*KK-2+2,1} = strcat(cnm{KK},' (Lavg)');
    end
    
    afcon = cat(1,afcon_pix,{'Total Pix'},afcon_per);
end

% ----- Write area fraction/HSL data files --------------------------------
if strcmpi(stain,'IF')
    
    % %     ADD CODE HERE!!
    
else
    af = afrac;    hsl = hslavg;
    for II = 1:length(fileform);
        
        if strcmpi(fileform{II},'.txt')
            
            afcontxt = afcon;    hslcontxt = hslcon;
                        
            % Setup correct filepaths based on type of analysis
            if indiv && (~layer || (layer && isempty(LDX))) && ~iscell(af)
                addpart = false;    addname = false;    addtot = false;    type = 1;
                fid1 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_HSL.txt'),'w');
                
            elseif indiv && (layer && ~isempty(LDX)) && ~iscell(af)
                addpart = false;    addname = false;    addtot = false;    type = 1;
                fid1 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF_',layernm{LDX,:},'.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_HSL_',layernm{LDX,:},'.txt'),'w');
                
            elseif indiv && (~layer || (layer && isempty(LDX))) && iscell(af)
                addpart = true;    addname = false;    addtot = true;    type = 2;
                fid1 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF_',char(vartype),'.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_HSL_',char(vartype),'.txt'),'w');
                
            elseif indiv && (layer && ~isempty(LDX)) && iscell(af)
                addpart = true;    addname = false;    addtot = true;    type = 2;
                fid1 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF_',layernm{LDX,:},'-',char(vartype),'.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_HSL_',layernm{LDX,:},'-',char(vartype),'.txt'),'w');
                
            elseif ~indiv && (~layer || (layer && isempty(LDX))) && ~iscell(af)
                addpart = false;    addname = true;    addtot = false;    type = 3;
                fid1 = fopen(strcat(path,groupnm,'\All_',groupnm,'_AF.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\All_',groupnm,'_HSL.txt'),'w');
                
            elseif ~indiv && (layer && ~isempty(LDX)) && ~iscell(af)
                addpart = false;    addname = true;    addtot = false;    type = 3;
                fid1 = fopen(strcat(path,groupnm,'\All_',groupnm,'_AF_',layernm{LDX,:},'.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\All_',groupnm,'_HSL_',layernm{LDX,:},'.txt'),'w');
                
            elseif ~indiv && (~layer || (layer && isempty(LDX))) && iscell(af)
                addpart = true;    addname = true;    addtot = true;    type = 4;
                fid1 = fopen(strcat(path,groupnm,'\All_',groupnm,'_AF_',char(vartype),'.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\All_',groupnm,'_HSL_',char(vartype),'.txt'),'w');
                
            elseif ~indiv && (layer && ~isempty(LDX)) && iscell(af)
                addpart = true;    addname = true;    addtot = true;    type = 4;
                fid1 = fopen(strcat(path,groupnm,'\All_',groupnm,'_AF_',layernm{LDX,:},'-',char(vartype),'.txt'),'w');
                fid2 = fopen(strcat(path,groupnm,'\All_',groupnm,'_HSL_',layernm{LDX,:},'-',char(vartype),'.txt'),'w');
            end
            
            % Build formatted strings for text headers
            if addpart
                aftxt1 = '''%s\t%s\t%s\t';                                      hsltxt1 = aftxt1;
                aftxt2 = '''Image Number'',''Circ. Part.'',''Rad. Part.''';     hsltxt2 = aftxt2;
                aftxt3 = '''%2.0f\t%2.0f\t%2.0f\t';                             hsltxt3 = aftxt3;
            else
                aftxt1 = '''%s\t';                                              hsltxt1 = aftxt1;
                aftxt2 = '''Image Number''';                                    hsltxt2 = aftxt2;
                aftxt3 = '''%2.0f\t';                                           hsltxt3 = aftxt3;
            end
            
            % Build formatted strings for data storage (area fraction)
            for KK = 1:length(afcontxt);
                
                aftxt1 = strcat(aftxt1,'%s\t');
                aftxt2 = strcat(aftxt2,',''',afcontxt{KK},'''');
                
                if KK <= ceil(length(afcontxt)/2);
                    aftxt3 = strcat(aftxt3,'%8.0f\t');
                else
                    aftxt3 = strcat(aftxt3,'%7.4f\t');
                end
            end
            
            % Build formatted strings for data storage (HSL)
            for KK = 1:length(hslcontxt);
                hsltxt1 = strcat(hsltxt1,'%s\t');
                hsltxt2 = strcat(hsltxt2,',''',hslcontxt{KK},'''');
                hsltxt3 = strcat(hsltxt3,'%7.4f\t');
            end
            
            % Adjust strings for partition outputs
            if addtot
                for KK = 1:length(cnm)
                    aftxt1 = strcat(aftxt1,'%s\t');
                    aftxt2 = strcat(aftxt2,',''',strcat(cnm{KK},' AF Tot.'),'''');
                    aftxt3 = strcat(aftxt3,'%9.6f\t');
                end
            end
            
            % Adjust strings for group outputs
            if addname
                aftxt1 = strcat(aftxt1,'%s\t\n''');             hsltxt1 = strcat(hsltxt1,'%s\t\n''');
                aftxt2 = strcat(aftxt2,',''Image Name''');      hsltxt2 = strcat(hsltxt2,',''Image Name''');
                aftxt3 = strcat(aftxt3,'%s\t\n''');             hsltxt3 = strcat(hsltxt3,'%s\t\n''');
            else
                aftxt1 = strcat(aftxt1,'\n''');                 hsltxt1 = strcat(hsltxt1,'\n''');
                aftxt3 = strcat(aftxt3,'\n''');                 hsltxt3 = strcat(hsltxt3,'\n''');
            end
            
            % Write area fraction / hsl data to file
            eval(strcat('fprintf(fid1,',aftxt1,',',aftxt2,');'))
            eval(strcat('fprintf(fid2,',hsltxt1,',',hsltxt2,');'))
            
            if type == 1;        % Individual, no partition
                
                hsl = hsl(:)';
                eval(strcat('fprintf(fid1,',aftxt3,',1,af);'));
                eval(strcat('fprintf(fid2,',hsltxt3,',1,hsl);'));
                
            elseif type == 2;    % Individual, with partition
                
                npartq = size(af,1);    npartr = size(af,2);
                
                for JJ = 1:npartq;
                    for KK = 1:npartr;
                        hslz = hsl{JJ,KK}(:)';  %#ok
                        eval(strcat('fprintf(fid1,',aftxt3,',1,JJ,KK,af{JJ,KK});'));
                        eval(strcat('fprintf(fid2,',hsltxt3,',1,JJ,KK,hslz);'));
                    end
                end
                
            elseif type == 3;    % Group, no partition
                fprintf('WRITE CODE FOR .TXT SAVING!!')
            elseif type == 4;    % Group, with partition
                fprintf('WRITE CODE FOR .TXT SAVING!!')
            end
            
            fclose(fid1);    fclose(fid2);
            
        elseif strcmpi(fileform{II},'.xls')
            
            afconxls = afcon;    hslconxls = hslcon;
            
            % Setup file paths based on type of analysis
            if indiv && (~layer || (layer && isempty(LDX))) && ~iscell(af)
                addpart = false;    addname = false;    addtot = false;    type = 1;
                fpath = strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'.xls');
                
            elseif indiv && (layer && ~isempty(LDX)) && ~iscell(af)
                addpart = false;    addname = false;    addtot = false;    type = 1;
                fpath = strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_',layernm{LDX,:},'.xls');
                
            elseif indiv && (~layer || (layer && isempty(LDX))) && iscell(af)
                addpart = true;    addname = false;    addtot = true;    type = 2;
                fpath = strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_',char(vartype),'.xls');
                
            elseif indiv && (layer && ~isempty(LDX)) && iscell(af)
                addpart = true;    addname = false;    addtot = true;    type = 2;
                fpath = strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_',layernm{LDX,:},'-',char(vartype),'.xls');
                
            elseif ~indiv && (~layer || (layer && isempty(LDX))) && ~iscell(af)
                addpart = false;    addname = true;    addtot = false;    type = 3;
                fpath = strcat(path,groupnm,'\All_',groupnm,'.xls');
                
            elseif ~indiv && (layer && ~isempty(LDX)) && ~iscell(af)
                addpart = false;    addname = true;    addtot = false;    type = 3;
                fpath = strcat(path,groupnm,'\All_',groupnm,'_',layernm{LDX,:},'.xls');
                
            elseif ~indiv && (~layer || (layer && isempty(LDX))) && iscell(af)
                addpart = true;    addname = true;    addtot = true;    type = 4;
                fpath = strcat(path,groupnm,'\All_',groupnm,'_',char(vartype),'.xls');
                
            elseif ~indiv && (layer && ~isempty(LDX)) && iscell(af)
                addpart = true;    addname = true;    addtot = true;    type = 4;
                fpath = strcat(path,groupnm,'\All_',groupnm,'_',layernm{LDX,:},'-',char(vartype),'.xls');
            end
            
            % Build formatted strings for text headers
            if addpart
                afconxls  = cat(1,{'Image Number'; 'Circ. Part.'; 'Rad. Part.'},afconxls);
                hslconxls = cat(1,{'Image Number'; 'Circ. Part.'; 'Rad. Part.'},hslconxls);
            else
                afconxls  = cat(1,{'Image Number'},afconxls);
                hslconxls = cat(1,{'Image Number'},hslconxls);
            end
            
            % Adjust headers for partition outputs
            if addtot
                for KK = 1:length(cnm)
                    afconxls = cat(1,afconxls,{strcat(cnm{KK},' AF Tot.')});
                end
            end
            
            % Adjust headers for group outputs
            if addname
                afconxls = cat(1,afconxls,{'Image Name'});      hslconxls = cat(1,hslconxls,{'Image Name'});
            end
            
            % Define cell widths for each column
            cwidth_af = 20*ones(1,length(afconxls));    cwidth_hsl = 20*ones(1,length(hslconxls));
            
            cwidth_af(1) = 15;    cwidth_hsl(1) = 15;
            if addpart;    cwidth_af(2:3) = 15;    cwidth_hsl(2:3) = 15;    end
            if addname;    cwidth_af(end) = 30;    cwidth_hsl(end) = 30;    end
            
            % Write data to excel file
            xlsformat(fpath, afconxls, hslconxls, cwidth_af, cwidth_hsl, af, hsl, indiv)
            
        elseif strcmpi(fileform{II},'.mat')
            
            % Setup correct filepaths based on type of analysis
            if indiv && (~layer || (layer && isempty(LDX))) && ~iscell(af)
                save(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF-HSL.mat'),'af','hsl');
                
            elseif indiv && (layer && ~isempty(LDX)) && ~iscell(af)
                save(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF-HSL_',layernm{LDX,:},'.mat'),'af','hsl');
                
            elseif indiv && (~layer || (layer && isempty(LDX))) && iscell(af)
                save(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF-HSL_',char(vartype),'.mat'),'af','hsl');
                
            elseif indiv && (layer && ~isempty(LDX)) && iscell(af)
                save(strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_AF-HSL_',layernm{LDX,:},'-',char(vartype),'.mat'),'af','hsl');
                
            elseif ~indiv && (~layer || (layer && isempty(LDX))) && ~iscell(af)
                save(strcat(path,groupnm,'\All_',groupnm,'_AF-HSL.mat'),'af','hsl');
                
            elseif ~indiv && (layer && ~isempty(LDX)) && ~iscell(af)
                save(strcat(path,groupnm,'\All_',groupnm,'_AF-HSL_',layernm{LDX,:},'.mat'),'af','hsl');
                
            elseif ~indiv && (~layer || (layer && isempty(LDX))) && iscell(af)
                save(strcat(path,groupnm,'\All_',groupnm,'_AF-HSL_',char(vartype),'.mat'),'af','hsl');
                
            elseif ~indiv && (layer && ~isempty(LDX)) && iscell(af)
                save(strcat(path,groupnm,'\All_',groupnm,'_AF-HSL_',layernm{LDX,:},'-',char(vartype),'.mat'),'af','hsl');
            end
            
        end
        
    end
    
end


% ===== Save area fraction/HSL data =======================================
if imsave && indiv && isempty(localpart);
    
    % Save individual constituent images
    cnm = cat(1,cnm,{'Total'});
    
    for II = 1:length(cnm);
        if (layer && ~isempty(LDX))
            imwrite(imstack(:,:,:,II),strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_',strcat(lower(cnm{II,1}(1)),cnm{II,1}(2:end)),'-',layernm{LDX,:},'.tif'),'tif');
        else
            imwrite(imstack(:,:,:,II),strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_',strcat(lower(cnm{II,1}(1)),cnm{II,1}(2:end)),'.tif'),'tif');
        end
    end

    % Build pseudocolored image with each constituent
    Ip = imstack(:,:,:,1);    Ip1 = Ip(:,:,1);    Ip2 = Ip(:,:,2);       Ip3 = Ip(:,:,3);
    
    for KK = 2:(size(imstack,4)-1);

        Imsk = ~((imstack(:,:,1,KK) == 1) & (imstack(:,:,2,KK) == 1) & (imstack(:,:,3,KK) == 1));
        
        I1 = imstack(:,:,1,KK);    Ip1(Imsk == 1) = I1(Imsk == 1);  
        I2 = imstack(:,:,2,KK);    Ip2(Imsk == 1) = I2(Imsk == 1);             
        I3 = imstack(:,:,3,KK);    Ip3(Imsk == 1) = I3(Imsk == 1); 
    end
    
    Ip = cat(3,Ip1,Ip2,Ip3);
    
    % Save pseudocolored image
    if (layer && ~isempty(LDX))
        imwrite(Ip,strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_combined-',layernm{LDX,:},'.tif'),'tif');
    else
        imwrite(Ip,strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_combined.tif'),'tif');
    end
       
elseif imsave && indiv && ~isempty(localpart);
    
    % Plot overlay of partitions on segmented image
    for II = 1:length(cnm);
        
        figure; imshow(imstack(:,:,:,II)); hold on
        
        for JJ = 1:npartq;
            for KK = 1:npartr;
                plot(localpart{JJ,KK}(:,2),localpart{JJ,KK}(:,1),'LineWidth',3,'Color',[0.0,0.0,0.0]);
            end
        end
        
        % Extract color data from figure and write partitioned image to directory
        set(gcf,'color','w');    f = getframe(gcf);
        
        if (layer && ~isempty(LDX))
            imwrite(f.cdata,strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_',strcat(lower(cnm{II,1}(1)),cnm{II,1}(2:end)),'-',layernm{LDX,:},'-partition-',char(vartype),'.tif'),'tif');
        else
            imwrite(f.cdata,strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_',strcat(lower(cnm{II,1}(1)),cnm{II,1}(2:end)),'-partition-',char(vartype),'.tif'),'tif');
        end
        
        close all force
        
    end
        
    % Build pseudocolored image
    Ip = imstack(:,:,:,1);    Ip1 = Ip(:,:,1);    Ip2 = Ip(:,:,2);       Ip3 = Ip(:,:,3);
    
    for KK = 2:(size(imstack,4)-1);

        Imsk = ~((imstack(:,:,1,KK) == 255) & (imstack(:,:,2,KK) == 255) & (imstack(:,:,3,KK) == 255));
        
        I1 = imstack(:,:,1,KK);    Ip1(Imsk == 1) = I1(Imsk == 1);  
        I2 = imstack(:,:,2,KK);    Ip2(Imsk == 1) = I2(Imsk == 1);             
        I3 = imstack(:,:,3,KK);    Ip3(Imsk == 1) = I3(Imsk == 1); 
    end
    
    Ip = cat(3,Ip1,Ip2,Ip3);
    
    % Plot overlay of partitions on segmented image
    figure; imshow(Ip); hold on
    
    for JJ = 1:npartq;
        for KK = 1:npartr;
            plot(localpart{JJ,KK}(:,2),localpart{JJ,KK}(:,1),'LineWidth',3,'Color',[0.0,0.0,0.0]);
        end
    end
    
    % Extract color data from figure and write partitioned image to directory
    set(gcf,'color','w');    f = getframe(gcf);
    
    if (layer && ~isempty(LDX))
        imwrite(f.cdata,strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_combined-',layernm{LDX,:},'-partition-',char(vartype),'.tif'),'tif');
    else
        imwrite(f.cdata,strcat(path,groupnm,'\',fname,'\',stain,'_',fname,'_combined-partition-',char(vartype),'.tif'),'tif');
    end
    
    close all force
        
end



% _________________________ nested functions ______________________________
% -------------------------------------------------------------------------
    function xlsformat(fpath, header_af, header_hsl, cwidth_af, cwidth_hsl, af, hsl, indiv)   % Write formatted text file
        
        % Arrange area fraction / HSL data
        if indiv && ~iscell(af)
            afdata = num2cell([1 af]);    hsldata = num2cell([1 hsl(:)']);
            
        elseif indiv && iscell(af)
            
            npq = size(af,1);    npr = size(af,2);    afdata = [];    hsldata = [];
            for ii = 1:npq;
                for jj = 1:npr;
                    afdata  = cat(1,afdata,num2cell([1 ii jj af{ii,jj}]));
                    hsltemp = hsl{ii,jj};    hsldata = cat(1,hsldata,num2cell([1 ii jj hsltemp(:)']));
                end
            end
            
        elseif ~indiv && ~iscell(af)
            
            im = size(af,1);
            for ii = 1:im;
                afdata = num2cell([ii af(ii,:)]);    hsldata = num2cell([ii hsl(ii,:)]);
            end
            
        elseif ~indiv && iscell(af)
            
            im = size(af,1);    npq = size(af{1,1},1);    npr = size(af{1,1},2);
            for ii = 1:im;
                for jj = 1:npq;
                    for kk = 1:npr;
                        afdata = num2cell([ii jj kk af{ii}{jj,kk}]);    hsldata = num2cell([ii jj kk hsl{ii}{jj,kk}]);
                    end
                end
            end
            
        end
        
        % Compile into single array
        afdata_write  = [header_af'; afdata];      Na = size(afdata_write,1);
        hsldata_write = [header_hsl'; hsldata];    Nh = size(hsldata_write,1);
        
        % Write data to file
        if exist(fpath,'file') ~= 0
            delete(fpath);    xlswrite(fpath,afdata_write,'Sheet1');    xlswrite(fpath,hsldata_write,'Sheet2');
        else
            xlswrite(fpath,afdata_write,'Sheet1');    xlswrite(fpath,hsldata_write,'Sheet2');
        end
        
        % Initiate connection with Excel
        hExcel = actxserver('Excel.Application');
        
        % Open workbook and define sheet names
        hWorkbook = hExcel.Workbooks.Open(fpath);
        hWorksheetAF  = hWorkbook.Sheets.Item('Sheet1');    hWorksheetAF.Name  = 'Area Fraction';
        hWorksheetHSL = hWorkbook.Sheets.Item('Sheet2');    hWorksheetHSL.Name = 'HSL';
        
        % Set column widths
        for ii = 1:length(cwidth_af);      hWorksheetAF.Columns.Item(ii).columnWidth = cwidth_af(ii);     end
        for ii = 1:length(cwidth_hsl);    hWorksheetHSL.Columns.Item(ii).columnWidth = cwidth_hsl(ii);    end
        
        % Set column justification
        for ii = 1:length(cwidth_af);      hWorksheetAF.Columns.Item(ii).HorizontalAlignment = -4108;     end
        for ii = 1:length(cwidth_hsl);    hWorksheetHSL.Columns.Item(ii).HorizontalAlignment = -4108;    end
        
        % Set cell formatting based on type of analysis - AF
        for ii = 1:length(header_af);
            
            cells = get(hWorksheetAF.cells,'item',1,ii);
            hB = cells.borders;    hBL = get(hB,'item',4);    set(hBL,'linestyle',3);
            
            if ~iscell(af)
                
                % First column
                if ii == 1;
                    for jj = 1:Na;
                        cells = get(hWorksheetAF.cells,'item',jj,ii);
                        hB = cells.borders;    hBR = get(hB,'item',2);      set(hBR,'linestyle',3);
                    end
                end
                
                % Formatting column index
                fidx = floor(1+((((length(header_af)-1)-1)/2)+1)); 
                
                % Pix-to-AF column and Last column
                if indiv;    chk = (ii == fidx+1);    else    chk = ((ii == fidx+1) || (ii == length(header_af)));    end
                
                if chk
                    for jj = 1:Na;
                        cells = get(hWorksheetAF.cells,'item',jj,ii);
                        hB = cells.borders;    hBL = get(hB,'item',1);      set(hBL,'linestyle',3);
                    end
                end
                
                % AF formatting
                if indiv;     chk = (ii > fidx);     else     chk = ((ii > fidx) && (ii < length(header_af)));    end
                
                if chk;
                    for jj = 2:Na;
                        cells = get(hWorksheetAF.cells,'item',jj,ii);    set(cells,'NumberFormat','0.0000');
                    end
                end
                
                
            elseif iscell(af)
                
                % First column
                if ii == 1 || ii == 3;
                    for jj = 1:Na;
                        cells = get(hWorksheetAF.cells,'item',jj,ii);
                        hB = cells.borders;    hBR = get(hB,'item',2);      set(hBR,'linestyle',3);
                    end
                end
                
                % Formatting column index
                fidx = floor(3+((((length(header_af)-3)-1)/3)+1));
                
                % Pix-to-AF column and Last column
                if indiv;    chk = ((ii == fidx+1) || (ii == fidx+1+((length(header_af)-3)-1)/3));    else    chk = ((ii == fidx+1) || (ii == fidx+1+((length(header_af)-3)-1)/3) || (ii == length(header_af)));    end
                
                if chk
                    for jj = 1:Na;
                        cells = get(hWorksheetAF.cells,'item',jj,ii);
                        hB = cells.borders;    hBL = get(hB,'item',1);      set(hBL,'linestyle',3);
                    end
                end
                
                % AF formatting
                if indiv;     chk1 = (ii > fidx);    chk2 = (ii > fidx+((length(header_af)-3)-1)/3);     else     chk1 = ((ii > fidx) && (ii < length(header_af)));    chk2 = ((ii > fidx+((length(header_af)-3)-1)/3) && (ii < length(header_af)));    end
                
                if chk1;
                    for jj = 2:Na;    cells = get(hWorksheetAF.cells,'item',jj,ii);    set(cells,'NumberFormat','0.0000');      end
                elseif chk2;
                    for jj = 2:Na;    cells = get(hWorksheetAF.cells,'item',jj,ii);    set(cells,'NumberFormat','0.000000');    end
                end

            end
        end
        
        % Set cell formatting based on type of analysis - HSL
        for ii = 1:length(header_hsl);
            
            cells = get(hWorksheetHSL.cells,'item',1,ii);
            hB = cells.borders;    hBL = get(hB,'item',4);    set(hBL,'linestyle',3);
            
            if ~iscell(af)
                
                % First column
                if ii == 1;
                    for jj = 1:Nh;
                        cells = get(hWorksheetHSL.cells,'item',jj,ii);
                        hB = cells.borders;    hBR = get(hB,'item',2);      set(hBR,'linestyle',3);
                    end
                end
                
                % Formatting column index
                nidx = (length(header_hsl)-1)/3;
                for jj = 1:nidx;    fidx(jj,1) = 3*jj + 1;    end
                
                % Pix-to-AF column and Last column
                if indiv;    chk = any(ii == fidx(1:end-1)+1);    else    chk = (any(ii == fidx(1:end-1)+1) || ii == length(header_hsl));    end
                
                if chk
                    for jj = 1:Nh;
                        cells = get(hWorksheetHSL.cells,'item',jj,ii);
                        hB = cells.borders;    hBL = get(hB,'item',1);      set(hBL,'linestyle',3);
                    end
                end
                
                % HSL formatting
                if indiv;       chk = (ii > 1);      else      chk = ((ii > 1) && (ii < length(header_hsl)));      end
                
                if chk;
                    for jj = 2:Nh;
                        cells = get(hWorksheetHSL.cells,'item',jj,ii);    set(cells,'NumberFormat','0.0000');
                    end
                end
                
            elseif iscell(af)
                
                % First column
                if ii == 1 || ii == 3;
                    for jj = 1:Nh;
                        cells = get(hWorksheetHSL.cells,'item',jj,ii);
                        hB = cells.borders;    hBR = get(hB,'item',2);      set(hBR,'linestyle',3);
                    end
                end
                
                % Formatting column index
                nidx = (length(header_hsl)-3)/3;
                for jj = 1:nidx;    fidx(jj,1) = 3*jj + 3;    end
                
                % Pix-to-AF column and Last column
                if indiv;    chk = any(ii == fidx(1:end-1)+1);    else    chk = (any(ii == fidx(1:end-1)+1) || ii == length(header_hsl));    end
                
                if chk
                    for jj = 1:Nh;
                        cells = get(hWorksheetHSL.cells,'item',jj,ii);
                        hB = cells.borders;    hBL = get(hB,'item',1);      set(hBL,'linestyle',3);
                    end
                end
                
                % HSL formatting
                if indiv;       chk = (ii > 3);      else      chk = ((ii > 3) && (ii < length(header_hsl)));      end
                
                if chk;
                    for jj = 2:Nh;
                        cells = get(hWorksheetHSL.cells,'item',jj,ii);    set(cells,'NumberFormat','0.0000');
                    end
                end
                
            end
        end
        
        % Save and close Excel file
        hWorkbook.Save;    hWorkbook.Close;    hExcel.Quit;
        
    end


    function [Havg, Savg, Lavg] = HSLavg(A, H)
        % Average H,S,and L values for current region in image
        
        % Isolate constituent pixels (and their location in the image) to compute average H-S-L values
        val = A(:,:,3)-1;    Lind = find(val);
        
        % Separate HSL image into individual channels
        Hval = H(:,:,1);    Sval = H(:,:,2);    Lval = H(:,:,3);
        
        % Initialize storage and extract constituent specific pixels from each H-S-L channel
        Hpt = zeros(length(Lind(:,1)),1);    Spt = zeros(length(Lind(:,1)),1);    Lpt = zeros(length(Lind(:,1)),1);
        
        % Constituent specific Hue, Saturation, and Lightness matrices
        for i = 1:length(Lind)
            Hpt(i,1) = Hval(Lind(i,1));    Spt(i,1) = Sval(Lind(i,1));    Lpt(i,1) = Lval(Lind(i,1));
        end
        
        % Periodically offset H values...helps with RED averages
        Hpt(Hpt > 180) = Hpt(Hpt > 180) - 360;
        
        % Compute average H-S-L value for constituent
        Havg = mean(Hpt) + 360;    Savg = mean(Spt);    Lavg = mean(Lpt);
        
    end

    function choice = genquestdlg(dlgOptions, dlgTitle, defOption, qStr, bttnsOredring)
        % genquestdlg
        % Create and open a button dialog box with many buttons.
        
        % Default params
        if nargout > 1;    error('MATLAB:genquestdlg:WrongNumberOutputs','Wrong number of output arguments for QUESTDLG');    end
        if nargin  < 1;    error('MATLAB:genquestdlg:WrongNumberInputs','Wrong number of input arguments for genquestdlg');   end
        if nargin == 1;           dlgTitle = ' ';         end
        if nargin <= 2;     defOption = dlgOptions{1};    end
        if nargin <= 3;    qStr = [];    titleSps = 0;    end
        if nargin <= 4;         bttnsOredring = [];       end
        if nargin >  5;    error('MATLAB:genquestdlg:TooManyInputs', 'Too many input arguments');    end
        
        % internal params
        bttnFontSize = 0.7;    btntxtH = 2;
        
        % Buttons ordering definition
        nButtons = length(dlgOptions);    nLongestOption = max(cellfun(@length, dlgOptions));
        
        % Set buttons ordering- N Columns and N Rows
        if isempty(bttnsOredring)
            bttnsOredring    = zeros(1,2);
            bttnsOredring(1) = ceil(sqrt(nButtons));
            bttnsOredring(2) = ceil(nButtons/bttnsOredring(1));
        end
        
        if bttnsOredring(1) > 1 && bttnsOredring(1) <= nButtons
            bttnRows = bttnsOredring(1);
            bttnCols = ceil(nButtons/bttnRows);
        else
            if bttnsOredring(2) > 1 && bttnsOredring(2) <= nButtons
                bttnCols = bttnsOredring(2);
            else
                bttnCols = floor(sqrt(nButtons));
            end
            bttnRows = ceil(nButtons/bttnCols);
        end
        
        if exist('titleSps','var') ~= 1
            titleSps = 1.25*btntxtH;  % Title gets more space then buttons.
        end
        spaceH = 0.5;    spaceW = 2;
        
        % ----- Dialog Figure definition ----------------------------------
        % Open a figure about screen center
        menuFigH=figure('Units', 'normalized', 'Position', [.5, .5, .1, .1], 'MenuBar', 'none',...
            'NumberTitle', 'off', 'Name', dlgTitle, 'CloseRequestFcn', 'uiresume(gcbf)');
        % 'CloseRequestFcn' override figure closing
        
        % make sure figure form allows good text representation Get screen resolution in characters
        getRootUnits = get(0,'Units');         set(0,'Units','characters');
        ScreenSize   = get(0,'ScreenSize');    set(0,'Units',getRootUnits);
        
        set(menuFigH,'Units','characters');    FigPos = get(menuFigH,'Position');
        figH = titleSps + (btntxtH+spaceH)*bttnRows + spaceH;
        figW = max(spaceW + btntxtH*bttnCols*(nLongestOption+spaceW),titleSps*length(qStr) + 2*spaceW);
        
        if figW > ScreenSize(3);    figW = ScreenSize(3);    end    % raise some flag, figure Width will be too small, use textwrap...
        if figH > ScreenSize(4);    figH = ScreenSize(4);    end    % raise some flag, figure Height will be too small
        
        FigL = FigPos(1) - figW/2;    FigB = FigPos(2) - figH/2;    FigPos = [FigL, FigB, figW, figH];    set(menuFigH, 'Position', FigPos);
        
        % Button group definition
        iDefOption = 1;
        buttonGroup = uibuttongroup('Parent',menuFigH,'Position',[0 0 1 1]);
        if exist('defOption','var') == 1 && ~isempty(defOption)
            if ischar(defOption)
                iDefOption = strcmpi(defOption, dlgOptions);
            elseif isnumeric(defOption) && (defOption > 1) && (defOption < nButtons)
                iDefOption = defOption; % defOption is an IDX of the default option
            end
        end
        
        % Question definition
        if titleSps > 0
            titleH = uicontrol( buttonGroup, 'Style', 'text', 'String', qStr,'FontUnits', 'normalized',  'FontSize', bttnFontSize,...
                'HorizontalAlignment', 'center', 'Units', 'characters', 'Position', [spaceW , figH-titleSps, figW-2*spaceW, titleSps] );
        end
        
        % Radio buttons definition
        bttnH = max(1,(figH-titleSps-spaceH)/bttnRows-spaceH);    interBttnStpH = bttnH + spaceH;
        bttnW = nLongestOption*btntxtH;                           interBttnStpW = max((figW-spaceW)/(bttnCols), bttnW+spaceW);
        
        buttnHndl = zeros(nButtons,1);
        for iBttn = 1:nButtons
            [iRow,iCol] = ind2sub([bttnRows, bttnCols],iBttn);
            currBttnHeigth = figH - titleSps-iRow*interBttnStpH;
            currBttnLeft = (iCol-1)*interBttnStpW + spaceW + (interBttnStpW-bttnW-spaceW)/2;
            
            buttnHndl(iBttn) = uicontrol( buttonGroup, 'Style', 'pushbutton', 'FontUnits', 'normalized', 'Units', 'characters', 'FontSize', bttnFontSize,...
                'String', dlgOptions{iBttn}, 'Callback', @my_bttnCallBack, 'Position', [currBttnLeft , currBttnHeigth, bttnW, bttnH] );
        end
        set(buttnHndl(iDefOption), 'Value', 1); % set default option
        set(cat(1, buttnHndl, titleH, buttonGroup),'Units', 'normalized')
        
        uiwait(menuFigH); % wait untill user makes his choise, or closes figure
        delete(menuFigH);
        
        function my_bttnCallBack(hObject, ~)
            choice = find(strcmpi(dlgOptions,get(hObject,'String')));
            uiresume(gcbf)
        end
        
    end

end