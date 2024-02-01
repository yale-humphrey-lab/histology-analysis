% ========================================================================
%  HISTOLOGY ANALYSIS PACKAGE
%  Matthew Bersi, PhD
%
%  Updated: September 2017
%  Changed: March 2022
% ========================================================================

clc; clearvars; close all force; warning('off','all');

% Names of stains and components that can be quantified in analysis:
% --- HISTOLOGY ANALYSIS
%   - General:
%       VVG - Verhoeff's van Gieson - Elastin (Black)
%       MTC - Masson's Trichrome - SMC/Cytoplasm (Red), Collagen (Blue)
%       PSR - Picrosirius Red - Collagen [Gradation: Thick (Red) to Thin (Green)]
%       MOV - Movat's Pentachrome - Elastin (Black), SMC/Cytoplasm (Deep Red-Purple), Collagen (Yellow-Brown), Ground Substance (Blue-Green), Fibrin (Bright Red)
%       POL - Polymer (Birefringence) - Polymer (Non-black)
%       IHC - Immunohistochemistry
%
%   - Artery / Vasculature:
%       VVG - Verhoeff's van Gieson - Elastin (Black)
%       MTC - Masson's Trichrome - SMC/Cytoplasm (Red), Collagen (Blue)
%       PSR - Picrosirius Red - Collagen [Gradation: Thick (Red) to Thin (Green)]
%       MOV - Movat's Pentachrome - Elastin (Black), SMC/Cytoplasm (Deep Red-Purple), Collagen (Yellow-Brown), Ground Substance (Blue-Green), Fibrin (Bright Red)
%       POL - Polymer (Birefringence) - Polymer (Non-black)
%       IHC - Immunohistochemistry
%
%   - Atherosclerosis:
%       MTC - Masson's Trichrome - SMC/Cytoplasm (Red), Collagen (Blue)
%       ORO - Oil Red-O - Lipid (Red)
%       EnF - Whole-mount Aorta Enface Staining - Lipid (Red)
%
%   - Myocardial Infarction
%       MTC - Masson's Trichrome - SMC/Cytoplasm (Red), Collagen (Blue) // Scar Morphology
%
% --- IMMUNOFLUORESCENCE ANALYSIS

% =========================================================================
% Define global variables
% =========================================================================

global stain tissue proctype % is proctype needed?

global layer radial circum layernm rembleb enhance

global path ext SF

global groupnm outnm rednm mergenm dapinm

global imsave vartype

% =========================================================================
% Select type of analysis to perform
% =========================================================================
% Load histology processing path...
addpath(genpath('histology-analysis'))

% Determine type of analysis to perform
[inputs,cancel] = histodlg([],'selection');

% Check if histology or immunofluorescence...
histoproc  = inputs.histoproc;          if histoproc;       proctype = {'histo'};       tissue = inputs.tissuehist;      end
immunoproc = inputs.immunoproc;         if immunoproc;      proctype = {'immuno'};      tissue = inputs.tissueif;        end

% If needed, load immunofluorescence directory
if immunoproc;
    addpath(genpath('immuno'));
end


% !!!!!!!!!!!!!!!!!!!!!!!!!!
% ADD FILENAME ORGANIZATION, STITCHING, REGISTRATION, ETC. FOR IMMUNOFLUORESCENCE ANALYSIS
% !!!!!!!!!!!!!!!!!!!!!!!!!!

% =========================================================================
% Input parameters for function files
% =========================================================================

% Begin image processing...
while cancel < 1;
    
    % Load correct initialization interface
    if histoproc;
        [inputs,cancel] = histodlg([],'initialize-histo');
    elseif immunoproc;
        [inputs,cancel] = histodlg([],'initialize-immuno');
    end
    
    % Terminate if 'cancel' is selected
    if cancel;      return;      end
    
    % Input/Output preferences for image analysis.
    % -- List of images to process
    name = inputs.name;
    
    % -- List of image properties and processing parameters
    SF = [inputs.SFr inputs.SFm];   % Scale factor
    ext = char(inputs.ext);         % File extension
    path = inputs.path;             % Filepath for current images
    stain = char(inputs.stain);     % Histological stain
    outnm = inputs.outnm;           % Output file names
    layer = inputs.layer;           % Layer segmentation: 0 - no, 1 - yes
    radial = inputs.radvar;         % Radial partition: 0 - no, 1 - yes
    circum = inputs.circvar;        % Circumferential partition: 0 - no, 1 - yes
    npartr = inputs.npr;            % Number of rad. partitions
    npartq = inputs.npq;            % Number of circ. partitions
    nlayer = inputs.nlayer;         % Number of layers
    imsave = inputs.imsave;         % Image save: 0 - no, 1 - yes
    rembleb = 'No';                 % Remove blebs (applies to VVG and MOV images of Vasculature)
    vartype = inputs.vartype;       % Type of processing: 'none', 'rad', 'circ', 'both'
    groupnm = inputs.groupnm;       % Output group name
    imthresh = inputs.imthresh;     % Type of thresholding: 'auto', 'manual'
    fileform = inputs.fileform;     % Format to save output files
    
    % -- Image Enhancement
    if histoproc
        enhance = questdlg('Enhance color of images prior to processing?','Image Enhancement','Yes','No','No');
    end
    
    % -- PSR selection
    if strcmpi(stain,'PSR');
        qpsr = questdlg('Select type of PSR image to process...','PSR Selection','Brightfield','Darkfield','Darkfield');
        if strcmpi(qpsr,'Brightfield');    stain = 'bPSR';    else;    stain = 'dPSR';    end
    end
    
    % -- List of immunofluorescence-specific image parameters
    if immunoproc
        rednm = inputs.rednm;                       % Name of red fluorophore target
        redelastin = inputs.redelastin;             % Red elastin flag: 0 - no, 1 - yes
        mfilename = cell(length(name(:,1)),1);      % Array to store modified filenames
    end
    
    % Save current path for next run...
    save('path.mat','path')
    
    
    % =========================================================================
    % Histological / Area Fraction analysis
    % =========================================================================
    % -- Background subtraction / Global thresholding ---------------------
    proc = 'Yes';    procLV = 'Yes';    LV = cell(length(name(:,1)),1);
    for IDX = 1:length(name(:,1));
        
        filename = name{IDX,:};     % Current image name
        
        %         if strcmpi(stain,'IF');
        %
        %             % !!!!!!!!!!!!!!!!!!!!!!!!!!
        %             % CHECK IF THIS FILENAME MANIPULATION IS STILL APPROPRIATE
        %             !!!!!!!!!!!!!!!!!!!!!!!!!!
        %
        %             sp = strsplit(filename,'_');
        %             strold = sp{1};
        %             for i = 2:4;
        %                 strnew = strcat(strold,'_',sp{i});
        %                 strold = strnew;
        %             end
        %
        %             strnew = strcat(strnew,'_');
        %
        %             fcomp = dir(fullfile(strcat(path,strnew,'*stitch*')));
        %             fcomp = {fcomp.name}';
        %
        %             comp = regexpi(fcomp,'merged');
        %
        %             for i = 1:length(comp);
        %                 if ~isempty(comp{i})
        %                     mfilename{j} = fcomp{i};
        %                     mergenm = mfilename{j}(comp{i}:comp{i}+length('merged')-1);
        %                     pathmerge = strcat(path,mfilename{j});
        %                 end
        %             end
        %
        %         end
        
                
        % Preprocessing of current image
        if ~strcmpi(proc,'No to All') || (strcmpi(proc,'No to All') && ~exist(strcat(path,groupnm,'/',outnm{IDX,:}),'dir'))
            
            % Determine if current image needs to be processed
            if exist(strcat(path,groupnm,'/',outnm{IDX,:}),'dir')
                I = imread(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'.',ext));
                figure; imshow(I); set(gcf,'Name',char(strcat(outnm{IDX,:},' //',{' '},'Previously Processed Image')));
                proc = questdlg('Adjust processed image?',char(strcat(stain,' //',{' '},'Image',{' '},num2str(IDX))),'Yes','No','No to All','No');      pause(0.1);     close all
            else
                proc = 'Yes';
            end
            
            % Apply preprocessing and threshold background objects
            if strcmpi(proc,'Yes');
                [I,cancel] = threshold(IDX,imthresh);
            end
                        
            % Stop current loop, if requested
            if cancel;   return;   end
            
        end
        
        % Select left ventricle, if required...
        if strcmpi(tissue,'Myocardial Infarction         ') && (~strcmpi(procLV,'No to All') || (strcmpi(procLV,'No to All') && ~exist(strcat(path,groupnm,'/',outnm{IDX,:},'/ventricle-center.mat'),'file')));
            
            % Determine if current image needs to be processed
            if exist(strcat(path,groupnm,'/',outnm{IDX,:},'/ventricle-center.mat'),'file')
                procLV = questdlg('Re-select point inside ventricle?',char(strcat(stain,' //',{' '},'Image',{' '},num2str(IDX))),'Yes','No','No to All','No');      pause(0.1);
            else
                procLV = 'Yes';
            end
            
            % Make and store selection
            if strcmpi(procLV,'Yes');
                I = imread(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'.',ext));
                figure; imshow(I);  set(gcf,'Name','Select lumen of LEFT VENTRICLE');
                jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);
                
                [lv_x,lv_y] = ginput(1);    close(gcf);    LV{IDX,1} = [lv_x,lv_y];
                save(strcat(path,groupnm,'/',outnm{IDX,:},'/ventricle-center.mat'),'LV');
            else
                LVo = importdata(strcat(path,groupnm,'/',outnm{IDX,:},'/ventricle-center.mat'));    LV{IDX,1} = LVo{1,1};
            end
            
        elseif strcmpi(tissue,'Myocardial Infarction         ') && (strcmpi(procLV,'No to All') && exist(strcat(path,groupnm,'/',outnm{IDX,:},'/ventricle-center.mat'),'file'));
            LVo = importdata(strcat(path,groupnm,'/',outnm{IDX,:},'/ventricle-center.mat'));    LV{IDX,1} = LVo{1,1};
        end
        
    end
    
    % Terminate code execution, if requested
    if cancel;   break;   end
    
    % -- Partitioning for variational analysis ----------------------------
    localpart = cell(length(name(:,1)),1);    proc = 'Yes';
    if ~strcmpi(vartype,'none') || (strcmpi(vartype,'none') && strcmpi(tissue,'Myocardial Infarction         '));
        
        for IDX = 1:length(name(:,1));
            
            if IDX == 1;    fprintf(char(strcat('===== Variational Analysis - ',{' '},outnm{IDX,:},' =====\n')));    else    fprintf(char(strcat('\n===== Variational Analysis - ',{' '},outnm{IDX,:},' =====\n')));    end
            
            % Update current filepath and load image
            fpstr = strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:});         I = imread(strcat(fpstr,'.',ext));
            
            % Determine if current image needs to be processed
            if ~strcmpi(proc,'No to All') || (strcmpi(proc,'No to All') && ~exist(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'_thickness.mat'),'file'))
                
                if exist(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'_thickness.mat'),'file');
                    proc = questdlg('Recompute wall thickness?',char(strcat(stain,' //',{' '},'Image',{' '},num2str(IDX))),'Yes','No','No to All','No');      pause(0.1);     close all
                else
                    proc = 'Yes';
                end
                
                % Extract inner and outer contours of cross-section in image
                if strcmpi(proc,'Yes');
                    
                    fprintf(strcat('  Extracting Boundaries'));    [Ibw, bprops] = boundaries(I,LV{IDX,1});
                    
                    % Compute local thickness profile
                    fprintf(strcat('    Computing Thickness'));    thick = thickness(IDX,Ibw,bprops);
                    
                    % Partition tissue section to allow for variational analysis
                    if npartr > 1 || npartq > 1
                        
                        % Setup partition labeling
                        if strcmpi(vartype,'circ');    vname = {'      CIRC'};    elseif strcmpi(vartype,'rad');    vname = {'       RAD'};    elseif strcmpi(vartype,'both');    vname = {'  CIRC&RAD'};    end
                        fprintf(char(strcat(vname,{' '},'Partitioning')));
                        
                        % Partition image based on number of r/theta segments
                        localpart{IDX,1} = partition(IDX, Ibw, I, bprops, npartr, npartq);
                        
                    else
                        localpart{IDX,1} = [];
                    end
                    
                else
                    
                    % Load previously saved partitions
                    if npartr > 1 || npartq > 1
                        if strcmpi(tissue,'Myocardial Infarction         ');
                            part = importdata(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'_partition_',char(vartype),'.mat'));
                            localpart{IDX,1} = part.localpart;
                        else
                            localpart{IDX,1} = importdata(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'_partition_',char(vartype),'.mat'));
                        end
                    else
                        localpart{IDX,1} = [];
                    end
                    
                end
                
            else
                
                % Load previously saved partitions
                if npartr > 1 || npartq > 1
                    if strcmpi(tissue,'Myocardial Infarction         ');
                        part = importdata(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'_partition_',char(vartype),'.mat'));
                        localpart{IDX,1} = part.localpart;
                    else
                        localpart{IDX,1} = importdata(strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:},'_partition_',char(vartype),'.mat'));
                    end
                else
                    localpart{IDX,1} = [];
                end
                
            end
        end
    end
    
    
    % -- Layer segmentation -----------------------------------------------
    if layer
        
        % Define names of each layer based on tissue type
        if strcmpi(tissue,'Artery / Vasculature') && nlayer == 2;
            layernm = {'media'; 'adventitia'};
            
        elseif strcmpi(tissue,'Artery / Vasculature') && nlayer == 3;
            layernm = {'intima'; 'media'; 'adventitia'};
            
        elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 2;
            layernm = {'layer1'; 'layer2'};
            
        elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 3;
            layernm = {'layer1'; 'layer2'; 'layer3'};
            
        elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 4;
            layernm = {'layer1'; 'layer2'; 'layer3'; 'layer4'};
            
        elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 5;
            layernm = {'layer1'; 'layer2'; 'layer3'; 'layer4'; 'layer5'};
        end
        
        % Initialize single segmentation modification flag array
        singlechk = zeros(length(name(:,1)),length(layernm)-1);
        
        
        % Attempt automatic segmentation, if applicable -------------------
        procopt = 1;      chkseg = 'no';
        for IDX = 1:length(name(:,1));
            
            fprintf(char(strcat('\n======== Layer Analysis - ',{' '},outnm{IDX,:},' ========\n')));
            
            % Update current filepath
            fpstr = strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:});
            
            % Define conditions in which to perform segmentation
            if nlayer == 2;
                cond = (~exist(strcat(fpstr,'_',layernm{1},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{2},'.',ext),'file'));
            elseif nlayer == 3;
                cond = (~exist(strcat(fpstr,'_',layernm{1},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{2},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{3},'.',ext),'file'));
            elseif nlayer == 4;
                cond = (~exist(strcat(fpstr,'_',layernm{1},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{2},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{3},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{4},'.',ext),'file'));
            elseif nlayer == 5;
                cond = (~exist(strcat(fpstr,'_',layernm{1},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{2},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{3},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{4},'.',ext),'file') || ~exist(strcat(fpstr,'_',layernm{5},'.',ext),'file'));
            end
            
            % Conditional segmentation based on number of layers
            for numseg = 1:nlayer - 1;
                if cond
                    
                    % Load images to perform segmentation
                    if strcmpi(tissue,'Artery / Vasculature');
                        
                        if numseg == 2;
                            fprintf(strcat('Segmenting: Intima/Media'));
                            if exist(strcat(fpstr,'_',layernm{numseg},'.',ext),'file')
                                I = imread(strcat(fpstr,'_',layernm{numseg},'.',ext));
                            else
                                I = imread(strcat(fpstr,'.',ext));
                            end
                        elseif numseg == 1;
                            fprintf(strcat('Segmenting: Media/Adventitia'));
                            I = imread(strcat(fpstr,'.',ext));
                        end
                        
                    else
                        fprintf(strcat('Segmenting: Layer',num2str(numseg),'/Layer',num2str(numseg+1)));
                        I = imread(strcat(fpstr,'.',ext));
                    end
                    
                    % Segment layers based on microstructure
                    [singlechk(IDX,numseg),~,~] = layer_segmentation(IDX,numseg,I,chkseg,[],procopt,[],[],[]);
                    
                end
            end
            
        end
        
        
        % Manually perform (or adjust) segmentation -----------------------
        procopt = 1;    allno = 'Yes';    chkseg = 'yes';   mode = {'Adjusting'; 'Compiling'};
        for IDX = 1:length(name(:,1));

            % Adjust results of 'automatic' segmentation
            for numseg = 1:nlayer-1;
                
                if strcmpi(tissue,'Artery / Vasculature');
                    if numseg == nlayer-1;
                        fprintf(strcat(mode{1},': Media/Adventitia'));
                    else
                        fprintf(strcat(mode{1},': Intima/Media'));
                    end
                else
                    fprintf(char(strcat(mode{1},':',{' '},upper(layernm{end-numseg}(1)),layernm{end-numseg}(2:end),'/',upper(layernm{end-numseg+1}(1)),layernm{end-numseg+1}(2:end))));
                end
                
                [~,allno,~] = layer_segmentation(IDX,numseg,I,chkseg,[],procopt,[],singlechk(IDX,numseg),allno);
                
            end
            
        end
        
        
        % Compile and save results of 'automatic' segmentation ------------
        for IDX = 1:length(name(:,1));
            
            % Parameters for saving segmentations
            procopt = 2;    allno = 'No to All';
            
            % Update current filepath and load current image
            fpstr = strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:});     I = imread(strcat(fpstr,'.',ext));
            
            % Initialize array for storage of calculated wall percentages
            wptemp = zeros(nlayer-1,5);
            
            % Compile results of adjusted segmentation
            for numseg = 1:nlayer-1;
                
                if strcmpi(tissue,'Artery / Vasculature');
                    if numseg == nlayer-1;
                        fprintf(strcat(mode{2},': Media/Adventitia'));
                    else
                        fprintf(strcat(mode{2},': Intima/Media'));
                    end
                else
                    fprintf(char(strcat(mode{2},':',{' '},upper(layernm{end-numseg}(1)),layernm{end-numseg}(2:end),'/',upper(layernm{end-numseg+1}(1)),layernm{end-numseg+1}(2:end))));
                end
                
                [~,~,wptemp(numseg,:)] = layer_segmentation(IDX,numseg,I,chkseg,[],procopt,[],singlechk(IDX,numseg),allno);
                
            end
            
            % Organize wall percentage arrays
            if strcmpi(tissue,'Artery / Vasculature') && nlayer == 2;
                wp = wptemp;
                
            elseif strcmpi(tissue,'Artery / Vasculature') && nlayer == 3;
                wp = cat(2,wptemp(1,1),wptemp(2,1:3),wptemp(1,4),wptemp(2,4:5));
                
            elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 2;
                wp = wptemp;
                
            elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 3;
                wp = cat(2,wptemp(2,1),wptemp(1,1:3),wptemp(2,4),wptemp(1,4:5));
                
            elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 4;
                wp = cat(2,wptemp(3,1:2),wptemp(1,1:3),wptemp(3,4:5),wptemp(1,4:5));
                
            elseif ~strcmpi(tissue,'Artery / Vasculature') && nlayer == 5;
                wp = cat(2,wptemp(4,1),wptemp(3,1:2),wptemp(1,1:3),wptemp(4,4),wptemp(3,4:5),wptemp(1,4:5));
            end
            
            fprintf(strcat('Saving Wall Percentages'));
            
            % Save wall percentages for each processed image
            [~,~,~] = layer_segmentation(IDX,[],I,[],wp,procopt,fileform,[],[]);
            
        end
        
        clear wp wptemp
        
        % Calculate regional variation in layer percentages ---------------
        if ~strcmpi(vartype,'none');
            
            chkseg = 'var';    th = 1;
            for IDX = 1:length(name(:,1));
                
                % Parameters for manual segmentation
                procopt = 3;    allno = 'Yes';
                
                % Update current filepath and load image
                fpstr = strcat(path,groupnm,'/',outnm{IDX,:},'/',outnm{IDX,:});         I = imread(strcat(fpstr,'.',ext));
                
                % Setup partition labeling
                if strcmpi(vartype,'circ');    vname = {'CIRC'};    elseif strcmpi(vartype,'rad');    vname = {'RAD'};    elseif strcmpi(vartype,'both');    vname = 'CIRC&RAD';    end
                
                fprintf(char(strcat(vname,{' '},'Wall Percentage')));
                
                % Compile local wall percentages for each image
                [~,~,wp{IDX,1}] = layer_segmentation(IDX,[],I,chkseg,[],procopt,[],[],[]);
                
                fprintf(strcat('Saving Local Wall Percentages'));
                
                % Save local wall percentages
                [~,~,~] = layer_segmentation(IDX,[],I,[],wp,procopt,fileform,[],[]);
                
            end
            
        end
        
    end
    
    
    % -- Colorimetric analyses ---------------------------------------------
    parfor IDX = 1:length(name(:,1));
        
        fprintf(char(strcat('\n======== Color Analysis - ',{' '},outnm{IDX,:},' ========\n')));
        
        [allstack, afdata, HSLav] = color_segmentation(IDX,[],[],fileform,1, path, groupnm, outnm, tissue, stain, ext, SF, layer, layernm, rembleb, imsave, vartype, enhance);
        
        if layer
            for LDX = 1:nlayer
                [allstack, afdata, HSLav] = color_segmentation(IDX,LDX,[],fileform,1, path, groupnm, outnm, tissue, stain, ext, SF, layer, layernm, rembleb, imsave, vartype, enhance);
            end
        end
        
        if ~isempty(localpart{IDX,1})
            [allstack, afdata, HSLav] = color_segmentation(IDX,[],localpart{IDX,1},fileform,1, path, groupnm, outnm, tissue, stain, ext, SF, layer, layernm, rembleb, imsave, vartype, enhance);
            
            if layer
                for LDX = 1:nlayer
                    [allstack, afdata, HSLav] = color_segmentation(IDX,LDX,localpart{IDX,1},fileform,1, path, groupnm, outnm, tissue, stain, ext, SF, layer, layernm, rembleb, imsave, vartype, enhance);
                end
            end
        end
        
        fprintf('Completed Processing...%s\n',outnm{IDX,:});  pause(0.1);  close all force;
        
    end
    
    
    % -- Scar morphological analysis ----------------------------
    if npartq > 1 && strcmpi(tissue,'Myocardial Infarction         ');
        
        for IDX = 1:length(name(:,1));
            
            fprintf(char(strcat('\n======== Scar Morphological Analysis - ',{' '},outnm{IDX,:},' ========\n')));
            
            scar_segmentation(IDX,fileform);
        end
        
    end
    
        
    
    % ------- THE COMMENTED SECTION BELOW STILL NEEDS TO BE UPDATED!!! --------
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    %
    % % GROUP VARIATIONAL SAVING...PROCOPT = 4
    %
    %
    %     % Compile and save results from all images defined in batch
    %     if ~cancel
    %         fold = dir(strcat(path,groupnm));
    %         fold = {fold([fold(:).isdir]).name}';
    %         fold = fold(3:end);
    %
    %         for i = 1:2%:type;
    %
    %             if layer
    %
    %                 for n = 1:length(layers);
    %
    %                     wp = [];
    %                     eval(strcat(typestr{i},'stack = [];'))
    %                     eval(strcat(typestr{i},'af = [];'))
    %                     eval(strcat(typestr{i},'hsl = [];'))
    %
    %                     for j = 1:length(fold);
    %
    %                         afpth = strcat(path,groupnm,'\',fold{j},'\',stain,'_',fold{j},'_',typestr{i},'_AF_',layers{n},'.mat');
    %                         hslpth = strcat(path,groupnm,'\',fold{j},'\',stain,'_',fold{j},'_',typestr{i},'_HSL_',layers{n},'.mat');
    %                         wppth = strcat(path,groupnm,'\',fold{j},'\',stain,'_',fold{j},'_wall percentage.mat');
    %
    %                         if exist(afpth,'file') > 0
    %                             afdt = importdata(afpth);
    %                         else
    %                             afdt = [];
    %                         end
    %
    %                         if exist(hslpth,'file') > 0
    %                             hsldt = importdata(hslpth);
    %                         else
    %                             hsldt = [];
    %                         end
    %
    %                         if n == 1 && (exist(wppth,'file') > 0)
    %                             wpdt = importdata(wppth);
    %                         elseif n == 1 && (exist(wppth,'file') == 0)
    %                             wpdt = [];
    %                         end
    %
    %                         if i == 1
    %
    %                             if strcmpi(stain,'POL')
    %                                 sz = 1;
    %                             elseif strcmpi(stain,'IHC')
    %                                 sz = 2;
    %                             elseif strcmpi(stain,'VVG') || strcmpi(stain,'MTC')
    %                                 sz = 3;
    %                             elseif strcmpi(stain,'PSR')
    %                                 sz = 5;
    %                             elseif strcmpi(stain,'MOV')
    %                                 sz = 6;
    %                             end
    %
    %                             for k = 1:sz;
    %                                 eval(strcat(typestr{i},'stack = cat(4,',typestr{i},'stack,NaN);'));
    %                             end
    %
    %                         elseif i == 2
    %
    %                             tmp = importdata(strcat(path,groupnm,'\',fold{j},'\',stain,'_',fold{j},'_Extracted_Total_',layers{n},'.mat'));
    %
    %                             eval(strcat(typestr{i},'stack.image',num2str(j),' = tmp;'));
    %
    %                         end
    %
    %                         if n == 1
    %                             wp = cat(1,wp,wpdt);
    %                         end
    %
    %                         eval(strcat(typestr{i},'af = cat(1,',typestr{i},'af,afdt);'));
    %                         eval(strcat(typestr{i},'hsl = cat(2,',typestr{i},'hsl,hsldt);'));
    %
    %                     end
    %
    %                     ncol = eval(strcat('size(',typestr{i},'hsl,2);'));
    %                     colsp = ncol/length(fold);
    %
    %                     eval(strcat(typestr{i},'hsltmp = ',typestr{i},'hsl;'));
    %                     for j = 1:colsp;
    %
    %                         tmp = eval(strcat(typestr{i},'hsltmp(:,j);'));
    %
    %                         for k = j+colsp:colsp:ncol;
    %                             tmp = cat(2,tmp,eval(strcat(typestr{i},'hsltmp(:,k);')));
    %                         end
    %
    %                         eval(strcat(typestr{i},'hsl(:,[1:ncol/colsp]+(j-1)*ncol/colsp) = tmp;'));
    %                     end
    %
    %                     outnm = fold;
    %                     analysis = typestr{i};
    %                     layernm = layers{n};
    %
    %                     if n == 1
    %                         media_adventitia([],wp,fileform,0,[]);
    %                     end
    %
    %                     eval(strcat('outputdata(',analysis,'stack,',analysis,'af,',analysis,'hsl,fileform,0);'));
    %
    %                 end
    %
    %             else
    %
    %                 eval(strcat(typestr{i},'stack = [];'))
    %                 eval(strcat(typestr{i},'af = [];'))
    %                 eval(strcat(typestr{i},'hsl = [];'))
    %
    %                 for j = 1:length(fold);
    %
    %                     afpth = strcat(path,groupnm,'\',fold{j},'\',stain,'_',fold{j},'_',typestr{i},'_AF.mat');
    %                     hslpth = strcat(path,groupnm,'\',fold{j},'\',stain,'_',fold{j},'_',typestr{i},'_HSL.mat');
    %
    %                     if exist(afpth,'file') > 0
    %                         afdt = importdata(afpth);
    %                     end
    %
    %                     if exist(hslpth,'file') > 0
    %                         hsldt = importdata(hslpth);
    %                     end
    %
    %                     if i == 1
    %
    %                         if strcmpi(stain,'POL')
    %                             sz = 1;
    %                         elseif strcmpi(stain,'IHC')
    %                             sz = 2;
    %                         elseif strcmpi(stain,'VVG') || strcmpi(stain,'MTC')
    %                             sz = 3;
    %                         elseif strcmpi(stain,'PSR')
    %                             sz = 5;
    %                         elseif strcmpi(stain,'MOV')
    %                             sz = 6;
    %                         end
    %
    %                         for k = 1:sz;
    %                             eval(strcat(typestr{i},'stack = cat(4,',typestr{i},'stack,NaN);'));
    %                         end
    %
    %                     elseif i == 2
    %
    %                         tmp = importdata(strcat(path,groupnm,'\',fold{j},'\',stain,'_',fold{j},'_Extracted_Total.mat'));
    %
    %                         eval(strcat(typestr{i},'stack.image',num2str(j),' = tmp;'));
    %
    %                     end
    %
    %                     eval(strcat(typestr{i},'af = cat(1,',typestr{i},'af,afdt);'));
    %                     eval(strcat(typestr{i},'hsl = cat(2,',typestr{i},'hsl,hsldt);'));
    %
    %                 end
    %
    %                 ncol = eval(strcat('size(',typestr{i},'hsl,2);'));
    %                 colsp = ncol/length(fold);
    %
    %                 eval(strcat(typestr{i},'hsltmp = ',typestr{i},'hsl;'));
    %                 for j = 1:colsp;
    %
    %                     tmp = eval(strcat(typestr{i},'hsltmp(:,j);'));
    %
    %                     for k = j+colsp:colsp:ncol;
    %                         tmp = cat(2,tmp,eval(strcat(typestr{i},'hsltmp(:,k);')));
    %                     end
    %
    %                     eval(strcat(typestr{i},'hsl(:,[1:ncol/colsp]+(j-1)*ncol/colsp) = tmp;'));
    %                 end
    %
    %                 outnm = fold;
    %                 analysis = typestr{i};
    %                 eval(strcat('outputdata(',analysis,'stack,',analysis,'af,',analysis,'hsl,fileform,0);'));
    %             end
    %
    %         end
    %
    %         clc
    %
    %     else
    %         close all force
    %     end
    
    
end
