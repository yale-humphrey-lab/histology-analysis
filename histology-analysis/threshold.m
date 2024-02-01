function [Ith,cancel] = threshold(IDX,imthresh)

% !!!!!!!!!!!!!!!!!!!!!!
% ADD OPTIONS FOR PSEUDOCOLORING OF IF IMAGES...IF CROP, MUST CROP ALL CHANNELS
% !!!!!!!!!!!!!!!!!!!!!!

global stain                                                                  

global path ext                                                                                     

global groupnm outnm rednm mergenm                                                                        

% Open raw images to be analyzed
fname = outnm{IDX,:};
I = imread(strcat(path,fname,'.',ext));      I = I(:,:,1:3);  % In case of alpha channel

% dPSR: Load corresponding brightfield PSR image
if strcmpi(stain,'dPSR');
    pidx = strfind(fname,'polarized');   bfname = [fname(1:pidx-1) 'brightfield'];
    
    if exist(strcat(path,bfname,'.',ext),'file');
        Ibf = imread(strcat(path,bfname,'.',ext));   Ibf = Ibf(:,:,1:3);
    end
end

% Automatic thresholding based on top/bottom-hat operation
if strcmpi(imthresh,'auto');

    % Morphological top/bottom-hat
    if strcmpi(stain,'dPSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
        Ibh1 = imtophat(I(:,:,1),strel('disk',50));        Ibh1(Ibh1<15) = 0;
        Ibh2 = imtophat(I(:,:,2),strel('disk',50));        Ibh2(Ibh2<15) = 0;
        Ibh3 = imtophat(I(:,:,3),strel('disk',50));        Ibh3(Ibh3<15) = 0;
        Ibh = Ibh1 + Ibh2 + Ibh3;                                bground = 0;
    else
        Ibh1 = imbothat(I(:,:,1),strel('disk',50));        Ibh1(Ibh1<15) = 0;
        Ibh2 = imbothat(I(:,:,2),strel('disk',50));        Ibh2(Ibh2<15) = 0;
        Ibh3 = imbothat(I(:,:,3),strel('disk',50));        Ibh3(Ibh3<15) = 0;
        Ibh = Ibh1 + Ibh2 + Ibh3;                              bground = 255;
    end
    
    % Create mask from top/bottom-hat
    Im = im2bw(Ibh,0.01);

    % Remove background in each image channel
    I1 = I(:,:,1);    I1(Im == 0) = bground;
    I2 = I(:,:,2);    I2(Im == 0) = bground;
    I3 = I(:,:,3);    I3(Im == 0) = bground;

    % Recompile thresholded image
    I = cat(3,I1,I2,I3);
    
    % dPSR: Process brightfield PSR image
    if exist('Ibf','var')
        Ibh1 = imbothat(Ibf(:,:,1),strel('disk',50));        Ibh1(Ibh1<15) = 0;
        Ibh2 = imbothat(Ibf(:,:,2),strel('disk',50));        Ibh2(Ibh2<15) = 0;
        Ibh3 = imbothat(Ibf(:,:,3),strel('disk',50));        Ibh3(Ibh3<15) = 0;
        Ibh = Ibh1 + Ibh2 + Ibh3;                              bground = 255;
        
        % Create mask from top/bottom-hat
        Im = im2bw(Ibh,0.01);
        
        % Remove background in each image channel
        I1 = Ibf(:,:,1);    I1(Im == 0) = bground;
        I2 = Ibf(:,:,2);    I2(Im == 0) = bground;
        I3 = Ibf(:,:,3);    I3(Im == 0) = bground;
        
        % Recompile thresholded image
        Ibf = cat(3,I1,I2,I3);
    end
    
end

% Show current image
figure; imshow(I)
jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);


% Check for additional image pre-processing...
done = 0;
while ~done;
    
    if strcmpi(stain,'IF')
        [process,cancel] = histodlg([],'processing-options-immuno');
    else
        [process,cancel] = histodlg([],'processing-options-histo');
    end
    
    % Crop image...
    if process.crop;
         
        % Crop images based on cell size
        hr = imrect(gca,[1 1 0.25*size(I,2) 0.25*size(I,1)]);    set(gcf,'Name','Resize RECTANGLE to crop image');     
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));    setPositionConstraintFcn(hr,fcn);    wait(hr);
        
        % Get coordinates of crop location
        rectpos = round(hr.getPosition);
        
        % Apply cropping rectangle
        I = imcrop(I,[rectpos(1) rectpos(2) rectpos(3)-1 rectpos(4)-1]);
        
        % dPSR: Process brightfield PSR image
        if exist('Ibf','var')
            Ibf = imcrop(Ibf,[rectpos(1) rectpos(2) rectpos(3)-1 rectpos(4)-1]);
        end
        
        if process.inner || process.outer
            close(gcf);    figure;    imshow(I);
            jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);
        end
        
    end
    
    % Manually threshold image...
    if process.manual;
        
        close all force
        
        % Set image size parameters
        scrsz = get(0,'ScreenSize');
        pos = [0.25*scrsz(3) 0.01*scrsz(4) 0.708*scrsz(3) 0.915*scrsz(4)];
        
        figure('Position',pos)
        h1 = subplot('Position',[0.05 0.29 0.92 0.68]);
        imshow(I)
        
        % Initialize colors for lower/upper bound
        Lcolor = [0.35, 0.35, 0.35];        Ucolor = [0, 0, 0];
        
        % RGB color thresholds for bright/darkfield images
        if strcmpi(stain,'dPSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
            Lthresh = 0;        Uthresh = 9; % 20;
        else
            Lthresh = 180;      Uthresh = 220;
        end
        
        % Generate histogram of RGB color values in each image channel
        histRGB(I,Lthresh,Uthresh,Lcolor,Ucolor);
        
        [~,cancelth] = histodlg(I,'threshold-manual');
        
        I = getimage(h1);
        
        if cancelth
            close all force;
        end
        
        % bPSR: manually threshold brightfield image...
        if exist('Ibf','var')
            
            close all force
            
            % Set image size parameters
            scrsz = get(0,'ScreenSize');
            pos = [0.25*scrsz(3) 0.01*scrsz(4) 0.708*scrsz(3) 0.915*scrsz(4)];
            
            figure('Position',pos)
            h1 = subplot('Position',[0.05 0.29 0.92 0.68]);
            imshow(Ibf)
            
            % Initialize colors for lower/upper bound
            Lcolor = [0.35, 0.35, 0.35];        Ucolor = [0, 0, 0];
            
            % RGB color thresholds for bright/darkfield images
            Lthresh = 180;      Uthresh = 220;
            
            % Generate histogram of RGB color values in each image channel
            stain = 'bPSR';   %#ok
            
            histRGB(Ibf,Lthresh,Uthresh,Lcolor,Ucolor);
            
            [~,cancelth] = histodlg(Ibf,'threshold-manual');
            
            Ibf = getimage(h1);
            
            if cancelth
                close all force;
            end
            
            stain = 'dPSR';

        end
        
    end
    
    % Remove area INSIDE of region...
    if process.inner;

        if exist('Ibf','var');   bfmask = true(size(Ibf,1),size(Ibf,2));   end

        rmvchk = 1;
        regionuc = {' INNER'};    regionlc = {' inner'};
        
        remreg = 'Yes';
        while strcmpi(remreg,'Yes');
            
            % Isolate region of interest in thresholded image
            set(gcf,'Name',strcat('Select',regionuc{1},' region to remove'))
            
            % Create mask from traced polygon
%             h = impoly(gca);    wait(h);    mask = h.createMask();
            h = imfreehand(gca);    wait(h);    mask = h.createMask();
            
            % Update mask for brightfield analysis, if dPSR
            if exist('Ibf','var');    bfmask(mask) = false;    end
            
            % Set background color
            if strcmpi(stain,'dPSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
                bground = 0;
            else
                bground = 255;
            end
            
            % Remove background in each image channel
            I1 = I(:,:,1);  I1(mask == rmvchk) = bground;
            I2 = I(:,:,2);  I2(mask == rmvchk) = bground;
            I3 = I(:,:,3);  I3(mask == rmvchk) = bground;
            
            % Recompile cleaned image
            I = cat(3,I1,I2,I3);
            
            % dPSR: Remove background in each brightfield image channel
            if exist('Ibf','var');
                % Remove background in each image channel
                I1 = Ibf(:,:,1);  I1(~bfmask) = 255;
                I2 = Ibf(:,:,2);  I2(~bfmask) = 255;
                I3 = Ibf(:,:,3);  I3(~bfmask) = 255;
                
                % Recompile cleaned image
                Ibf = cat(3,I1,I2,I3);
            end
            
            % Show cleaned image
            close;                          figure;         imshow(I);
            jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);
            
            % Clean image as much as necessary...
            remreg = questdlg(strcat('Remove more of ',regionlc{1},' region?'),'Remove?','Yes','No','No');
            
        end
        
    end
    
    % Remove area OUTSIDE of region...
    if process.outer;
        
        if exist('Ibf','var');   bfmask = true(size(Ibf,1),size(Ibf,2));   end

        rmvchk = 0;
        regionuc = {' OUTER'};    regionlc = {' outer'};
        
        remreg = 'Yes';
        while strcmpi(remreg,'Yes');
            
            % Isolate region of interest in thresholded image
            set(gcf,'Name',strcat('Select',regionuc{1},' region to remove'))
            
            % Create mask from traced polygon
%             h = impoly(gca);    wait(h);    mask = h.createMask();
            h = imfreehand(gca);    wait(h);    mask = h.createMask();
            
            % Update mask for brightfield analysis, if dPSR
            if strcmpi(stain,'dPSR');   bfmask(~mask) = false;   end
            
            % Set background color
            if strcmpi(stain,'dPSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
                bground = 0;
            else
                bground = 255;
            end
            
            % Remove background in each image channel
            I1 = I(:,:,1);  I1(mask == rmvchk) = bground;
            I2 = I(:,:,2);  I2(mask == rmvchk) = bground;
            I3 = I(:,:,3);  I3(mask == rmvchk) = bground;
            
            % Recompile cleaned image
            I = cat(3,I1,I2,I3);
            
            % dPSR: Remove background in each brightfield image channel
            if exist('Ibf','var');
                % Remove background in each image channel
                I1 = Ibf(:,:,1);  I1(~bfmask) = 255;
                I2 = Ibf(:,:,2);  I2(~bfmask) = 255;
                I3 = Ibf(:,:,3);  I3(~bfmask) = 255;
                
                % Recompile cleaned image
                Ibf = cat(3,I1,I2,I3);
            end
            
            % Show cleaned image
            close;                          figure;         imshow(I);
            jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);
            
            % Clean image as much as necessary...
            remreg = questdlg(strcat('Remove more of ',regionlc{1},' region?'),'Remove?','Yes','No','No');
            
        end
                        
    end
       
    % Restore original image...
    if process.reset;
        I = imread(strcat(path,fname,'.',ext));      I = I(:,:,1:3);  % In case of alpha channel
        
        if strcmpi(stain,'dPSR') && exist('Ibf','var')
            Ibf = imread(strcat(path,bfname,'.',ext));   Ibf = Ibf(:,:,1:3);
        end
        
    end
    
    % Close all open figure windows...
    close all force;
    
    % Show processed image as long as processing will continue...
    if ~cancel && any([process.reset, process.manual, process.crop, process.inner, process.outer])
        figure;                         imshow(I)
        jf = get(gcf,'JavaFrame');      pause(0.1);     set(jf,'Maximized',1);
    end
    
    % Either cancel processing or store completed output...
    if cancel
        Ith = NaN;  return
    elseif ~any([process.reset, process.manual, process.crop, process.inner, process.outer])
        Ith = I;    done = 1;
    end
    
end


% Extract name of image
fname = outnm{IDX,:};

% Compile path to save image
savepath = strcat(path,groupnm,'/',fname,'/');

% Create new directory, if needed
if ~isdir(savepath);    mkdir(savepath);    end

% Update filename for saving
if strcmpi(stain,'IF');    mfname = strrep(fname,rednm,mergenm);    else    mfname = fname;    end

% Write thresholded image to directory
if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
    imwrite(Ith,strcat(savepath,mfname,'.',ext),ext,'Quality',100);
else
    imwrite(Ith,strcat(savepath,mfname,'.',ext),ext);
end

% dPSR: Write thresholded brightfield image to directory
if exist('Ibf','var')
    
    % Compile path to save image
    bfsavepath = strcat(path,groupnm,'/',bfname,'/');
    
    % Create new directory, if needed
    if ~isdir(bfsavepath);    mkdir(bfsavepath);    end

    if strcmpi(ext,'jpg') || strcmpi(ext,'jpeg')
        imwrite(Ibf,strcat(bfsavepath,bfname,'.',ext),ext,'Quality',100);
    else
        imwrite(Ibf,strcat(bfsavepath,bfname,'.',ext),ext);
    end
end

close all force

end
