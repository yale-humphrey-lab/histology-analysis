function [Answer,Cancelled] = histodlg(Im,method)

global stain tissue

% Define dialog boxes to be used in different parts of histology analysis code.
%
% Script is built using: inputsdlg.m (https://www.mathworks.com/matlabcentral/fileexchange/25862-inputsdlg--enhanced-input-dialog-box)

if strcmpi(method,'selection');
    
    Title = 'Processing Selection';
    
    Options.Resize = 'off';
    Options.CancelButton = 'on';
    Options.ButtonNames = 'Continue';
    Options.AlignControls = 'on';
    
    Prompt = {};
    Formats = {};
    DefAns = struct([]);
    
    Prompt(end+1,:) = {'                            Type of Images to Analyze...','',''};
    Formats(1,1).type = 'text';
    Formats(1,1).style = 'text';
    Formats(1,1).span = [1 2];
    
    Prompt(end+1,:) = {'Histology (Brightfield/Darkfield)','histoproc',[]};
    Formats(2,1).type = 'check';
    Formats(2,1).style = 'checkbox';
    Formats(2,1).span = [1 1];
    Formats(2,1).labelloc = 'topleft';
    Formats(2,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k+1),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+2),'Enable',switchEnableVar(h(k),0)); @(h,k)set(h(k+3),'Enable',switchEnableVar(h(k),1)) }));
    DefAns(1).histoproc = false;
    
    Prompt(end+1,:) = {'Immunofluorescence','immunoproc',[]};
    Formats(2,2).type = 'check';
    Formats(2,2).style = 'checkbox';
    Formats(2,2).span = [1 1];
    Formats(2,2).labelloc = 'topleft';
    Formats(2,2).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-1),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+2),'Enable',switchEnableVar(h(k),0)); @(h,k)set(h(k+1),'Enable',switchEnableVar(h(k),1)) }));
    DefAns.immunoproc = false;
    
    Prompt(end+1,:) = {'Tissue Type','tissuehist',[]};
    tissueval = {'General'; 'Artery / Vasculature'; 'Atherosclerosis'; 'Myocardial Infarction         '};
    Formats(3,1).type = 'list';
    Formats(3,1).style = 'popupmenu';
    Formats(3,1).format = 'text';
    Formats(3,1).required = 'on';
    Formats(3,1).items = tissueval;
    Formats(3,1).labelloc = 'topleft';
    Formats(3,1).enable = 'off';
    Formats(3,1).size = [0 24];
    Formats(3,1).span = [1 1];
    DefAns.tissuehist = {'Artery / Vasculature'};    %%%% EDIT
    
    Prompt(end+1,:) = {'Tissue Type','tissueif',[]};
    tissueval = {'General'; 'Artery / Vasculature          '};
    Formats(3,2).type = 'list';
    Formats(3,2).style = 'popupmenu';
    Formats(3,2).format = 'text';
    Formats(3,2).required = 'on';
    Formats(3,2).items = tissueval;
    Formats(3,2).labelloc = 'topleft';
    Formats(3,2).enable = 'off';
    Formats(3,2).size = [0 24];
    Formats(3,2).span = [1 1];
    DefAns.tissueif = {'General'};
    
    [Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options,method);
    
    % ==========================================================================================================================
    
elseif strcmpi(method,'processing-options-histo');
    
    Title = 'Processing Options';
    
    Options.Resize = 'off';
    Options.CancelButton = 'on';
    Options.ButtonNames = 'Continue';
    Options.AlignControls = 'on';
    
    Prompt = {};
    Formats = {};
    DefAns = struct([]);
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(1,1).type = 'text';
    Formats(1,1).style = 'text';
    Formats(1,1).span = [1 1];
    
    Prompt(end+1,:) = {'Reset Image','reset',[]};
    Formats(2,1).type = 'check';
    Formats(2,1).style = 'checkbox';
    Formats(2,1).span = [1 1];
    Formats(2,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k+2),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+4),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+5),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+6),'Enable',switchEnableVar(h(k),1)) }));
    DefAns(1).reset = false;
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(3,1).type = 'text';
    Formats(3,1).style = 'text';
    Formats(3,1).span = [1 1];
    
    Prompt(end+1,:) = {'Manual Threshold','manual',[]};
    Formats(4,1).type = 'check';
    Formats(4,1).style = 'checkbox';
    Formats(4,1).span = [1 1];
    Formats(4,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-2),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+2),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+3),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+4),'Enable',switchEnableVar(h(k),1)) }));
    DefAns.manual = false;
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(5,1).type = 'text';
    Formats(5,1).style = 'text';
    Formats(5,1).span = [1 1];
    
    Prompt(end+1,:) = {'Crop Image','crop',[]};
    Formats(6,1).type = 'check';
    Formats(6,1).style = 'checkbox';
    Formats(6,1).span = [1 1];
    Formats(6,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-4),'Enable',switchEnableVar3(h(k),h(k+1),h(k+2),1)); @(h,k)set(h(k-2),'Enable',switchEnableVar3(h(k),h(k+1),h(k+2),1)) }));
    DefAns.crop = false;
    
    Prompt(end+1,:) = {'Remove INSIDE Region','inner',[]};
    Formats(7,1).type = 'check';
    Formats(7,1).style = 'checkbox';
    Formats(7,1).span = [1 1];
    Formats(7,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-5),'Enable',switchEnableVar3(h(k-1),h(k),h(k+1),1)); @(h,k)set(h(k-3),'Enable',switchEnableVar3(h(k-1),h(k),h(k+1),1)) }));
    DefAns.inner = false;
    
    Prompt(end+1,:) = {'Remove OUTSIDE Region','outer',[]};
    Formats(8,1).type = 'check';
    Formats(8,1).style = 'checkbox';
    Formats(8,1).span = [1 1];
    Formats(8,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-6),'Enable',switchEnableVar3(h(k-2),h(k-1),h(k),1)); @(h,k)set(h(k-4),'Enable',switchEnableVar3(h(k-2),h(k-1),h(k),1)) }));
    DefAns.outer = false;
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(9,1).type = 'text';
    Formats(9,1).style = 'text';
    Formats(9,1).span = [1 1];
    
    [Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options,method);
    
    % ==========================================================================================================================
    
elseif strcmpi(method,'processing-options-immuno');
    
    Title = 'Processing Options';
    
    Options.Resize = 'off';
    Options.CancelButton = 'on';
    Options.ButtonNames = 'Continue';
    Options.AlignControls = 'on';
    
    Prompt = {};
    Formats = {};
    DefAns = struct([]);
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(1,1).type = 'text';
    Formats(1,1).style = 'text';
    Formats(1,1).span = [1 1];
    
    Prompt(end+1,:) = {'Reset Image','reset',[]};
    Formats(2,1).type = 'check';
    Formats(2,1).style = 'checkbox';
    Formats(2,1).span = [1 1];
    Formats(2,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k+2),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+4),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+5),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+6),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+8),'Enable',switchEnableVar(h(k),1)) }));
    DefAns(1).reset = false;
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(3,1).type = 'text';
    Formats(3,1).style = 'text';
    Formats(3,1).span = [1 1];
    
    Prompt(end+1,:) = {'Manual Threshold','manual',[]};
    Formats(4,1).type = 'check';
    Formats(4,1).style = 'checkbox';
    Formats(4,1).span = [1 1];
    Formats(4,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-2),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+2),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+3),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+4),'Enable',switchEnableVar(h(k),1)); @(h,k)set(h(k+6),'Enable',switchEnableVar(h(k),1)) }));
    DefAns.manual = false;
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(5,1).type = 'text';
    Formats(5,1).style = 'text';
    Formats(5,1).span = [1 1];
    
    Prompt(end+1,:) = {'Crop Image','crop',[]};
    Formats(6,1).type = 'check';
    Formats(6,1).style = 'checkbox';
    Formats(6,1).span = [1 1];
    Formats(6,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-4),'Enable',switchEnableVar3(h(k),h(k+1),h(k+2),1)); @(h,k)set(h(k-2),'Enable',switchEnableVar3(h(k),h(k+1),h(k+2),1)); @(h,k)set(h(k+4),'Enable',switchEnableVar3(h(k),h(k+1),h(k+2),1)) }));
    DefAns.crop = false;
    
    Prompt(end+1,:) = {'Remove INSIDE Region','inner',[]};
    Formats(7,1).type = 'check';
    Formats(7,1).style = 'checkbox';
    Formats(7,1).span = [1 1];
    Formats(7,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-5),'Enable',switchEnableVar3(h(k-1),h(k),h(k+1),1)); @(h,k)set(h(k-3),'Enable',switchEnableVar3(h(k-1),h(k),h(k+1),1)); @(h,k)set(h(k+3),'Enable',switchEnableVar3(h(k-1),h(k),h(k+1),1)) }));
    DefAns.inner = false;
    
    Prompt(end+1,:) = {'Remove OUTSIDE Region','outer',[]};
    Formats(8,1).type = 'check';
    Formats(8,1).style = 'checkbox';
    Formats(8,1).span = [1 1];
    Formats(8,1).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-6),'Enable',switchEnableVar3(h(k-2),h(k-1),h(k),1)); @(h,k)set(h(k-4),'Enable',switchEnableVar3(h(k-2),h(k-1),h(k),1)); @(h,k)set(h(k+2),'Enable',switchEnableVar3(h(k-2),h(k-1),h(k),1)) }));
    DefAns.outer = false;
    
    Prompt(end+1,:) = {'----------------------------------------------','',''};
    Formats(9,1).type = 'text';
    Formats(9,1).style = 'text';
    Formats(9,1).span = [1 1];
        
    [Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options,method);
    
    % ==========================================================================================================================
    
elseif strcmpi(method,'threshold-manual');
    
    Title = 'Image Threshold';
    
    Options.Resize = 'off';
    Options.CancelButton = 'on';
    Options.ButtonNames = 'Continue';
    Options.AlignControls = 'on';
    
    Prompt = {};
    Formats = {};
    DefAns = struct([]);
    
    Prompt(end+1,:) = {' ', 'Lcolor',[]};
    Formats(2,1).type = 'color';
    Formats(2,1).style = 'pushbutton';
    Formats(2,1).callback.ButtonDownFcn = @(~,~,h,k)histRGB(Im,str2double(get(h(k+2),'String')),str2double(get(h(k+5),'String')),get(h(k),'BackgroundColor'),get(h(k+3),'BackgroundColor'));
    Formats(2,1).size = [35 35];
    Formats(2,1).span = [1 1];
    DefAns(1).Lcolor = [0.35 0.35 0.35];
    
    Prompt(end+1,:) = {'Lower Bound', 'LBraw',[]};
    Formats(2,2).type = 'range';
    Formats(2,2).style = 'slider';
    Formats(2,2).limits = [0 255];
    Formats(2,2).labelloc = 'topleft';
    Formats(2,2).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k+1),'String',num2str(round(get(h(k),'Value')*10)/10)), @(h,k)histRGB(Im,str2double(get(h(k+1),'String')),str2double(get(h(k+4),'String')),get(h(k-1),'BackgroundColor'),get(h(k+2),'BackgroundColor'))}));
    Formats(2,2).span = [1 1];
    
    if sum(strcmp(stain,{'dPSR','POL','EnF','IF'})) > 0
        DefAns.LBraw = 0;
    else
        DefAns.LBraw = 180;
    end
    
    Prompt(end+1,:) = {' ','LB',[]};
    Formats(2,3).type = 'edit';
    Formats(2,3).format = 'float';
    Formats(2,3).size = [100 18];
    Formats(2,3).limits = [0 255];
    Formats(2,3).labelloc = 'topleft';
    Formats(2,3).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-1),'Value',str2double(get(h(k),'String'))), @(h,k)histRGB(Im,str2double(get(h(k),'String')),str2double(get(h(k+3),'String')),get(h(k-2),'BackgroundColor'),get(h(k+1),'BackgroundColor'))}));
    Formats(2,3).span = [1 1];
    DefAns.LB = DefAns.LBraw;
    
    Prompt(end+1,:) = {' ', 'Ucolor',[]};
    Formats(3,1).type = 'color';
    Formats(3,1).style = 'pushbutton';
    Formats(3,1).callback.ButtonDownFcn = @(~,~,h,k)histRGB(Im,str2double(get(h(k-1),'String')),str2double(get(h(k+2),'String')),get(h(k-3),'BackgroundColor'),get(h(k),'BackgroundColor'));
    Formats(3,1).size = [35 35];
    Formats(3,1).span = [1 1];
    DefAns.Ucolor = [0 0 0];
    
    Prompt(end+1,:) = {'Upper Bound', 'UBraw',[]};
    Formats(3,2).type = 'range';
    Formats(3,2).style = 'slider';
    Formats(3,2).limits = [0 255];
    Formats(3,2).labelloc = 'topleft';
    Formats(3,2).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k+1),'String',num2str(round(get(h(k),'Value')*10)/10)), @(h,k)histRGB(Im,str2double(get(h(k-2),'String')),str2double(get(h(k+1),'String')),get(h(k-4),'BackgroundColor'),get(h(k-1),'BackgroundColor'))}));
    Formats(3,2).span = [1 1];
    
    if sum(strcmp(stain,{'dPSR','POL','EnF','IF'})) > 0
        DefAns.UBraw = 9; % 20
    else
        DefAns.UBraw = 220;
    end
    
    Prompt(end+1,:) = {' ','UB',[]};
    Formats(3,3).type = 'edit';
    Formats(3,3).format = 'float';
    Formats(3,3).size = [100 18];
    Formats(3,3).limits = [0 255];
    Formats(3,3).labelloc = 'topleft';
    Formats(3,3).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-1),'Value',str2double(get(h(k),'String'))), @(h,k)histRGB(Im,str2double(get(h(k-3),'String')),str2double(get(h(k),'String')),get(h(k-5),'BackgroundColor'),get(h(k-2),'BackgroundColor'))}));
    Formats(3,3).span = [1 1];
    DefAns.UB = DefAns.UBraw;
    
    Prompt(end+1,:) = {'Apply Threshold','',''};
    Formats(5,1).type = 'button';
    Formats(5,1).style = 'pushbutton';
    Formats(5,1).span = [1 3];
    Formats(5,1).callback = @(~,~,h,k)RGBfilter(Im,str2double(get(h(k-4),'String')),str2double(get(h(k-1),'String')));
    
    [Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options,method);
    
    % ==========================================================================================================================
    
elseif strcmpi(method,'initialize-histo');
    
    Title = 'Initialize values for file input/output';
    
    Options.Resize = 'off';
    Options.CancelButton = 'on';
    Options.ButtonNames = 'Continue';
    Options.AlignControls = 'off';
    
    Prompt = {};
    Formats = {};
    DefAns = struct([]);
    
    if exist(strcat(cd,'\path.mat'),'file') == 0;
        upath = userpath;
        
        if strcmpi(upath(end),';')
            upath = strcat(upath(1:end-1),'\');
        else
            upath = strcat(upath,'\');
        end
        
        ext = 'tif';
        
        files = cellstr(cdfiles(upath,ext));
        
        if isempty(files) == 0;
            DefAns(1).name = files(1,:);
        else
            DefAns(1).name = files;
        end
        
        DefAns.path = upath;
        
    else
        upath = importdata(strcat(cd,'\path.mat'));
        ext = 'tif';
        
        if isdir(upath) == 0;
            upath = userpath;
            
            if strcmpi(upath(end),';')
                upath = strcat(upath(1:end-1),'\');
            else
                upath = strcat(upath,'\');
            end
        end
        
        files = cellstr(cdfiles(upath,ext));
        
        
        if isempty(files) == 0;
            DefAns(1).name = files(1,:);
        else
            DefAns(1).name = files;
        end
        
        DefAns.path = upath;
    end
    
    Prompt(end+1,:) = {'Image Filenames','name',[]};
    Formats(1,2).type = 'list';
    Formats(1,2).style = 'listbox';
    Formats(1,2).format = 'text';
    Formats(1,2).limits = [0 9];
    Formats(1,2).required = 'on';
    Formats(1,2).size = [250 255];
    Formats(1,2).items = files;
    Formats(1,2).labelloc = 'topleft';
    Formats(1,2).callback = @(~,~,h,k)set(h(k+3),'String',fileadjust(h(k)));
    Formats(1,2).span = [6 1];
    
    Prompt(end+1,:) = {'File Path to Location of Images','path',[]};
    Formats(1,3).type = 'edit';
    Formats(1,3).format = 'dir';
    Formats(1,3).required = 'on';
    Formats(1,3).size = [0 22];
    Formats(1,3).labelloc = 'topleft';
    Formats(1,3).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k),'String',pathadjust(h(k))); @(h,k)set(h(k-1),'Value',1); @(h,k)set(h(k+2),'String',{' '}); @(h,k)set(h(k-1),'String',cdfiles(get(h(k),'String'),handlevalue(h(k+7))))}));
    Formats(1,3).span = [1 6];
    
    Prompt(end+1,:) = {'Output Group Name','groupnm',[]};
    Formats(2,3).type = 'edit';
    Formats(2,3).format = 'text';
    Formats(2,3).required = 'on';
    Formats(2,3).size = [174 22];
    Formats(2,3).labelloc = 'topleft';
    Formats(2,3).span = [1 2];
    
    Prompt(end+1,:) = {'Output File Names','outnm',[]};
    Formats(2,5).type = 'edit';
    Formats(2,5).format = 'text';
    Formats(2,5).limits = [0 9];
    Formats(2,5).required = 'on';
    Formats(2,5).size = [205 165];
    Formats(2,5).labelloc = 'topleft';
    Formats(2,5).callback = @(~,~,h,k)filecheck(h(k),h(k-3));
    Formats(2,5).span = [4 1];
    
    Prompt(end+1,:) = {'Output Format','fileform',[]};
    outval = {'.mat';'.txt';'.xls'};
    Formats(2,6).type = 'list';
    Formats(2,6).style = 'listbox';
    Formats(2,6).format = 'text';
    Formats(2,6).limits = [0 2];
    Formats(2,6).required = 'on';
    Formats(2,6).items = outval;
    Formats(2,6).size = [3*25.2 6.5*18];
    Formats(2,6).labelloc = 'topleft';
    Formats(2,6).span = [3 1];
    DefAns.fileform = {'.mat';'.xls'};    %%%% EDIT
    
    Prompt(end+1,:) = {'Magnification','mago',[]};
    magval = {'60x';'40x';'20x';'10x';'4x';'1x';'Other'};
    Formats(2,7).type = 'list';
    Formats(2,7).style = 'popupmenu';
    Formats(2,7).format = 'text';
    Formats(2,7).required = 'on';
    Formats(2,7).size = [0 22];
    Formats(2,7).items = magval;
    Formats(2,7).labelloc = 'topleft';
    Formats(2,7).callback = @(~,~,h,k)set(h(k+1),'String',checkOtherMag(h(k)));
    Formats(2,7).span = [1 1];
    DefAns.mago = {'20x'};
    
    Prompt(end+1,:) = {' ','mag',[]};
    Formats(2,8).type = 'edit';
    Formats(2,8).format = 'text';
    Formats(2,8).size = [0 22];
    Formats(2,8).labelloc = 'topcenter';
    Formats(2,8).span = [1 1];
    Formats(2,8).enable = 'off';
    DefAns.mag = DefAns.mago;
    
    Prompt(end+1,:) = {'Histological Stain','stain',[]};
    
    if strcmpi(tissue,'General') 
        stainval = {'H&E';'VVG';'MTC';'PSR';'MOV';'IHC';'VnK';'AzR';'POL'};
    elseif strcmpi(tissue,'Artery / Vasculature')
        stainval = {'VVG';'MTC';'PSR';'MOV';'IHC'};
    elseif strcmpi(tissue,'Atherosclerosis')
        stainval = {'EnF';'ORO';'MTC'};
    elseif strcmpi(tissue,'Myocardial Infarction         ')
        stainval = {'MTC';'PSR'};
    end
    
    Formats(3,3).type = 'list';
    Formats(3,3).style = 'listbox';
    Formats(3,3).format = 'text';
    Formats(3,3).limits = [0 1];
    Formats(3,3).required = 'on';
    Formats(3,3).items = stainval;
    Formats(3,3).size = 7.*[11.5 16];
    Formats(3,3).labelloc = 'topleft';
    
    if strcmpi(tissue,'Artery / Vasculature') || strcmpi(tissue,'Myocardial Infarction         ');
        Formats(3,3).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k+6),'Enable',checkEnableStain(h(k),handlevalue(h(k)),1)); @(h,k)set(h(k+7),'Enable',checkEnableStain(h(k),handlevalue(h(k)),1)); @(h,k)set(h(k+11),'Enable',checkEnableStain(h(k+6),handlevalue(h(k)),0)); @(h,k)set(h(k+12),'Enable',checkEnableStain(h(k+7),handlevalue(h(k)),0)) }));
    end
    
    Formats(3,3).span = [3 1];
    
    if strcmpi(tissue,'General') 
        DefAns.stain = {'MOV'};    %%%% EDIT
    elseif strcmpi(tissue,'Artery / Vasculature')
        DefAns.stain = {'MOV'};    %%%% EDIT
    elseif strcmpi(tissue,'Atherosclerosis')
        DefAns.stain = {'EnF'};
    elseif strcmpi(tissue,'Myocardial Infarction         ')
        DefAns.stain = {'MTC'};
    end
    
    Prompt(end+1,:) = {'File Extension','ext',[]};
    extval = {'tif';'tiff';'jpg';'jpeg';'bmp';'png'};
    Formats(3,4).type = 'list';
    Formats(3,4).style = 'listbox';
    Formats(3,4).format = 'text';
    Formats(3,4).limits = [0 1];
    Formats(3,4).required = 'on';
    Formats(3,4).items = extval;
    Formats(3,4).size = 7.*[11.5 16];
    Formats(3,4).labelloc = 'topleft';
    Formats(3,4).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-8),'String',cdfiles(get(h(k-7),'String'),handlevalue(h(k)))); @(h,k)set(h(k-5),'String',{' '})}));
    Formats(3,4).span = [3 1];
    DefAns.ext = {'tif'};
    
    Prompt(end+1,:) = {'Resolution','reso',[]};
    resval = {'4080 x 3072';'2040 x 1536';'1600 x 1200';'1360 x 1024';'680 x 512';'Other'};
    Formats(3,7).type = 'list';
    Formats(3,7).style = 'popupmenu';
    Formats(3,7).format = 'text';
    Formats(3,7).required = 'on';
    Formats(3,7).items = resval;
    Formats(3,7).size = [0 22];
    Formats(3,7).labelloc = 'topleft';
    Formats(3,7).callback = @(~,~,h,k)set(h(k+1),'String',checkOtherRes(h(k)));
    Formats(3,7).span = [1 1];
    DefAns.reso = {'2040 x 1536'};    %%%% EDIT
    
    Prompt(end+1,:) = {' ','res',[]};
    Formats(3,8).type = 'edit';
    Formats(3,8).format = 'text';
    Formats(3,8).labelloc = 'topleft';
    Formats(3,8).size = [0 22];
    Formats(3,8).span = [1 1];
    Formats(3,8).enable = 'off';
    DefAns.res = DefAns.reso;
    
    Prompt(end+1,:) = {' ----------  Variational Analysis  ----------','',''};
    Formats(4,7).type = 'text';
    Formats(4,7).style = 'text';
    Formats(4,7).span = [1 2];
    
    Prompt(end+1,:) = {'Layer Seg.','layer',[]};
    Formats(5,6).type = 'check';
    Formats(5,6).style = 'checkbox';
    Formats(5,6).callback = @(~,~,h,k)set(h(k+5),'Enable',checkEnableVar(h(k)));
    Formats(5,6).span = [1 1];
    Formats(5,6).enable = 'on';
    DefAns.layer = false;
    
    Prompt(end+1,:) = {'Radial','radvar',[]};
    Formats(5,7).type = 'check';
    Formats(5,7).style = 'checkbox';
    Formats(5,7).callback = @(~,~,h,k)set(h(k+5),'Enable',checkEnableVar(h(k)));
    Formats(5,7).span = [1 1];
    
    if strcmpi(tissue,'General') || strcmpi(tissue,'Atherosclerosis')
        Formats(5,7).enable = 'off';
    end
    
    DefAns.radvar = false;
    
    Prompt(end+1,:) = {'Circumferential','circvar',[]};
    Formats(5,8).type = 'check';
    Formats(5,8).style = 'checkbox';
    Formats(5,8).callback = @(~,~,h,k)set(h(k+5),'Enable',checkEnableVar(h(k)));
    Formats(5,8).span = [1 1];
    
    if strcmpi(tissue,'General') || strcmpi(tissue,'Atherosclerosis')
        Formats(5,8).enable = 'off';
    end
    
    DefAns.circvar = false;
    
    Prompt(end+1,:) = {'Save Individual Constituent Images','imsave',[]};
    Formats(6,3).type = 'check';
    Formats(6,3).style = 'checkbox';
    Formats(6,3).span = [1 2];
    Formats(6,3).enable = 'off';
    Formats(6,3).labelloc = 'topleft';
    DefAns.imsave = true;
    
    Prompt(end+1,:) = {'Background Subtraction','imthresh',[]};
    thval = {'Automatic Thresholding                      '; 'Manual Thresholding'};
    Formats(6,5).type = 'list';
    Formats(6,5).style = 'popupmenu';
    Formats(6,5).format = 'text';
    Formats(6,5).limits = [0 1];
    Formats(6,5).required = 'on';
    Formats(6,5).items = thval;
    Formats(6,5).span = [1 1];
    Formats(6,5).size = [0 22];
    Formats(6,5).labelloc = 'topleft';
    DefAns.imthresh = {'Manual Thresholding'};    %%%% EDIT
    
    Prompt(end+1,:) = {'# Layers','nlayer',[]};
    Formats(6,6).type = 'edit';
    Formats(6,6).format = 'integer';
    
    if strcmpi(tissue,'Artery / Vasculature');
        Formats(6,6).limits = [0 3];
    else
        Formats(6,6).limits = [0 5];
    end
    
    Formats(6,6).limits = [0 5];
    Formats(6,6).size = [77 18];
    Formats(6,6).labelloc = 'topleft';
    Formats(6,6).span = [1 1];
    Formats(6,6).enable = 'off';
    DefAns.nlayer = 2;    %%%% EDIT
    
    Prompt(end+1,:) = {'# Rad. Partitions','npr',[]};
    Formats(6,7).type = 'edit';
    Formats(6,7).format = 'integer';
    Formats(6,7).limits = [0 20];
    Formats(6,7).size = [77 18];
    Formats(6,7).labelloc = 'topleft';
    Formats(6,7).span = [1 1];
    Formats(6,7).enable = 'off';
    DefAns.npr = 1;
    
    Prompt(end+1,:) = {'# Circ. Partitions','npq',[]};
    Formats(6,8).type = 'edit';
    Formats(6,8).format = 'integer';
    Formats(6,8).limits = [0 1000];
    Formats(6,8).size = [77 18];
    Formats(6,8).labelloc = 'topleft';
    Formats(6,8).span = [1 1];
    Formats(6,8).enable = 'off';
    DefAns.npq = 1;
    
    [Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options,method);
    
    if ~Cancelled
        
        % Update output names array...
        Answer.outchk = Answer.outnm;    nmnum = length(Answer.outchk(:,1));    fidx = NaN(nmnum,1);
        for II = 1:nmnum;
            fidx(II,1) = length(Answer.outchk(II,:));    fchk = Answer.outnm(II,fidx(II,1));
            while strcmpi(fchk,' ')
                fidx(II,1) = fidx(II,1) - 1;    fchk = Answer.outchk(II,fidx(II,1));
            end
        end
        
        Answer.outnm = cell(nmnum,1);
        for II = 1:nmnum;
            Answer.outnm{II,1} = Answer.outchk(II,1:fidx(II,1));       
        end
        
        Answer = rmfield(Answer,'outchk');
        
        % Type of analysis to perform...
        if Answer.radvar && Answer.circvar;
            Answer.vartype = {'both'};
        elseif ~Answer.radvar && Answer.circvar;
            Answer.vartype = {'circ'};
        elseif Answer.radvar && ~Answer.circvar;
            Answer.vartype = {'rad'};
        elseif ~Answer.radvar && ~Answer.circvar;
            Answer.vartype = {'none'};
        end
        
        % Type of thresholding to perform...
        autochk = strfind(Answer.imthresh,'Automatic');
        manuchk = strfind(Answer.imthresh,'Manual');
        
        if ~isempty(autochk{1});
            Answer.imthresh = {'auto'};
        elseif ~isempty(manuchk{1});
            Answer.imthresh = {'manual'};
        end
        
        % Always save .mat files...
        if ~ismember('.mat',Answer.fileform)
            Answer.fileform{end+1,1} = '.mat';
        end
        
        % Extract magnification and resolution scale factors...
        resnm = char(Answer.res);
        xind = strfind(resnm,'x');
        ycur = str2double(resnm(1:xind-2));
        xcur = str2double(resnm(xind+2:end));
        
        ydef = 1360;
        xdef = 1024;
        Answer.SFr = (xcur/xdef)*(ycur/ydef);
        
        magnm = char(Answer.mag);
        xind = strfind(magnm,'x');
        mcur = str2double(magnm(1:xind-1));
        
        mdef = 20;
        Answer.SFm = mcur/mdef;
        
    end
    
    % ==========================================================================================================================
    
elseif strcmpi(method,'initialize-immuno');
    
    Title = 'Initialize values for file input/output';
    
    Options.Resize = 'off';
    Options.CancelButton = 'on';
    Options.ButtonNames = 'Continue';
    Options.AlignControls = 'off';
    
    Prompt = {};
    Formats = {};
    DefAns = struct([]);
    
    if exist(strcat(cd,'\path.mat'),'file') == 0;
        upath = userpath;
        
        if strcmpi(upath(end),';')
            upath = strcat(upath(1:end-1),'\');
        else
            upath = strcat(upath,'\');
        end
        
        ext = 'tif';
        
        files = cellstr(cdfiles_if(upath,[],ext));
        
        if isempty(files) == 0;
            DefAns(1).name = files(1,:);
        else
            DefAns(1).name = files;
        end
        
        DefAns.path = upath;
        
    else
        upath = importdata(strcat(cd,'\path.mat'));
        ext = 'tif';
        
        if isdir(upath) == 0;
            upath = userpath;
            
            if strcmpi(upath(end),';')
                upath = strcat(upath(1:end-1),'\');
            else
                upath = strcat(upath,'\');
            end
        end
        
        files = cellstr(cdfiles_if(upath,[],ext));
        
        
        if isempty(files) == 0;
            DefAns(1).name = files(1,:);
        else
            DefAns(1).name = files;
        end
        
        DefAns.path = upath;
    end
    
    Prompt(end+1,:) = {'Image Filenames','name',[]};
    Formats(1,2).type = 'list';
    Formats(1,2).style = 'listbox';
    Formats(1,2).format = 'text';
    Formats(1,2).limits = [0 9];
    Formats(1,2).required = 'on';
    Formats(1,2).size = [250 255];
    Formats(1,2).items = files;
    Formats(1,2).labelloc = 'topleft';
    Formats(1,2).callback = @(~,~,h,k)set(h(k+9),'String',fileadjust(h(k)));
    Formats(1,2).span = [6 1];
    
    Prompt(end+1,:) = {'File Path to Location of Images','path',[]};
    Formats(1,3).type = 'edit';
    Formats(1,3).format = 'dir';
    Formats(1,3).required = 'on';
    Formats(1,3).labelloc = 'topleft';
    Formats(1,3).size = [0 22];
    Formats(1,3).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k),'String',pathadjust(h(k))); @(h,k)set(h(k-1),'Value',1); @(h,k)set(h(k+8),'String',{' '}); @(h,k)set(h(k-1),'String',cdfiles_if(get(h(k),'String'),get(h(k+2),'String'),handlevalue(h(k+7))))}));
    Formats(1,3).span = [1 6];
    
    Prompt(end+1,:) = {'Output Group Name','groupnm',[]};
    Formats(2,3).type = 'edit';
    Formats(2,3).format = 'text';
    Formats(2,3).required = 'on';
    Formats(2,3).size = [174 22];
    Formats(2,3).labelloc = 'topleft';
    Formats(2,3).span = [1 2];
    
    Prompt(end+1,:) = {'Red Channel Name','rednm',[]};
    Formats(2,5).type = 'edit';
    Formats(2,5).format = 'text';
    Formats(2,5).required = 'on';
    Formats(2,5).size = [205 22];
    Formats(2,5).labelloc = 'topleft';
    Formats(2,5).callback = @(~,~,h,k)set(h(k-3),'String',cdfiles_if(get(h(k-2),'String'),get(h(k),'String'),handlevalue(h(k+5))));
    Formats(2,5).span = [1 1];
    
    Prompt(end+1,:) = {'Output Format','fileform',[]};
    outval = {'.mat';'.txt';'.xls'};
    Formats(2,6).type = 'list';
    Formats(2,6).style = 'listbox';
    Formats(2,6).format = 'text';
    Formats(2,6).limits = [0 2];
    Formats(2,6).required = 'on';
    Formats(2,6).items = outval;
    Formats(2,6).size = [3*25.2 4*18];
    Formats(2,6).labelloc = 'topleft';
    Formats(2,6).span = [2 1];
    DefAns.fileform = {'.mat'};
    
    Prompt(end+1,:) = {'Magnification','mago',[]};
    magval = {'60x';'40x';'20x';'10x';'4x';'Other'};
    Formats(2,7).type = 'list';
    Formats(2,7).style = 'popupmenu';
    Formats(2,7).format = 'text';
    Formats(2,7).required = 'on';
    Formats(2,7).size = [0 22];
    Formats(2,7).items = magval;
    Formats(2,7).labelloc = 'topleft';
    Formats(2,7).callback = @(~,~,h,k)set(h(k+1),'String',checkOtherMag(h(k)));
    Formats(2,7).span = [1 1];
    DefAns.mago = {'20x'};
    
    Prompt(end+1,:) = {' ','mag',[]};
    Formats(2,8).type = 'edit';
    Formats(2,8).format = 'text';
    Formats(2,8).labelloc = 'topcenter';
    Formats(2,8).size = [0 22];
    Formats(2,8).span = [1 1];
    Formats(2,8).enable = 'off';
    DefAns.mag = DefAns.mago;
    
    Prompt(end+1,:) = {'Histological Stain','stain',[]};
    stainval = {'IF'};
    Formats(3,3).type = 'list';
    Formats(3,3).style = 'listbox';
    Formats(3,3).format = 'text';
    Formats(3,3).limits = [0 1];
    Formats(3,3).required = 'on';
    Formats(3,3).items = stainval;
    Formats(3,3).size = 7.*[11.5 16];
    Formats(3,3).labelloc = 'topleft';
    Formats(3,3).span = [3 1];
    DefAns.stain = {'IF'};
    
    Prompt(end+1,:) = {'File Extension','ext',[]};
    extval = {'tif';'tiff';'jpg';'jpeg';'bmp';'png'};
    Formats(3,4).type = 'list';
    Formats(3,4).style = 'listbox';
    Formats(3,4).format = 'text';
    Formats(3,4).limits = [0 1];
    Formats(3,4).required = 'on';
    Formats(3,4).items = extval;
    Formats(3,4).size = 7.*[11.5 16];
    Formats(3,4).labelloc = 'topleft';
    Formats(3,4).callback = @(~,~,h,k)(cellfun(@(x)feval(x,h,k),{@(h,k)set(h(k-8),'String',cdfiles_if(get(h(k-7),'String'),get(h(k-5),'String'),handlevalue(h(k)))); @(h,k)set(h(k+1),'String',{' '})}));
    Formats(3,4).span = [3 1];
    DefAns.ext = {'tif'};
    
    Prompt(end+1,:) = {'Output File Names','outnm',[]};
    Formats(3,5).type = 'edit';
    Formats(3,5).format = 'text';
    Formats(3,5).limits = [0 9];
    Formats(3,5).required = 'on';
    Formats(3,5).size = [205 120];
    Formats(3,5).labelloc = 'topleft';
    Formats(3,5).callback = @(~,~,h,k)filecheck(h(k),h(k-9));
    Formats(3,5).span = [3 1];
    
    Prompt(end+1,:) = {'Resolution','reso',[]};
    resval = {'4080 x 3072';'2040 x 1536';'1360 x 1024';'680 x 512';'Other'};
    Formats(3,7).type = 'list';
    Formats(3,7).style = 'popupmenu';
    Formats(3,7).format = 'text';
    Formats(3,7).required = 'on';
    Formats(3,7).size = [0 22];
    Formats(3,7).items = resval;
    Formats(3,7).labelloc = 'topleft';
    Formats(3,7).callback = @(~,~,h,k)set(h(k+1),'String',checkOtherRes(h(k)));
    Formats(3,7).span = [1 1];
    DefAns.reso = {'1360 x 1024'};
    
    Prompt(end+1,:) = {' ','res',[]};
    Formats(3,8).type = 'edit';
    Formats(3,8).format = 'text';
    Formats(3,8).labelloc = 'topleft';
    Formats(3,8).size = [0 22];
    Formats(3,8).span = [1 1];
    Formats(3,8).enable = 'off';
    DefAns.res = DefAns.reso;
    
    Prompt(end+1,:) = {'Red Elastin?','redelastin',[]};
    Formats(4,6).type = 'check';
    Formats(4,6).style = 'checkbox';
    Formats(4,6).labelloc = 'topleft';
    DefAns.redelastin = false;
    
    Prompt(end+1,:) = {' ----------  Variational Analysis  ----------','',''};
    Formats(4,7).type = 'text';
    Formats(4,7).style = 'text';
    Formats(4,7).span = [1 2];
    
    Prompt(end+1,:) = {'Layer Seg.','layer',[]};
    Formats(5,6).type = 'check';
    Formats(5,6).style = 'checkbox';
    Formats(5,6).callback = @(~,~,h,k)set(h(k+5),'Enable',checkEnableVar(h(k)));
    Formats(5,6).span = [1 1];
    Formats(5,6).enable = 'on';
    DefAns.layer = false;
    
    Prompt(end+1,:) = {'Radial','radvar',[]};
    Formats(5,7).type = 'check';
    Formats(5,7).style = 'checkbox';
    Formats(5,7).callback = @(~,~,h,k)set(h(k+5),'Enable',checkEnableVar(h(k)));
    Formats(5,7).span = [1 1];
    
    if strcmpi(tissue,'General');
        Formats(5,7).enable = 'off';
    end
    
    DefAns.radvar = false;
    
    Prompt(end+1,:) = {'Circumferential','circvar',[]};
    Formats(5,8).type = 'check';
    Formats(5,8).style = 'checkbox';
    Formats(5,8).callback = @(~,~,h,k)set(h(k+5),'Enable',checkEnableVar(h(k)));
    Formats(5,8).span = [1 1];
    
    if strcmpi(tissue,'General');
        Formats(5,8).enable = 'off';
    end
    
    DefAns.circvar = false;
    
    Prompt(end+1,:) = {'Save Individual Constituent Images','imsave',[]};
    Formats(6,3).type = 'check';
    Formats(6,3).style = 'checkbox';
    Formats(6,3).span = [1 2];
    Formats(6,3).enable = 'off';
    Formats(6,3).labelloc = 'topleft';
    DefAns.imsave = true;
    
    Prompt(end+1,:) = {'Background Subtraction','imthresh',[]};
    thval = {'Automatic Thresholding                      '; 'Manual Thresholding'};
    Formats(6,5).type = 'list';
    Formats(6,5).style = 'popupmenu';
    Formats(6,5).format = 'text';
    Formats(6,5).limits = [0 1];
    Formats(6,5).required = 'on';
    Formats(6,5).items = thval;
    Formats(6,5).span = [1 1];
    Formats(6,5).size = [0 22];
    Formats(6,5).labelloc = 'topleft';
    DefAns.imthresh = {'Automatic Thresholding                      '};
    
    Prompt(end+1,:) = {'# Layers','nlayer',[]};
    Formats(6,6).type = 'edit';
    Formats(6,6).format = 'integer';
    
    if strcmpi(tissue,'Artery / Vasculature');
        Formats(6,6).limits = [0 3];
    else
        Formats(6,6).limits = [0 5];
    end
    
    Formats(6,6).size = [77 18];
    Formats(6,6).labelloc = 'topleft';
    Formats(6,6).span = [1 1];
    Formats(6,6).enable = 'off';
    DefAns.nlayer = 1;
    
    Prompt(end+1,:) = {'# Rad. Partitions','npr',[]};
    Formats(6,7).type = 'edit';
    Formats(6,7).format = 'integer';
    Formats(6,7).limits = [0 20];
    Formats(6,7).size = [77 18];
    Formats(6,7).labelloc = 'topleft';
    Formats(6,7).span = [1 1];
    Formats(6,7).enable = 'off';
    DefAns.npr = 1;
    
    Prompt(end+1,:) = {'# Circ. Partitions','npq',[]};
    Formats(6,8).type = 'edit';
    Formats(6,8).format = 'integer';
    Formats(6,8).limits = [0 1000];
    Formats(6,8).size = [77 18];
    Formats(6,8).labelloc = 'topleft';
    Formats(6,8).span = [1 1];
    Formats(6,8).enable = 'off';
    DefAns.npq = 1;
    
    [Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options,method);
    
    if ~Cancelled
        
        % Update output names array...
        Answer.outchk = Answer.outnm;    nmnum = length(Answer.outchk(:,1));    fidx = NaN(nmnum,1);
        for II = 1:nmnum;
            fidx(II,1) = length(Answer.outchk(II,:));    fchk = Answer.outnm(II,fidx(II,1));
            while strcmpi(fchk,' ')
                fidx(II,1) = fidx(II,1) - 1;    fchk = Answer.outchk(II,fidx(II,1));
            end
        end
        
        Answer.outnm = cell(nmnum,1);
        for II = 1:nmnum;
            Answer.outnm{II,1} = Answer.outchk(II,1:fidx(II,1));
        end
        
        Answer = rmfield(Answer,'outchk');
        
        % Type of analysis to perform...
        if Answer.radvar && Answer.circvar;
            Answer.vartype = {'both'};
        elseif ~Answer.radvar && Answer.circvar;
            Answer.vartype = {'circ'};
        elseif Answer.radvar && ~Answer.circvar;
            Answer.vartype = {'rad'};
        elseif ~Answer.radvar && ~Answer.circvar;
            Answer.vartype = {'none'};
        end
        
        % Type of thresholding to perform...
        autochk = strfind(Answer.imthresh,'Automatic');
        manuchk = strfind(Answer.imthresh,'Manual');
        
        if ~isempty(autochk{1});
            Answer.imthresh = {'auto'};
        elseif ~isempty(manuchk{1});
            Answer.imthresh = {'manual'};
        end
        
        % Always save .mat files...
        if ~ismember('.mat',Answer.fileform)
            Answer.fileform{end+1,1} = '.mat';
        end
        
        % Extract magnification and resolution scale factors...
        resnm = char(Answer.res);
        xind = strfind(resnm,'x');
        ycur = str2double(resnm(1:xind-2));
        xcur = str2double(resnm(xind+2:end));
        
        ydef = 1360;
        xdef = 1024;
        Answer.SFr = (xcur/xdef)*(ycur/ydef);
        
        magnm = char(Answer.mag);
        xind = strfind(magnm,'x');
        mcur = str2double(magnm(1:xind-1));
        
        mdef = 20;
        Answer.SFm = mcur/mdef;
        
    end
    
end


% ---------------------- Nested Functions ---------------------------------
    function files = cdfiles(upath,ext)
        
        % Extract names of files with correct extension in directory
        files = dir(fullfile(strcat(upath,'*.',ext)));
        files = {files.name}';
        
    end


    function files = cdfiles_if(upath,rednm,ext)
        
        % Extract names of files with correct extension in directory
        if isempty(rednm)
            files = dir(fullfile(strcat(upath,'*',rednm,'stitch*.',ext)));
        else
            files = dir(fullfile(strcat(upath,'*',rednm,'*stitch*.',ext)));
        end
        files = {files.name}';
        
    end


    function path = pathadjust(h)
        
        upath = get(h,'String');
        
        if strcmp(upath(end),'\') == 0
            path = strcat(upath,'\');
        else
            path = upath;
        end
        
    end


    function val = handlevalue(h)
        
        string = get(h,'String');
        value = get(h,'Value');
        
        val = char(string(value,:));
    end


    function enable = checkEnableStain(h,stain,type)
        
        val = get(h,'Value');
        
        if type == 0;
            if strcmpi(stain,'POL');
                enable = 'off';
            elseif ~strcmpi(stain,'POL');
                if val == 1
                    enable = 'on';
                elseif val == 0;
                    enable = 'off';
                end
            end
        elseif type == 1;
            if strcmpi(stain,'POL');
                enable = 'off';
            elseif ~strcmpi(stain,'POL');
                enable = 'on';
            end
        end
        
    end


    function enable = checkEnableVar(h)
        
        val = get(h,'Value');
        
        if val == 1;
            enable = 'on';
        elseif val == 0;
            enable = 'off';
        end
        
    end


    function enable = switchEnableVar(h,dir)
        
        val = get(h,'Value');
        
        if val == 1;
            if dir == 1;
                enable = 'off';
            elseif dir == 0;
                enable = 'on';
            end
        elseif val == 0;
            if dir == 1;
                enable = 'on';
            elseif dir == 0;
                enable = 'off';
            end
        end
        
    end


    function enable = switchEnableVar3(h1,h2,h3,dir)
        
        val1 = get(h1,'Value');
        val2 = get(h2,'Value');
        val3 = get(h3,'Value');
        
        if any([val1,val2,val3]);
            if dir == 1;
                enable = 'off';
            elseif dir == 0;
                enable = 'on';
            end
        elseif ~any([val1,val2,val3]);
            if dir == 1;
                enable = 'on';
            elseif dir == 0;
                enable = 'off';
            end
        end
        
    end


    function flnm = fileadjust(h)
        
        val = handlevalue(h);
        
        if isempty(val) == 0;
            ind = strfind(cellstr(val),'.');
            
            for i = 1:length(val(:,1));
                valtemp = char(val(i,:));
                
                indvec = cell2mat(ind(i));
                index = length(indvec);
                flnm(i,:) = cellstr(valtemp(1:indvec(index)-1));
            end
        else
            flnm = {' '};
        end
        
        flnm = char(flnm);
    end


    function filecheck(h1,h2)
        
        val = get(h1,'String');
        num1 = length(val(:,1));
        
        val = handlevalue(h2);
        num2 = length(val(:,1));
        
        if num1 ~= num2
            h = errordlg('Must have same number of ''Output File Names'' as selected ''Image Filenames''','Missing Output File Name(s)','modal');
            uiwait(h)
        end
        
    end


    function label = checkOtherMag(h)
        
        str = handlevalue(h);
        
        if strcmpi(str,'other') == 1;
            label = strcat(inputdlg('Magnification:','Custom',1),'x');
        else
            label = str;
        end
        
    end


    function label = checkOtherRes(h)
        
        str = handlevalue(h);
        
        if strcmpi(str,'other') == 1;
            rsz = inputdlg({'X-Resolution:';'Y-Resolution'},'Custom',1);
            label = strcat(rsz(1),{' x '},rsz(2));
        else
            label = str;
        end
        
    end


    function [Answer,Canceled] = inputsdlg(Prompt, Title, Formats, DefAns, Options, Step)
        %INPUTSDLG Enhanced input dialog box supporting multiple data types
        % ANSWER = INPUTSDLG(PROMPT) creates a modal dialog box that returns user
        % input for multiple prompts in the cell array ANSWER. PROMPT is a 1-D
        % cell array containing the PROMPT strings.
        %
        % Alternatively, PROMPT can have up to 4 columns. The first column
        % sppecifies the prompt string. The second column to specify the struct
        % field names to output ANSWER as a structure. The third column specifies
        % units (i.e., post-fix labels to the right of controls) to display. The
        % fourth column specifies the tooltip string. The tooltip string is ignored
        % for text type.
        %
        % INPUTSDLG uses UIWAIT to suspend execution until the user responds.
        %
        % ANSWER = INPUTSDLG(PROMPT,NAME) specifies the title for the dialog.
        %
        % Note that INPUTSDLG(PROMPT) & INPUTSDLG(PROMPT,NAME) are similar to the
        % standard INPUTDLG function, except for the dialog layout.
        %
        % ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS) can be used to specify the type
        % of parameters to display with FORMATS matrix of structures. The dimension
        % of FORMATS defines how PROMPT items are laid out in the dialog box. For
        % example, if PROMPT has 6 elements and the size of FORMATS is 2x3 then,
        % the items are shown in 2 rows, 3 columns format. The items in PROMPT
        % correspond to a horizontal traversal of FORMATS.
        %
        % The fields in FORMATS structure are:
        %
        %   type     - Type of control ['check',{'edit'},'list','range','text',
        %                               'color','table','button','none']
        %   style    - UI control type used. One of:
        %               {'checkbox'},                for 'check' type
        %               {'edit'}                     for 'edit' type
        %               {'listbox','popupmenu','radiobutton','togglebutton'}
        %                                            for 'list' type
        %               {'slider'}                   for 'range' type
        %               {'text'}                     for 'text' type
        %               {'edit'}                     for 'color' type
        %               {'pushbutton'}               for 'button' and 'color' types
        %               {'table'}                    for 'table' type
        %   format   - Data format: ['string','date','float','integer','logical',
        %                            'vector','file','dir']
        %   limits   - [min max] (see below for details)
        %   required -  'on'   - control must have an answer
        %              {'off'} - control may return empty answer
        %   items    - Type 'edit', Format 'file': File flter spec
        %              Type 'list': Selection items (cell of strings)
        %              Type 'table': Names of columns (cell of strings)
        %   size     - [width height] in pixels. Set to <=0 to auto-size.
        %   enable   - Defines how to respond to mouse button clicks, including which
        %              callback routines execute. One of:
        %               {'on'}      - UI control is operational.
        %                'inactive' � UI control is not operational, but looks the
        %                             same as when Enable is on.
        %                'off'      � UI uicontrol is not operational and its image
        %                             is grayed out.
        %   margin  -  A scalar or 2-element vector specifying the margin between control
        %              and its labels in pixels.
        %   labelloc - Prompt label location:
        %               {'lefttop'}   - left of control, aligned to top
        %                'leftmiddle' - left of control, aligned to middle
        %                'leftbottom' - left of control, aligned to bottom
        %                'topleft'    - above control, aligned left
        %                'topcenter'  - above control, aligned center
        %                'topright'   - above control, aligned right
        %   unitsloc - Units label location:
        %               {'righttop'}     - right of control, aligned to top
        %                'rightmiddle'   - right of control, aligned to middle
        %                'rightbottom'   - right of control, aligned to bottom
        %                'bottomleft'    - below control, aligned left
        %                'bottomcenter'  - below control, aligned center
        %                'bottomright'   - below control, aligned right
        %   callback - Defines callback funcion, a routine that executes whenever
        %              you activate the uicontrol object. For the controls with
        %              separate dialog, their callback functions are executed after
        %              the dialog is closed. The callback function must be given as
        %              a function handle with following syntax:
        %
        %                 my_callbackfcn(hobj,evt,handles,k)
        %
        %              where hobj and evt are the passed-through standard MATLAB
        %              callback arguments, handles is a Nx3 array of dialog
        %              objects. Here, the n-th row corresponds to the n-th PROMPT,
        %              and handles(n,1) is the calling object handle (i.e., same as
        %              hobj). handles(n,2) are the prompt texts and handles(n,3)
        %              are the prompt unit texts.
        %
        %              For example, Formats(n,m).callback.ButtonDownFcn sets the
        %              the button-press callback function.
        %   span     - Defines size of objects in fields [rows columns]
        %
        % A missing field (either missing from FORMATS struct or the field value is
        % left empty for an element) will be filled with a default field value.
        %
        % FORMATS type field defines what type of prompt item to be shown.
        %
        %   type  Description
        %   -------------------------------------------------------------------
        %   edit     Standard edit box (single or multi-line mode)
        %   check    Check box for boolean item
        %   list     Chose from a list of items ('listbox' style allows multiple item
        %            selection)
        %   range    Use slider to chose a value over a range
        %   text     Static text (e.g., for instructions)
        %   color    Color selection using 'uisetcolor'
        %   button   Execute function defined in 'callback'
        %   table    Uitable
        %   none     A placeholder. May be used for its neighboring item to extend
        %            over multiple columns or rows (i.e., "to merge cells")
        %
        % The allowed data format depends on the type of the field:
        %
        %   type    allowed format
        %   --------------------------------------------
        %   check     {logical}, integer
        %   edit      {text}, date, float, integer, file, dir, vector
        %   list      {integer}, text
        %   range     {float}
        %   color     {float}, integer
        %   table     (a cell string to specify ColumnFormat of uiTable)
        %   button    any data format allowed
        %
        % Formats 'file' and 'dir' for 'edit' type uses the standard UIGETFILE,
        % UIPUTFILE, and UIGETDIR functions to retrieve a file or directory name.
        %
        % The role of limits field varies depending on other parameters:
        %
        %   style         role of limits
        %   ---------------------------------------------------
        %   checkbox      If data format is integer, limits(1) is the ANSWER value
        %                 if the check box is not selected  box is not selected and
        %                 limits(2) is the ANSWER if the check box is selected.
        %   edit:text
        %                 If diff(limits)>0, text aligns with the prompt label. If
        %                 diff(limits)<0, tet aligns with the control.
        %   edit::date
        %                 limits must be a free-format date format string or a
        %                 scalar value specifying the date format. Supported format
        %                 numbers are: 0,1,2,6,13,14,15,16,23. The default date
        %                 format is 2 ('mm/dd/yy'). See the tables in DATESTR help
        %                 for the format definitions. As long as the user entry is
        %                 a valid date/time expression, the dialog box
        %                 automatically converts to the assigned format.
        %   edit::float, edit::integer
        %                 This style defines the range of allowed values
        %   edit::vector
        %                 limits specifies the allowed number of elements in a
        %                 vector
        %   edit::file
        %                 If 0<diff(limits)<=1 uses UIGETFILE in single select
        %                 mode with single-line edit. If diff(limits)>1 uses
        %                 UIGETFILE in multi-select mode with multi-line edit. If
        %                 diff(limits)<=0 usees UIPUTFILE with single-line edit
        %   list::listbox If diff(limits)>1, multiple items can be selected. If
        %                 auto-height and limits(1)>0, at most limits(1) lines will
        %                 be shown.
        %   slider        limits(1) defines the smallest value while
        %                 limits(2) defines the largest value
        %   color         'float' format limit must be [0 1] and 'integer' format
        %                 limit must be [0 255]. If not, the behavior is not
        %                 determined.
        %   table         cell array defining ColumnWidths
        %   none          If diff(limits)==0 space is left empty (default)
        %                 If diff(limits)>0 : lets the item from left to extend
        %                 If diff(limits)<0 : lets the item from above to extend
        %                 NOTE: Use 'span' field of a control to automatically set
        %                 the extension mode.
        %
        % Similar to how PROMPT strings are laid out, when FORMATS.style is set to
        % either 'radiobutton' or 'togglebutton', FORMATS.items are laid out
        % according to the dimension of FORMATS.items.
        %
        % There are two quick format options as well:
        %
        %  Quick Format Option 1 (mimicing INPUTDLG behavior):
        %   FORMATS can specify the number of lines for each edit-type prompt in
        %   FORMATS. FORMATS may be a constant value or a column vector having
        %   one element per PROMPT that specifies how many lines per input field.
        %   FORMATS may also be a matrix where the first column specifies how
        %   many rows for the input field and the second column specifies how
        %   many columns wide the input field should be.
        %
        %  Quick Format Option 2:
        %   FORMATS can specify the types of controls and use their default
        %   configurations. This option, however, cannot be used to specify
        %   'list' control as its items are not specified. To use this option,
        %   provide a string (if only 1 control) or a cell array of strings. If
        %   a cell array is given, its dimension is used for the dialog
        %   layout.
        %
        % ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS,DEFAULTANSWER) specifies the
        % default answer to display for each PROMPT. For a non-tiled layout,
        % DEFAULTANSWER must contain the same number of elements as PROMPT (that
        % are not of 'none' style). If PROMPT does not provide ANSWER structure
        % fields, DEFAULTANSWER should be a cell array with element type
        % corresponding to FORMATS.format. Leave the cell element empty for a
        % prompt with 'text' type. If ANSWER is a structure, DEFAULTANSWER must be
        % a struct with the specified fields. (If additional fields are present in
        % DEFAULTANSWER, they will be returned as parts of ANSWER.)
        %
        % For edit::file controls, a default answer that does not correspond to an
        % existing file will be used as a default path and/or file name in the
        % browse window.  It is passed as the DefaultName parameter to UIGETFILE or
        % UIPUTFILE.
        %
        % To enable Tiled Mode, FORMATS must be given as a vector and DEFAULTANSWER
        % must be given as a cell matrix or a struct vector. If FORMATS is a row
        % vector, the dialog controls are tiled vertically; conversely if it is a
        % column vector, the controls are tiled horizontally. If DEFAULTANSWER is
        % given as a cell matrix, the number of rows of DEFAULTANSWER must match
        % the number of elements of PROMPT. Each column of DEFAULTANSWER forms a
        % tile row/column. If DEFAULTANSWER is given as a struct vector, each
        % struct element forms a tile row/column.
        %
        % ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS,DEFAULTANSWER,OPTIONS) specifies
        % additional options. If OPTIONS is the string 'on', the dialog is made
        % resizable. If OPTIONS is a structure, the fields recognized are:
        %
        %  Option Field Description {} indicates the default value
        %  ----------------------------------------------------------------------
        %  Resize        Make dialog resizable: 'on' | {'off'}
        %  WindowStyle   Sets dialog window style: {'normal'} | 'modal'
        %  Interpreter   Label text interpreter: 'latex' | {'tex'} | 'none'
        %  CancelButton  Show Cancel button: {'on'} | 'off'
        %  ApplyButton   Adds Apply button: 'on' | {'off'}
        %  Sep           Space b/w prompts in pixels: {10}
        %  ButtonNames   Customize OK|Cancel|Apply button names: {up to 3 elements}
        %  AlignControls Align adjacent controls in the same column: 'on' | {'off'}
        %  FontSize      Customize font size. Default: get(0,'DefaultUicontrolFontSize')
        %  CreateFcn     Callback function executed right after dialog creation
        %                with syntax my_createfcn(hobj,evt,handles) with
        %                standard MATLAB callback arguments, hobj & evt, and Nx3
        %                array of handles. Here, the n-th row corresponds to the
        %                n-th PROMPT, and handles(n,1) is the calling object handle
        %                (i.e., same as hobj). handles(n,2) are the prompt texts
        %                and handles(n,3) are the prompt unit texts.
        %  DeleteFcn     Callback function executed just before deleting dialog.
        %                The function syntax is the same as CreateFcn.
        %
        % [ANSWER,CANCELED] = INPUTSDLG(...) returns CANCELED = TRUE if user
        % pressed Cancel button, closed the dialog, or pressed ESC. In such event,
        % the content of ANSWER is set to the default values.
        %
        % Note on Apply Button feature. Pressing the Apply button makes the current
        % change permanent. That is, pressing Cancel button after pressing Apply
        % button only reverts ANSWER back to the states when the Apply button was
        % pressed last. Also, if user pressed Apply button, CANCELED flag will not
        % be set even if user canceled out of the dialog box.
        %
        % Examples:
        %
        % prompt={'Enter the matrix size for x^2:';'Enter the colormap name:'};
        % name='Input for Peaks function';
        % formats(1) = struct('type','edit','format','integer','limits',[1 inf]);
        % formats(2) = struct('type','edit','format','text','limits',[0 1]);
        % defaultanswer={20,'hsv'};
        %
        % [answer,canceled] = inputsdlg(prompt,name,formats,defaultanswer);
        %
        % formats(2).size = -1; % auto-expand width and auto-set height
        % options.Resize='on';
        % options.WindowStyle='normal';
        % options.Interpreter='tex';
        %
        % answer = inputsdlg(prompt,name,formats,defaultanswer,options);
        %
        % prompt(:,2) = {'Ndim';'Cmap'};
        % defaultanswer = struct(defaultanswer,prompt(:,2),1);
        %
        % answer = inputsdlg(prompt,name,formats,defaultanswer,options);
        %
        % See also INPUTDLG, DIALOG, ERRORDLG, HELPDLG, LISTDLG, MSGBOX,
        %  QUESTDLG, UIGETFILE, UIPUTFILE, UIGETDIR, DATESTR.
        
        % Version 2.0 (July 17 2013)
        % Written by: Takeshi Ikuma
        % Contributors: Andreas Greuer, Luke Reisner, Florian Hatz
        % Created: Nov. 16, 2009
        % Revision History:
        %  v.1.1 (Nov. 19, 2009)
        %  * Fixed bugs (reported by AG):
        %   - not returning Canceled output
        %   - erroneous struct output behavior
        %   - error if all row elements of a column are auto-expandable
        %  * Added Apply button option
        %  * Added support for Units (label to the right of controls)
        %  * Updated the help text
        %  v.1.11 (Nov. 20, 2009)
        %  * Fixed bugs (reported by AG):
        %   - incorrect Canceled output when Cancel button is pressed
        %  v.1.12 (Nov. 20, 2009)
        %  * Fixed bugs (reported by AG):
        %   - again incorrect Canceled output behavior
        %  v.1.2 (May 20, 2010)
        %  * Fixed bugs (reported by AG & Jason):
        %   - Apply button->Canel button does not revert back to post-apply answers.
        %   - Line 265 handles.Figure -> handles.fig
        %  * Added edit::date support
        %  * Added formats.enable support
        %  * Added options.CancelButton support
        %  * Added options.ButtonNames support
        %  v.1.2.1 (June 11, 2010)
        %  * Fixed default option bug (reported by Jason)
        %  v.1.2.2 (July 15, 2010)
        %  * Rewritten checkoptions() (to correct issues reported by Jason)
        %  * Bug Fix: file & dir control enable config were interpreted backwards
        %  v.1.2.3 (July 19, 2010)
        %  * checkoptions() bug fix (to correct issues reported by Kevin)
        %  v.1.3 (August 13, 2010, by Luke Reisner)
        %  * Improved dialog layout:
        %   - Less wasted space, better control distribution, more consistent margins
        %   - Buttons are right-aligned per OS standards
        %  * Changed edit::date to return a simple date vector (see DATEVEC help)
        %  * Added support for free-form date format specifiers to edit::date
        %  * Added ability to limit the number of displayed lines for a listbox
        %  * Added ability to set default browse path/filename for edit::file controls
        %  * Added options.AlignControls to align adjacent controls in the same column
        %  * Added options.UnitsMargin to control spacing between controls and units
        %  * Fixed bugs:
        %   - Flickering or misplaced controls when dialog first appears
        %   - Radiobutton and togglebutton controls couldn't be disabled
        %   - Edit::integer controls allowed non-integer values
        %   - Slider controls didn't auto-size properly
        %   - Other minor miscellaneous bugs
        %  v.2.0 (July 17, 2013, by T. Ikuma & F. Hatz)
        %  * PROMPT(:,4) to specify tooltip strings
        %  * Enabled Tiled Mode with DEFAULTANSWER & FORMATS dimension specs. See
        %    inputsdlg_demo_struct
        %  * Added types: table, color, button
        %  * Added formats: logical and vector
        %  * Added 'text' format for 'list' type
        %  * Added 'callback' format field
        %  * Added 'CreateFcn' option field
        %  * Added 'DeleteFcn' option field
        %  * Added 'required' format field
        %  * Added 'labelloc' and 'unitsloc' format fields to customize location of
        %    the labels
        %  * removed UnitMargin option field and added 'margin' format field.
        %  * Added 'span' format field to make spanning across multiple rows and
        %    columns simpler
        %  * Improved handling of 'Formats' (new variable 'span' / see inputsdlg_demo)
        %  * A dialog with non-editable text only displays 'OK' button and does not
        %    return any argument
        %  * edit:file: diff(formats.limits)==0 => uiputfile
        %  v.2.0.1 (Aug 05, 2013, by T. Ikuma)
        %  * Bug fix on inputdlg compatible calls
        %  v.2.0.2 (Aug 06, 2013, by T. Ikuma)
        %  * Version compatibility fix for pre-R2013a
        %  v.2.0.3 (Sep 07, 2013)
        %  * Bug fix on parsing popup menu callback (reported by E. Morales)
        
        % to-do list
        % * Add edit:font for the built-in uisetfont dialog
        % * Support for DefaultXXX option fields to set default Formats field values
        % * Auto-layout given figure's desired width
        % * Limits-based option for positioning of text-type
        % * Better spanning support
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% # of argument Check %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        narginchk(0,6);
        nargoutchk(0,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Handle Input Args %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin<1, Prompt={}; end
        if nargin<2, Title = ''; end
        if nargin<3, Formats=struct([]); end
        if nargin<4, DefAns = {}; end
        if nargin<5, Options = struct([]); end
        
        % Check Prompt input
        [Prompt,FieldNames,Units,TooltipStrings,err] = checkprompt(Prompt);
        if ~isempty(err), error(err{:}); end
        NumQuest = numel(Prompt); % number of prompts
        
        if isempty(Title)
            Title = ' ';
        elseif iscellstr(Title)
            Title = Title{1}; % take the first entry
        elseif ~ischar(Title)
            error('inputsdlg:InvalidInput','Title must be a string of cell string.');
        end
        
        % make sure that the Options is valid
        [Options,FormatDefaultFields,err] = checkoptions(Options);
        if ~isempty(err), error(err{:}); end
        
        % make sure that the Formats structure is valid & fill it in default values
        % as needed
        [Formats,err] = checkformats(Formats,NumQuest,FormatDefaultFields);
        if ~isempty(err), error(err{:}); end
        
        % make sure that the DefAns is valid & set Answer using DefAns and default
        % values if DefAns not given
        [DefaultAnswer,TileMode,err] = checkdefaults(DefAns,Formats,FieldNames);
        if ~isempty(err), error(err{:}); end
        Answer = DefaultAnswer;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Create Dialog GUI %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % lay contents out on a dialog box
        [Formats,handles,sinfo] ...
            = buildgui(Title,Prompt,Units,TooltipStrings,FieldNames,Formats,TileMode,Options,Step);
        
        % fill in the default answers, setup callbacks,
        initgui();
        
        IsRequired = arrayfun(@(fmt)strcmp(fmt.required,'on'),Formats);
        DefaultReqMet = checkReq(DefaultAnswer);
        ReqMet = false; % also modified by the Apply button
        
        Applied = false; % set true by pressing Apply Button
        Canceled = ~ishghandle(handles.fig);
        
        % Go into uiwait if the figure handle is still valid.
        % This is mostly the case during regular use.
        try
            while ~(Canceled || ReqMet) % exit only if Canceled or all the rerequired fields are filled
                
                % Wait till uiresume is called
                uiwait(handles.fig);
                
                % Check handle validity again since figure could be deleted externally
                Canceled = strcmp(get(handles.fig,'UserData'),'Cancel');
                
                if Canceled % return the default answer
                    Answer = DefaultAnswer; % revert back to the default answer
                else
                    Answer = getAnswer(); % get the final answers
                    ReqMet = checkReq(Answer);
                    if ~ReqMet
                        h = errordlg('All required parameters must be filled.','Missing Required Value(s)','modal');
                        uiwait(h);
                    end
                    
                end
            end
            
            % if user deletefcn defined, call it now
            if ~isempty(Options.DeleteFcn)
                Options.DeleteFcn(handles.fig,[],handles.ctrls);
            end
            
            % Close the figure if it's still open
            delete(handles.fig);
            
        catch ME
            if ishghandle(handles.fig)
                delete(handles.fig);
                ME.getReport
                throw(ME);
            else
                error('Inputsdlg dialog window was closed externally');
            end
        end
        
        % If Canceled, convert Canceled to integer depending on Applied condition
        Canceled = Canceled * (Canceled + (Applied && DefaultReqMet));
        % 0 - OK pressed
        % 1 - Canceled
        % 2 - Canceled but prior to it Applied button has been pressed, filled in
        %     all the required answers
        
        % If Tiled, reshape the answer
        Answer = reshape(Answer,NumQuest,max(TileMode));
        
        % If FieldNames given, convert Answer to struct
        Answer = selectivecell2struct(Answer,FieldNames);
        
        % If dialog contains only non-editable texts, return w/o output argument
        if all(sinfo.istext)
            clear Answer
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% NESTED FUNCTIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function initgui
            
            fig = handles.fig;
            
            % set OK button as the default
            fh = handle(fig);
            fh.setDefaultButton(handles.btns(1));
            
            % Set callbacks of figure/control buttons
            set(fig, 'UserData', 'Cancel',...
                'WindowKeyPressFcn', @doFigureKeyPress,...
                'ResizeFcn', @(~,~)resizegui(handles,sinfo,Options),...
                'CloseRequestFcn', @(hObj,evd)doOKCancel(hObj,evd,false));
            for k = 1 : numel(handles.btns)
                hbtn = handles.btns(k);
                switch get(hbtn, 'UserData')
                    case 'OK'
                        set(hbtn, 'Callback', @(hObj,evd)doOKCancel(hObj,evd,true));
                    case 'Cancel'
                        set(hbtn, 'Callback', @(hObj,evd)doOKCancel(hObj,evd,false));
                    case 'Apply'
                        set(hbtn, 'Callback', @doApply);
                end
            end
            
            % Set callback functions and set default values
            for k = 1:numel(Answer)
                
                h = handles.ctrls(k,1);
                val = Answer{k};
                fmt = Formats(k);
                cbnames = fieldnames(fmt.callback);
                
                ena = {'Enable',fmt.enable};
                fcnname = 'Callback';
                fcn = {};
                aux = {};
                
                % get callback function handles for custom callbacks to be
                % called from inputsdlg's default callback
                [tf,I] = ismember('callback',lower(cbnames));
                if tf
                    cbfcn = fmt.callback.(cbnames{I});
                else
                    cbfcn = {};
                end
                [tf,I] = ismember('buttondownfcn',lower(cbnames));
                if tf
                    bdfcn = fmt.callback.(cbnames{I});
                else
                    bdfcn = {};
                end
                
                switch fmt.style
                    case 'edit'
                        ansname = 'String';
                        switch fmt.format
                            case {'integer','float'}
                                % for numeric edit box, check for the range & set mouse down behavior
                                fcn = @(hObj,evd)checkNumericRange(hObj,evd,k,cbfcn);
                                aux = {'UserData',val}; % save the numeric data
                                val = num2str(val);
                            case 'date'
                                fcn = @(hObj,evd)checkDate(hObj,evd,k,cbfcn);
                            case 'vector'
                                % for vector edit box, check for the range & set mouse down behavior
                                fcn = @(hObj,evd)checkVector(hObj,evd,k,cbfcn);
                                val = num2str(val(:));
                            case 'file'
                                mode = diff(fmt.limits);
                                fcnname = 'ButtonDownFcn';
                                if strcmp(ena{2},'on')
                                    ena{2} = 'inactive';
                                    fcn = @(hObj,evd)openFilePrompt(hObj,evd,k,bdfcn);
                                end
                                
                                val = cellstr(val);
                                if ~isempty(val)
                                    dirname = fileparts(val{1});
                                    if ~isdir(dirname), dirname = ''; end
                                else
                                    dirname = '';
                                end
                                aux = {'UserData',dirname};
                                
                                if mode <= 1 % single-file
                                    val = val{1};
                                end
                            case 'dir'
                                fcnname = 'ButtonDownFcn';
                                if strcmp(ena{2},'on')
                                    ena{2} = 'inactive';
                                    fcn = @(hObj,evd)openDirPrompt(hObj,evd,k,bdfcn);
                                end
                        end
                    case 'pushbutton'
                        if strcmp(fmt.type,'color')
                            if strcmp(ena{2},'on')
                                fcn = @(hObj,evd)openColorPrompt(hObj,evd,k,bdfcn);
                            end
                            ansname = 'BackgroundColor';
                        else
                            ansname = 'UserData';
                        end
                    case {'radiobutton', 'togglebutton'}
                        fcnname = 'SelectionChangeFcn';
                        ansname = 'SelectedObject';
                        if strcmp(fmt.format,'integer')
                            hbtn = get(h,'UserData'); % hButtons
                            val = hbtn(val);
                        else
                            hbtn = get(h,'UserData'); % hButtons
                            val = findobj(hbtn,'flat','String',val);
                        end
                        set(hbtn,ena{:});
                        ena = {};
                    case 'table'
                        ansname = 'Data';
                        fcnname = 'CellEditCallback';
                    case {'checkbox' 'listbox' 'popupmenu' 'slider'}
                        ansname = 'Value';
                end
                
                % Set control's properties
                if ~strcmp(fmt.type,'text')
                    set(h,ena{:},fcnname,fcn,ansname,val,aux{:});
                    
                    % Set custom callbacks
                    if ~isempty(fmt.callback)
                        for n = 1:numel(cbnames)
                            % set if the specfied callback function not already assigned
                            % as the part of INPUTSDLG functionality.
                            if ~strcmpi(cbnames{n},fcnname) || isempty(fcn)
                                cbfcn = @(hobj,evt)fmt.callback.(cbnames{n})(hobj,evt,handles.ctrls,k);
                                try
                                    set(h,cbnames{n},cbfcn);
                                catch
                                    error('Invalid callback function name.');
                                end
                            end
                        end
                    end
                end
            end
            
            % make sure we are on screen
            movegui(handles.fig)
            
            % if there is a figure out there and it's modal, we need to be modal too
            if ~isempty(gcbf) && strcmp(get(gcbf,'WindowStyle'),'modal')
                set(handles.fig,'WindowStyle','modal');
            end
            
            % if user createfcn defined, call it now
            if ~isempty(Options.CreateFcn)
                Options.CreateFcn(handles.fig,[],handles.ctrls);
            end
            
            set(handles.fig,'Visible','on');
            drawnow;
            
            % set focus on the first uicontol
            h = findobj(handles.ctrls(:,1),'flat','-not','Style','text');
            if ~isempty(h)
                h = h(1);
                switch get(h,'type')
                    case 'uicontrol', uicontrol(h);
                    case 'uipanel', uicontrol(get(h,'SelectedObject'));
                    case 'uitable', uitable(h);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function answer = getAnswer()
            
            answer = cell(size(Answer));
            
            % retrieve answer from controls
            for i = 1:numel(answer)
                
                h = handles.ctrls(i,1);
                fmt = Formats(i);
                
                switch fmt.style
                    case {'checkbox' 'popupmenu' 'listbox' 'slider'}
                        val = get(h,'Value');
                        
                        if ~isempty(val)
                            switch fmt.format
                                case 'text' % 'popupmenu' 'listbox'
                                    str = get(h,'String');
                                    if isempty(str) == 1;
                                        answer{i} = {};
                                    else
                                        answer{i} = str(val);
                                    end
                                case 'logical' % 'checkbox'
                                    answer{i} = val==get(h,'Max');
                                otherwise %{'float' 'integer'}
                                    answer{i} = val;
                            end
                        end
                    case 'edit'
                        str = get(h,'String');
                        switch fmt.format
                            case {'float','integer'}
                                if ~isempty(str)
                                    answer{i} = str2double(str);
                                else
                                    answer{i} = [];
                                end
                            case 'date'
                                if isempty(str)
                                    answer{i} = [];  % Return an empty date vector if no date was entered
                                else
                                    answer{i} = datevec(str, fmt.limits);
                                end
                            case 'file'
                                if diff(fmt.limits)>1
                                    answer{i} = cellstr(str);
                                else
                                    answer{i} = str;
                                end
                            case {'vector'}
                                answer{i} = str2num(str); %#ok
                            otherwise %case {'text' 'dir'}
                                answer{i} = str;
                        end
                    case {'radiobutton' 'togglebutton'} % uibuttongroup
                        hbtn = get(h,'SelectedObject');
                        if isempty(hbtn) && strcmp(fmt.required,'on')
                            disp('required but not given')
                        end
                        if strcmp(fmt.format,'text')
                            if isempty(hbtn)
                                answer{i} = '';
                            else
                                answer{i} = get(hbtn,'String');
                            end
                        else
                            if isempty(hbtn)
                                answer{i} = [];
                            else
                                answer{i} = find(hbtn==get(h,'UserData'));
                            end
                        end
                    case 'pushbutton'
                        if strcmp(fmt.type,'color')
                            answer{i} = get(h,'BackgroundColor');
                            if strcmp(fmt.format,'integer')
                                answer{i} = uint8(round(answer{i}*255));
                            end
                        else
                            answer{i} = get(h,'UserData');
                        end
                    case 'table'
                        answer{i} = get(h,'Data');
                end
            end
        end
        
        function reqmet = checkReq(answer)
            idx = cellfun(@isempty,answer);
            reqmet = ~any(IsRequired(:)&idx(:));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Control Button Callback callback functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function doFigureKeyPress(obj, evd)
            
            [tf,I] = ismember(evd.Key,{'return','space','escape'});
            if ~tf, return; end % nothing special to do
            
            if any(I==[1 2])
                
                % Ignore under a condition when GCO is an edit uicontrol
                % Known Potential Issue: GCO could be modified from outside
                hedit = findobj(gco,'flat','type','uicontrol','style','edit');
                
                if ~isempty(hedit)
                    
                    % check for the conditions
                    if I==1 % return
                        % resume only if not currently in a multi-line edit uicontrol
                        if get(hedit,'Max')-get(hedit,'Min')>1
                            return;
                        end
                    else % space
                        % don't resume if current focus is on an edit control
                        return;
                    end
                end
                
                % equivalent to pressing OK button
                set(obj,'UserData','OK');
            elseif I==3 % Cancel
                % equivalent to pressing Cancel button
                set(obj,'UserData','Cancel');
            else
                return;
            end
            
            % set focus on OK button (to update the value of the current control)
            uicontrol(handles.btns(1));
            
            % if reached this far, valid key press to close the dialog
            uiresume(obj);
            
        end
        
        function doOKCancel(~, ~, isok)
            if isok
                set(gcbf,'UserData','OK');
                uiresume(gcbf);
            else
                set(gcbf,'UserData','Cancel');
                uiresume(gcbf); % cancel
            end
        end
        
        function doApply(~,~)
            DefaultAnswer = getAnswer(); % retrieve the current answers from the controls
            DefaultReqMet = checkReq(DefaultAnswer);
            Applied = true;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UICONTROL ButtonDownFcn callback functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function checkNumericRange(hObj,evt,k,cbfcn)
            
            fmt = Formats(k);
            
            str = get(hObj,'String');
            if isempty(str)
                if fmt.required
                    % show an error dialog
                    h = errordlg('This parameter must be filled with a value.','Required Value','modal');
                    uiwait(h);
                    % revert back to the previous value
                    set(hObj,'String',num2str(Answer{k}));
                    return;
                else
                    Answer{k} = [];
                end
            else
                
                % convert to float/integer
                val = str2double(str);
                isint = strcmp(fmt.format, 'integer');
                if isint
                    val = round(val);  % Round to the nearest integer
                end
                
                lim = fmt.limits;
                if val>=lim(1) && val<=lim(2) % good result
                    % Re-format the control's text according to the value
                    set(hObj, 'String', num2str(val));
                    Answer{k} = val; % store the numeric answer to revert to later
                else % out-of-range value
                    if isint
                        msg = sprintf('%d, %d',lim(1),lim(2));
                    else
                        msg = sprintf('%f, %f',lim(1),lim(2));
                    end
                    h = errordlg(sprintf('This parameter must be within the range [%s].',msg),'Invalid Value','modal');
                    uiwait(h);
                    
                    % revert back to the previous value
                    set(hObj,'String',num2str(Answer{k}));
                    return;
                end
            end
            
            % run custom callback function
            if ~isempty(cbfcn)
                cbfcn(hObj,evt,handles.ctrls,k);
            end
            
        end
        
        function checkVector(hObj,evd,k,cbfcn)
            
            fmt = Formats(k);
            
            str = get(hObj,'String');
            if isempty(str)
                if fmt.required
                    % show an error dialog
                    h = errordlg('This parameter must be filled with a value.','Required Value','modal');
                    uiwait(h);
                    % revert back to the previous value
                    set(hObj,'String',num2str(Answer{k}));
                    return;
                else
                    Answer{k} = [];
                    return;
                end
            end
            
            % convert to float/integer
            val = str2num(str); %#ok
            N = numel(val);
            
            lim = fmt.limits;
            if N>=lim(1) && N<=lim(2) % good result
                % Re-format the control's text according to the value
                Answer{k} = val(:); % force to be a column vector and store to revert to later
                set(hObj, 'String', num2str(Answer{k}));
            else % out-of-range value
                h = errordlg(sprintf('This vector parameter must have %d to %d elements.',lim(1),lim(2)),'Invalid Number of Elements','modal');
                uiwait(h);
                
                % revert back to the previous value
                set(hObj,'String',num2str(Answer{k}));
            end
            
            % run custom callback function
            if ~isempty(cbfcn)
                cbfcn(hObj,evd,handles.ctrls,k);
            end
        end
        
        function checkDate(hObj,evd,k,cbfcn)
            
            format = Formats(k).limits;
            
            str = get(hObj,'string');
            if isempty(str)  % Avoid calling datenum() which prints a warning for empty strings
                Answer{k} = '';
                set(hObj, 'String', Answer{k});
                return;
            end
            try
                num = datenum(str, format);  % Check if the input matches the custom date format first
            catch
                try
                    num = datenum(str);  % Check if the input matches any other supported date format
                catch
                    h = errordlg(sprintf('Unsupported date format.'),'Invalid Value','modal');
                    uiwait(h);
                    set(hObj,'String',Answer{k});
                    return;
                end
            end
            Answer{k} = datestr(num,format);
            set(hObj,'String',Answer{k});
            
            % run custom callback function
            if ~isempty(cbfcn)
                cbfcn(hObj,evd,handles.ctrls,k);
            end
        end
        
        function openFilePrompt(hObj,evd,k,cbfcn)
            fmt = Formats(k);
            spec = fmt.items;
            opt = {};
            mode = diff(fmt.limits);
            filename = get(hObj,'String');
            if mode<=0 % uiputfile
                uifilefcn = @uiputfile;
                file = filename;
            else
                uifilefcn = @uigetfile;
                if mode>1 % multi-select
                    file = get(hObj,'UserData');
                    opt = {'MultiSelect','on'};
                else
                    file = filename;
                end
            end
            
            % open the file prompt
            [f,p] = uifilefcn(spec,'',file,opt{:});
            if ~(isequal(f,0) && isequal(p,0)) % canceled, no change
                
                % store & display the data
                file = fullfile(p,f); % form full path(es) to the selected files
                set(hObj,'String',file);
                if mode>1
                    set(hObj,'UserData',p);
                end
                
                % run custom callback function
                if ~isempty(cbfcn)
                    cbfcn(hObj,evd,handles.ctrls,k);
                end
            end
            
            % bring focus back to the control
            uicontrol(hObj);
        end
        
        function openDirPrompt(hObj,evd,k,cbfcn)
            
            p = uigetdir(get(hObj,'String'));
            if ~isequal(p,0)
                % display the selected directory
                set(hObj,'String',p);
                
                % run custom callback function
                if ~isempty(cbfcn)
                    cbfcn(hObj,evd,handles.ctrls,k);
                end
            end
            
            % bring focus back to the control
            uicontrol(hObj);
        end
        
        function openColorPrompt(hObj,evd,k,cbfcn)
            
            p = uisetcolor(get(hObj,'BackgroundColor'));
            if ~isempty(p)
                % display the selected color
                set(hObj,'BackgroundColor',p);
                
                % run custom callback function
                if ~isempty(cbfcn)
                    cbfcn(hObj,evd,handles.ctrls,k);
                end
            end
            
            % bring focus back to the control
            uicontrol(hObj);
        end
    end

    function S = selectivecell2struct(C,fields)
        idx = ~cellfun(@isempty,fields);
        if any(idx)
            idx = find(idx)';
            S = cell2struct(C(idx,:),fields(idx),1);
        else
            S = C;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKPROMPT :: Check Prompt input is valid & fill default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Prompt,FieldNames,Units,TooltipStrings,err] = checkprompt(Prompt)
        
        % default configuration
        FieldNames = {}; % answer in a cell
        Units = {}; % no units
        TooltipStrings = {};
        
        % standard error
        err = {'inputsdlg:InvalidInput','Prompt must be a cell string with up to four columns.'};
        
        if isempty(Prompt), Prompt = {'Input:'};
        elseif ~iscell(Prompt), Prompt = cellstr(Prompt);
        end
        
        [nrow,ncol] = size(Prompt);
        
        % prompt given in a row -> transpose
        if ncol>4
            if nrow<4, Prompt = Prompt.'; [nrow,ncol] = size(Prompt);
            else return; % too many columns given
            end
        end
        
        % struct fields defined
        if ncol>1
            idx = cellfun(@isempty,Prompt(:,2));
            FieldNames = Prompt(:,2);
            FieldNames(idx) = {''}; % make sure it is empty cellstr
            
            idx(:) = ~idx;
            if numel(unique(FieldNames(idx)))~=sum(idx)
                err{2} = 'Duplicate struct field name found.';
                return;
            end
        else
            FieldNames = repmat({''},nrow,1);
        end
        
        % unit labels defined
        if ncol>2, Units = Prompt(:,3);
        else       Units = repmat({''},nrow,1);
        end
        
        % tooltip strings defined
        if ncol>3, TooltipStrings = Prompt(:,3);
        else       TooltipStrings = repmat({''},nrow,1);
        end
        
        % return only the labels in Prompt argument
        Prompt(:,2:end) = [];
        
        err = {}; % all cleared
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKFORMATS :: Check Formats input is valid & fill default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Formats,err] = checkformats(Formats,NumQuest,fields)
        
        err = {};
        
        if isempty(Formats) % if Formats not defined, use the first entry (edit:text)
            Formats = fields(1:2,:);
            Formats(cellfun(@(c)iscell(c)&&isempty(c),Formats)) = {{{}}};
            Formats = repmat(struct(Formats{:}),NumQuest,1);
        end
        
        fnames = lower(fields(1,:));
        fields(1,:) = [];
        nfields = numel(fnames); % sans the first row
        
        % backward compatibility (NumLines)
        if isnumeric(Formats)
            [rw,cl]=size(Formats);
            ok = rw==1;
            if ok
                OneVect = ones(NumQuest,1);
                if cl == 2, NumLines=Formats(OneVect,:);
                elseif cl == 1, NumLines=Formats(OneVect);
                elseif cl == NumQuest, NumLines = Formats';
                else ok = false;
                end
            end
            if rw == NumQuest && any(cl == [1 2]), NumLines = Formats;
            elseif ~ok
                err = {'MATLAB:inputdlg:IncorrectSize', 'NumLines size is incorrect.'};
                return;
            end
            
            % set to default edit control (column stacked)
            fields(3:end,:) = []; % all to be edit boxes
            Formats = repmat(struct(fields{:}),NumQuest,1);
            
            % set limits according to NumLines(:,1)
            numlines = mat2cell([zeros(NumQuest,1) NumLines(:,1)],ones(NumQuest,1),2);
            [Formats.limits] = deal(numlines{:});
            
            % sets the width to be 10*NumLines(:,2)
            if (size(NumLines,2) == 2)
                sizes = mat2cell([zeros(NumQuest,1) NumLines(:,2)],ones(NumQuest,1),2);
                [Formats.size] = deal(sizes{:});
            end
            
            return;
        elseif ischar(Formats) || iscellstr(Formats) % given type
            if ischar(Formats), Formats = cellstr(Formats); end
            Formats = cell2struct(Formats,'type',3);
        elseif ~isstruct(Formats)
            err = {'inputsdlg:InvalidInput','FORMATS must be an array of structure.'};
            return
        end
        
        % Dialog grid dimension
        fdims = size(Formats);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % If span field is given, fill Format struct accordingly
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(Formats,'span')
            
            idx = ~arrayfun(@(s)isempty(s.span),Formats);
            if ~arrayfun(@(s,sz)isnumeric(s.span) && ...
                    (numel(s.span)==2 && all(s.span==floor(s.span) & s.span>0)),Formats(idx))
                err = {'inputsdlg:InvalidInput','FORMATS.span must be 2-element vectors of positive integers.'};
                return;
            end
            if ~all(arrayfun(@(s)all(s.span<=fdims),Formats(idx)))
                err = {'inputsdlg:InvalidInput','FORMATS.span extends the control out of the grid size specified by Formats matrix.'};
                return;
            end
            
            [i,j] = find(arrayfun(@(s)any(size(s.span)>[1 1]),Formats));
            for ind = [i j].' % for each spanning control
                span = Formats(ind(1),ind(2)).span;
                % extend from left on the top-most row
                for col = ind(2)+(1:span(2)-1)
                    Formats(ind(1),col).type = 'none';
                    Formats(ind(1),col).limits = [0 1]; % extend from left
                end
                % extend from above for all other rows
                for row = ind(1)+(1:span(1)-1)
                    for col = ind(2)+(0:span(2)-1)
                        Formats(row,col).type = 'none';
                        Formats(row,col).limits = [1 0]; % extend from above
                    end
                end
            end
            Formats = rmfield(Formats,'span');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Convert Formats to cell
        [~,I] = ismember(lower(fieldnames(Formats)),fnames);
        if any(I==0)
            err = {'inputsdlg:InvalidFormatsField','FORMATS contains invalid field name(s).'};
            return
        end
        fvals = cell([prod(fdims) nfields]);
        fvals(:,I) = struct2cell(Formats(:)).';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check type field (Column 1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if type field does not contain NumQuest many non-'none' types,
        % specify the first empty fields to be the default type
        Iempty = cellfun(@isempty,fvals(:,1));
        Inone = cellfun(@(s)strcmp(s,'none'),fvals(:,1));
        N = NumQuest-sum(~(Iempty|Inone));
        if N<0
            err = {'inputsdlg:InvalidInput','Not enough FORMATS size for the display.'};
            return
        elseif N>0
            idx = find(Iempty,N);
            fvals(idx) = fields(2,1);
            Iempty(idx) = false;
        end
        
        % set the rest of empty entries to none (spacer)
        fvals(Iempty,1) = {'none'};
        
        % identify unknown types
        fvals(:,1) = lower(fvals(:,1));
        if verLessThan('matlab','7.14') % pre-R2012a
            [tf,typeindex] = ismember(fvals(:,1),fields(end:-1:1,1)); % typeindex: first matching case in fields
            typeindex(tf) = size(fields,1)-typeindex(tf)+1;
        else
            [~,typeindex] = ismember(fvals(:,1),fields(:,1),'R2012a'); % typeindex: first matching case in fields
        end
        
        % check for invalid type field contents
        if any(typeindex==0)
            err = {'inputsdlg:InvalidInput','FORMATS.type must be one of {''check'',''edit'',''list'',''range'',''color'',''button'',''table'',''none''}.'};
            return
        end
        
        % check number of entries matching NumQuest (number of PROMPT elements)
        if sum(~strcmp('none',fvals(:,1)))~=NumQuest
            err = {'inputsdlg:InvalidInput',sprintf('%s\n%s',...
                'FORMATS must have matching number of elements to PROMPT (exluding ''none'' type).',...
                'If .span field is used, also check for overlapping controls.')};
            return
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fill empty format field (Column 2) to their type's defaults
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx = cellfun(@isempty,fvals(:,2));
        fvals(idx,2) = fields(typeindex(idx),2);
        
        Itable = strcmp('table',fvals(:,1));
        Iother = ~Itable;
        
        if ~all(iscellstr(fvals(Iother,2)))
            err = {'inputsdlg:InvalidInput','FORMATS.format must be given as a cell string for a non-table Type.'};
            return;
        end
        if ~(all(iscell(fvals(Itable,2))))
            err = {'inputsdlg:InvalidInput','FORMATS.format must be given as a cell array, specifying ColumnFormat of uitable.'};
            return;
        end
        
        fvals(Iother,2) = lower(fvals(Iother,2)); % make format strings lower case
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set all empty fields to type's defaults
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %[~,locb] = ismember(fvals(Iother,[1 2]),fields(1:end-1,[1 2]),'rows');
        a = fvals(Iother,[1 2]);
        if ~isempty(a)
            uB = fields(1:end-1,[1 2]); % guarantees to be unique, excludes 'table' entry
            [sortA,IndSortA] = sortrows(a);
            groupsSortA = [true;any(~strcmp(sortA(1:end-1,:),sortA(2:end,:)),2)]; % Finds the index of matching entriesuA = sortA(groupsSortA,:);
            uA = sortA(groupsSortA,:);
            icA = cumsum(groupsSortA);                             % Lists position, starting at 1.
            icA(IndSortA) = icA;                                  % Re-reference indC to indexing of sortA.
            [sortuAuB,IndSortuAuB] = sortrows([uA;uB]);
            d = all(strcmp(sortuAuB(1:end-1,:),sortuAuB(2:end,:)),2);     % d indicates matching entries
            ndx1 = IndSortuAuB(d);                          % NDX1 are locations of repeats in C
            szuA = size(uA,1);
            [lia,locb] = ismember(icA,ndx1,'R2012a');    % Find locb by using given indices
            d = find(d);
            newd = d(locb(lia));                    % NEWD is D for non-unique A
            where = IndSortuAuB(newd+1)-szuA;   % Index values of uB through UNIQUE
            locb(lia) = where;                      % Return first or last occurrence of A within B
        else
            locb = [];
        end
        
        if ~all(lia)
            err = {'inputsdlg:InvalidInput','Invalid FORMATS.format specified.'};
            return;
        end
        
        typefmtindex = zeros(size(typeindex));
        typefmtindex(Iother) = locb;
        typefmtindex(Itable) = size(fields,1);
        
        idx = cellfun(@isempty,fvals);
        [I,J] = find(idx);
        T = typefmtindex(I);
        
        fvals(idx) = fields(sub2ind(size(fields),T,J));
        % Set all characters in fixed-format string field to lower case
        fvals(:,[3 7]) = lower(fvals(:,[3 7])); % style & enable
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check for proper combination of format specs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [types,~,IB] = unique(fvals(:,1));
        for n = 1:numel(types)
            
            idx = find(IB==n); % index of the matching grid elements
            
            switch types{n}
                case 'text'
                    % check style
                    if ~all(strcmp(fvals(idx,3),'text'))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''text'' type must be ''text''.'};
                        return
                    end
                case 'check'
                    % check format
                    if ~any(ismember(fvals(idx,2),{'logical','integer'}))
                        err = {'inputsdlg:InvalidInput','FORMATS.format for ''check'' type must be ''logical'' or ''integer''.'};
                        return
                    end
                    % check style
                    if ~all(strcmp(fvals(idx,3),'checkbox'))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''check'' type must be ''checkbox''.'};
                        return
                    end
                case 'edit'
                    % check format
                    fmts = fvals(idx,2);
                    if ~any(ismember(fmts,{'text','data','float','integer','file','dir','vector'}))
                        err = {'inputsdlg:InvalidInput','FORMATS.format for ''edit'' type must be one of ''text'', ''data'', ''float'', ''integer'', ''file'', ''dir'', or ''vector''.'};
                        return
                    end
                    % check style
                    if ~any(strcmp(fvals(idx,3),'edit'))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''edit'' type must be ''edit''.'};
                        return
                    end
                case 'list'
                    % check format
                    if ~any(ismember(fvals(idx,2),{'integer','text'}))
                        err = {'inputsdlg:InvalidInput','FORMATS.format for ''list'' type must be ''integer'' or ''text''.'};
                        return
                    end
                    
                    % check style
                    if ~any(ismember(fvals(idx,3),{'listbox','popupmenu','radiobutton','togglebutton'}))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''list'' type must be one of ''listbox'', ''popupmenu'', ''radiobutton'', or ''togglebutton''.'};
                        return
                    end
                    
                    % check items - convert if string array or numeric array given
                    idx1 = idx(cellfun(@ischar,fvals(idx,4)));
                    fvals(idx1,4) = cellfun(@cellstr,fvals(idx1,4),'UniformOutput',false);
                    idx2 = idx(cellfun(@isnumeric,fvals(idx,4)));
                    fvals(idx2,4) = cellfun(@num2cell,fvals(idx2,4),'UniformOutput',false);
                    
                    if ~any(cellfun(@iscellstr,fvals(idx,4)))
                        err = {'inputsdlg:InvalidInput','FORMATS.items must be either a cell of strings or of numbers.'};
                        return
                    end
                case 'range'
                    % check format
                    if ~any(strcmp(fvals(idx,2),'float'))
                        err = {'inputsdlg:InvalidInput','FORMATS.format for ''range'' type must be ''float''.'};
                        return
                    end
                    % check style
                    if ~any(strcmp(fvals(idx,3),'slider'))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''range'' type must be ''slider''.'};
                        return
                    end
                case 'table'
                    if ~any(strcmp(fvals(idx,3),'table'))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''table'' type must be ''table''.'};
                        return
                    end
                case 'color'
                    if ~any(ismember(fvals(idx,2),{'float','integer'}))
                        err = {'inputsdlg:InvalidInput','FORMATS.format for ''range'' type must be ''float'' or ''integer''.'};
                        return
                    end
                    if ~any(strcmp(fvals(idx,3),'pushbutton'))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''color'' type must be ''pushbutton''.'};
                        return
                    end
                case 'button'
                    if ~any(strcmp(fvals(idx,3),'pushbutton'))
                        err = {'inputsdlg:InvalidInput','FORMATS.style for ''range'' type must be ''slider''.'};
                        return
                    end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check limits (Column 5)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % edit::date - limits specifies the date display format
        idx = strcmp(fvals(:,1),'edit') & strcmp(fvals(:,2),'date');
        if any(idx)
            lims = fvals(idx,5);
            idx1 = cellfun(@(l)ischar(l) && size(l,1)==1,lims);
            try
                cellfun(@(l)datestr(1,l),lims(idx1),'UniformOutput',false);
            catch
                err = {'inputsdlg:InvalidInput', 'Invalid free-form format string in FORMATS.limits for ''date'' control.'};
                return;
            end
            if any(cellfun(@(l)isnumeric(l)&&isscalar(l)&&any(l==[0 1 2 6 13 14 15 16 23]),lims(idx1)))
                err = {'inputsdlg:InvalidInput','FORMATS.limits for ''edit::date'' format must be one of 0,1,2,6,13,14,15,16,23.'};
                return;
            end
        end
        
        % table - limits specifies the Column widths
        I = idx;
        idx = strcmp(fvals(:,1),'table');
        if any(idx)
            lim = fvals(idx,5);
            
            if ~any(cellfun(@iscell,lim))
                err = {'inputsdlg:InvalidInput','FORMATS.limits for ''table'' type must be given as a cell vector.'};
                return;
            end
        end
        
        % all other contols - limits must be 2-element vector
        I(:) = ~(I | idx);
        if ~all(cellfun(@(l)isnumeric(l)&&numel(l)==2,fvals(I,5)))
            err = {'inputsdlg:InvalidInput','FORMATS.limits must be given as a two-element vector.'};
            return;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check size (Column 6)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~all(cellfun(@(sz)isnumeric(sz) && any(numel(sz)==[1 2]) && ~any(isnan(sz)),fvals(:,6)))
            err = {'inputsdlg:InvalidInput','FORMATS.size must be 1 or 2 element non-NaN vector.'};
            return
        end
        idx = cellfun(@numel,fvals(:,6))==1;
        fvals(idx,6) = cellfun(@(sz)[sz 0],fvals(idx,6),'UniformOutput',false);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check enable (Column 7)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~all(ismember(fvals(:,7),{'on','inactive','off'}))
            err = {'inputsdlg:InvalidInput','FORMATS.enable must be one of {''on'',''inactive'',''off''}.'};
            return;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check required (Column 8)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~all(ismember(fvals(:,8),{'on','off'}))
            err = {'inputsdlg:InvalidInput','FORMATS.required must be ''on'' or ''off''.'};
            return;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check callback (Column 9)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx0 = cellfun(@isstruct,fvals(:,9));
        idx1 = cellfun(@(cb)isa(cb,'function_handle'),fvals(:,9));
        if ~all(idx0|idx1|cellfun(@isempty,fvals(:,9)))
            err = {'inputsdlg:InvalidInput','FORMATS.callback must be given as a function handle or a struct containing function handles.'};
            return;
        end
        if any(idx0) % given as struct, make sure all its elements are function handles
            I = find(idx0);
            for n = 1:numel(I)
                cbvals = fvals{idx0(I(n)),9};
                if ~(isempty(cbvals) || all(structfun(@(cb)isa(cb,'function_handle'),cbvals)))
                    err = {'inputsdlg:InvalidInput','FORMATS.callback must be given as a function handle or a struct containing function handles.'};
                    return;
                end
            end
        end
        if any(idx1) % convert to the struct form
            idx1 = find(idx1);
            types = fvals(idx1,1);
            formats = fvals(idx1,2);
            styles = fvals(idx1,3);
            cbs = fvals(idx1,9);
            
            Jtb = cellfun(@(t)strcmp('table',t),types); % uitable
            if any(Jtb)
                fvals(idx1(Jtb),9) = cellfun(@(cb)struct('CellEditCallback',cb),cbs(Jtb),'UniformOutput',false);
            end
            
            Jbg = cellfun(@(t,s)strcmp('list',t) & ~ismember(s,{'listbox','popupmenu'}),types,styles); % uibuttongroup
            if any(Jbg)
                fvals(idx1(Jbg),9) = cellfun(@(cb)struct('SelectionChangeFcn',cb),cbs(Jbg),'UniformOutput',false);
            end
            
            Jbd = cellfun(@(s,t,f)strcmp('edit',s) && (strcmp('color',t) || ismember(f,{'file','dir','font'})),...
                styles,types,formats); % inactive edit uicontrol
            if any(Jbd)
                fvals(idx1(Jbd),9) = cellfun(@(cb)struct('ButtonDownFcn',cb),cbs(Jbd),'UniformOutput',false);
            end
            
            J = ~(Jtb|Jbg|Jbd); % uicontrol
            if any(J)
                fvals(idx1(J),9) = cellfun(@(cb)struct('Callback',cb),cbs(J),'UniformOutput',false);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check prompt label location (Column 10)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~all(ismember(fvals(:,10),{'lefttop','leftmiddle','leftbottom','topleft','topcenter','topright'}))
            err = {'inputsdlg:InvalidInput','FORMATS.labelloc must be ''lefttop'', ''leftmiddle'', ''leftbottom'', ''topleft'', ''topcenter'', or ''topright''.'};
            return;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check units label location (Column 11)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~all(ismember(fvals(:,11),{'righttop','rightmiddle','rightbottom','bottomleft','bottomcenter','bottomright'}))
            err = {'inputsdlg:InvalidInput','FORMATS.unitsloc must be ''righttop'', ''rightmiddle'', ''rightbottom'', ''bottomleft'', ''bottomcenter'', or ''bottomright''.'};
            return;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check label margins (Column 12)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~all(cellfun(@(v)isnumeric(v) && any(numel(v)==[1 2]) && ~any(isinf(v)|v<0),fvals(:,12)))
            err = {'inputsdlg:InvalidInput','FORMATS.margin must be 1 or 2 element positive vector.'};
            return
        end
        
        idx = cellfun(@isscalar,fvals(:,12));
        fvals(idx,12) = cellfun(@(v)repmat(v,[1 2]),fvals(idx,12),'UniformOutput',false);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather back as Formats struct
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Formats = reshape(cell2struct(fvals,fnames,2),fdims);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKDEFAULTS :: Check the specified default values are compatible
%%% with Formats and if one not given fill in an initial value
    function [DefAns,TileMode,err] = checkdefaults(DefAns,Formats,FieldNames)
        % AnsStr:
        
        % set the TileMode
        TileMode = [1 1]; % default: no tiling
        dims = size(Formats); % grid dimension
        if any(dims==1) % To tile, formats must be given as a vector
            tiledim = (dims(2)==1) + 1; % if 1, tile rows; if 2, tile columns
            ansdims = size(DefAns);
            tf = ansdims>1;
            if isstruct(DefAns) && any(tf) % if vector, tile
                TileMode(tiledim) = numel(DefAns);
            elseif iscell(DefAns) && all(tf) % if matrix, tile
                TileMode(tiledim) = ansdims(2);
            end
        end
        
        if all(TileMode==1) % no tiling
            DefAns = DefAns(:);
        end
        
        % reshape to a "row" vector & trim Formats to only include relevant entries (non-'none' types)
        Formats = Formats'; % go through row first
        Formats = Formats(~strcmp('none',{Formats.type}));
        len = numel(Formats); % if tiled, expects len*prod(TileMode) DefAns
        
        if isempty(DefAns) % if DefAns not given
            DefAns = cell(len,1); % will set DefAns to default values
        elseif isstruct(DefAns)
            if isempty(FieldNames) % FieldNames must be given via PROMPT
                err = {'inputsdlg:InvalidInput','DEFAULTANSWER is given in as a struct but its field names are not specified in PROMPT (the 2nd column).'};
                return;
            end
            
            % Convert Struct to cell according to FieldNames
            Inotempty = find(~cellfun(@isempty,FieldNames)); % ignore prompt w/o FieldName
            [~,I] = ismember(fieldnames(DefAns),FieldNames(Inotempty));
            if any(I==0)
                err = {'inputsdlg:InvalidFormatsField','DEFAULTANSWER contains invalid field name(s).'};
                return
            end
            DefStr = DefAns; % save the original input
            DefAns = cell(len,prod(TileMode)); % len==numel(FieldNames) is guaranteed
            DefAns(Inotempty(I),:) = struct2cell(DefStr(:));
            
        elseif ~iscell(DefAns)
            err = {'inputsdlg:InvalidInput','Default answer must be given in a cell array or a structure.'};
            return;
        elseif length(DefAns)~=len
            err = {'inputsdlg:InvalidInput','Default answer cell dimension disagrees with the number of prompt'};
            return;
        end
        
        % go through each default values
        for k = 1:len
            
            fmt = Formats(k);
            
            % if any of the tiled element is given empty & its answer is required
            Iempty = cellfun(@isempty,DefAns(k,:));
            if any(Iempty)
                
                switch fmt.type
                    case 'check' % off
                        if strcmp(fmt.format,'logical')
                            val = false;
                        else
                            val = fmt.limits(1);
                        end
                    case 'edit'
                        switch fmt.format
                            case {'float','integer'}
                                liminf = isinf(fmt.limits);
                                if all(liminf) % both limits inf
                                    val = 0;
                                elseif any(liminf) % 1 limit inf
                                    val = fmt.limits(~liminf);
                                else % neither is inf
                                    val = round(mean(fmt.limits));
                                end
                            case 'vector'
                                val = zeros(1,fmt.limits(1));
                            otherwise %{'text','date','file','dir'}
                                val = '';
                        end
                    case 'list' % first item
                        val = 1;
                    case 'range' % middle value
                        val = mean(fmt.limits);
                    case 'color'
                        val = [1 1 1];
                    case 'table'
                        val = {};
                    otherwise
                        val = [];
                end
                
                % set the default value to all empty controls on all tiles
                [DefAns{k,Iempty}] = deal(val);
            end
            
            % for all entries that are not empty check the validity of given values
            Ifilled = ~Iempty;
            if ~strcmp(fmt.type,'table') && any(Ifilled)
                
                vals = DefAns(k,Ifilled); % given default values
                
                switch fmt.format
                    case 'text'
                        % must be char or cellsctr
                        Jischar = cellfun(@ischar,vals);
                        Jiscellstr = cellfun(@iscellstr,vals);
                        if ~all(Jischar|Jiscellstr)
                            err = {'inputsdlg:InvalidInput','Default text data format must be char.'};
                            return;
                        end
                        
                        % for the list type, value must be one of the allowed item
                        if strcmp(fmt.type,'list')
                            isUiControl = ismember(fmt.style,{'listbox', 'popupmenu'});
                            if any(Jischar)
                                [tf,idx] = ismember(vals,fmt.items);
                                if ~all(tf)
                                    err = {'inputsdlg:InvalidInput','Default list item is not valid.'};
                                    return;
                                end
                                if isUiControl
                                    % convert to integer format
                                    vals(Jischar) = num2cell(idx);
                                end
                            end
                            if any(Jiscellstr)
                                Jiscellstr = find(Jiscellstr);
                                for j = Jiscellstr(:).'
                                    [tf,idx] = ismember(vals{j},fmt.items);
                                    if ~all(tf)
                                        err = {'inputsdlg:InvalidInput','Default list item is not valid.'};
                                        return;
                                    end
                                    % convert to integer format
                                    if isUiControl
                                        vals{j} = idx;
                                    end
                                end
                            end
                            
                            DefAns(k,Ifilled) = vals;
                        end
                        
                    case 'date' % given as a date vector
                        % must be a date string or date number or date vector
                        if ~all(cellfun(@(v)ischar(v)||isnumeric(v),vals))
                            err = {'inputsdlg:InvalidInput','Default date date must be a valid datenum, datestr, or datevec.'};
                            return;
                        end
                        % store the date as a string
                        try
                            DefAns(k,Ifilled) = cellfun(@(v)datestr(v, fmt.limits),vals,'UniformOutput',false);
                        catch % not a valid date input given
                            err = {'inputsdlg:InvalidInput','Default date date must be a valid datenum, datestr, or datevec.'};
                            return;
                        end
                    case 'float'
                        
                        if ~all(cellfun(@(v)isfloat(v)&&all(v>=fmt.limits(1)&v<=fmt.limits(2)),vals))
                            err = {'inputsdlg:InvalidInput','Default float data must be a numeric within specified limits.'};
                            return;
                        end
                        
                        if strcmp(fmt.type,'color') % must be a valid RGB tuple
                            if ~all(cellfun(@(v)numel(v)==3 && max(v)<=1 && min(v)>=0,vals))
                                err = {'inputsdlg:InvalidInput','Default coloar data must have 3 elements.'};
                                return;
                            end
                        else % for all other types, value must be scalar
                            if ~all(cellfun(@isscalar,vals))
                                err = {'inputsdlg:InvalidInput','Default float data must be scalar.'};
                                return;
                            end
                        end
                        
                    case 'integer' % can be multi-select if type=list
                        if ~all(cellfun(@(v)isnumeric(v)&&all(round(v)==v),vals))
                            err = {'inputsdlg:InvalidInput','Default integer data must integer.'};
                            return;
                        end
                        
                        switch fmt.type
                            case 'list'
                                isUiBtnGrp = ismember(fmt.style,{'togglebutton', 'radiobutton'});
                                msel = strcmp(fmt.style,'listbox') && diff(fmt.limits)>1;
                                
                                % must be a valid index to items and
                                % if multiple-selection is not enabled, must be scalar
                                if any(cellfun(@(v)any(v<1|v>numel(fmt.items)),vals)) ...
                                        || ~(msel || all(cellfun(@isscalar,vals)))
                                    err = {'inputsdlg:InvalidInput','Default list index data is invalid.'};
                                    return;
                                end
                                
                                if msel % if multiple-selection enabled, make sure values are unique
                                    DefAns(k,Ifilled) = cellfun(@unique,vals,'UniformOutput',false);
                                elseif isUiBtnGrp
                                    DefAns(k,Ifilled) = fmt.items(vals);
                                end
                                
                            case 'color' % must be a RGB tuple
                                if ~all(cellfun(@(v)numel(v)==3 && max(v)<=255 && min(v)>=0,vals))
                                    err = {'inputsdlg:InvalidInput','Default coloar data must have 3 elements.'};
                                    return;
                                end
                                % convert to float representation
                                DefAns(k,Ifilled) = cellfun(@(v)double(v)/255,vals,'UniformOutput',false);
                            case 'check' % must be a scalar and one of limits
                                if ~all(cellfun(@(v)isscalar(v) && any(v==fmt.limits)))
                                    err = {'inputsdlg:InvalidInput','Default integer check data must be a scalar value matching Format.Limits'};
                                    return;
                                end
                            otherwise %must be a scalar
                                % limits specifies the value range
                                if ~all(cellfun(@(v)v>=fmt.limits(1)&v<=fmt.limits(2),vals))
                                    err = {'inputsdlg:InvalidInput','Out-of-range default integer data.'};
                                    return;
                                end
                                
                                if ~all(cellfun(@isscalar,vals))
                                    err = {'inputsdlg:InvalidInput','Default integer data must be scalar.'};
                                    return;
                                end
                        end
                        
                    case 'logical'
                        
                        if ~all(cellfun(@(v)islogical(v)&&isscalar(v),vals))
                            err = {'inputsdlg:InvalidInput','Default logical data must be of logical or numeric scalar.'};
                            return;
                        end
                        
                        % convert to integer format
                        DefAns(k,Ifilled) = cellfun(@(v)fmt.limits(v+1),vals,'UniformOutput',false);
                        
                    case 'file'
                        % must be char or cellstring
                        if ~all(cellfun(@(v)ischar(v)||iscellstr(v),vals))
                            err = {'inputsdlg:InvalidInput','Default file data must be given as char or cellstr.'};
                            return;
                        end
                        
                        % make everything cellstr
                        vals = cellfun(@cellstr,vals,'UniformOutput',false);
                        
                        dlim = diff(fmt.limits);
                        
                        if dlim<=1 && ~all(cellfun(@isscalar,vals)) % single-file control
                            err = {'inputsdlg:InvalidInput','Multiple default files are given for single-file control.'};
                            return;
                        end
                        
                        if dlim>=0 % for uigetfile, either directory name or existing file name
                            if ~any(cellfun(@(v)all(cellfun(@(file)exist(file,'file'),v)),vals))
                                err = {'inputsdlg:InvalidInput','Default file for uigetfile must exist on the computer.'};
                                return;
                            end
                            
                            % resolve full file name, if only name is given for a file
                            % on matlab path
                            changed = false;
                            for n = 1:prod(TileMode)
                                files = cellfun(@which,vals{n},'UniformOutput',false);
                                I = ~cellfun(@isempty,files);
                                changed = changed || any(I);
                                vals{n}(I) = files(I);
                                
                                % multi-select files must be on a same directory
                                if numel(unique(cellfun(@fileparts,vals{n},'UniformOutput',false)))>1
                                    err = {'inputsdlg:InvalidInput','Default files for multi-select control must be from a same directory.'};
                                    return;
                                end
                            end
                            
                            if changed
                                DefAns(k,Ifilled) = vals;
                            end
                        end
                    case 'dir'
                        if ~all(cellfun(@(v)ischar(v)&&isrow(v)&&isdir(v),vals)) % directory must exist
                            err = {'inputsdlg:InvalidInput','Default dir must be a valid path.'};
                            return;
                        end
                    case 'vector'
                        if ~all(cellfun(@(v)isnumeric(v)&&numel(v)>=fmt.limits(1)&&numel(v)<=fmt.limits(2),vals))
                            err = {'inputsdlg:InvalidInput','Default vector data must be numeric and the number of elements must be within the limit.'};
                            return;
                        end
                end
            end
        end
        
        DefAns = DefAns(:);
        err = {}; % all good
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [Options,FormatFields,err] = checkoptions(UserOptions)
        
        err = {'inputsdlg:InvalidInput',''};
        
        % default options
        Fields = {
            'Resize'             'off'
            'WindowStyle'        'normal'
            'Interpreter'        'tex'
            'DialogSizeMinimum'	[inf inf] % automatically set
            'CancelButton'       'on'
            'ApplyButton'        'off'
            'ButtonNames'        {{'OK','Cancel','Apply'}}
            'Sep'                10
            'AlignControls'      'off'
            'FontSize'           get(0,'DefaultUicontrolFontSize')
            'CreateFcn'          {{}}
            'DeleteFcn'          {{}}
            }.';
        
        % default formats, given type & format
        FormatFields = [
            {'type'   'format'  'style'      'items' 'limits'   'size' 'enable' 'required' 'callback' 'labelloc' 'unitsloc' 'margin'}
            {'edit'   'text'    'edit'       {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]} % default if Formats or Formats.type not given
            {'edit'   'integer' 'edit'       {}      [-inf inf] [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'edit'   'float'   'edit'       {}      [-inf inf] [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'edit'   'vector'  'edit'       {}      [0 inf]    [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]}
            {'edit'   'date'    'edit'       {}      2          [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'edit'   'file'  'edit'  {'*.*' 'All Files'} [0 1] [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'edit'   'dir'     'edit'       {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]} % default if Formats or Formats.type not given
            {'check'  'logical' 'checkbox'   {}      [0 1]      [0 0]  'on'     'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'check'  'integer' 'checkbox'   {}      [0 1]      [0 0]  'on'     'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'list'   'integer' 'popupmenu'  {}      [0 1]      [0 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'list'   'text'    'popupmenu'  {}      [0 1]      [0 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'range'  'float'   'slider'     {}      [0 1]      [0 0]  'on'     'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'color'  'float'   'pushbutton' {}      [0 1]      [65 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'color'  'integer' 'pushbutton' {}      [0 255]    [65 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'button' ''        'pushbutton' {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
            {'text'   ''        'text'       {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]}
            {'none'   ''        ''           {}      [0 0]      [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]} % default if Formats.type empty
            {'table'  {}        'table'      {}      {'auto'}   [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]}];
        
        Options = struct(Fields{:});
        
        if isempty(UserOptions) % no option specified, use default
            err = {};
            return;
        elseif numel(UserOptions)~=1
            err{2} = 'Options struct must be a scalar.';
            return;
        end
        
        % check if User Resize Option is given as on/off string
        if ischar(UserOptions) && strcmpi(UserOptions,{'on','off'})
            UserOptions.Resize = UserOptions;
        elseif ~isstruct(UserOptions)
            err{2} = 'Options must be ''on'', ''off'', or a struct.';
            return;
        end
        
        % to-do: Separate overall options to format options
        % optfnames = fieldnames(UserOptions);
        % optnames = regexpi(optfnames,'default(\S+)','tokens');
        % Iopt = ~cellfun(@isempty,optnames);
        % optnames = lower(cellfun(@(c)c{1},[optnames{Iopt}],'UniformOutput',false));
        % fmtnames = lower(FormatFields(1,:));
        % nfmtfields = numel(fmtnames); % sans the first row
        
        % remove all the format options
        % Options = rmfield(Options,optfnames(Iopt)));
        
        % start with 'type'-'format' combo for the default control type
        % [tf,I] = ismember(optnames,fnames);
        
        
        % check UserOptions struct & update Options fields
        for fname_cstr = fieldnames(UserOptions)' % for each user option field
            
            fname = char(fname_cstr); % use plain char string (not cellstr)
            val = UserOptions.(fname);
            
            % make sure string value is given as a cellstr
            if ischar(val), val = cellstr(val); end
            
            % if field not filled, use default value
            if isempty(UserOptions.(fname)), continue; end
            
            
            switch lower(fname)
                case 'resize'
                    if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
                        err{2} = 'Resize option must be ''on'' or ''off''.';
                        return;
                    end
                    Options.Resize = char(val);
                case 'windowstyle'
                    if numel(val)~=1 || ~any(strcmpi(val,{'normal','modal','docked'}))
                        err{2} = 'WindowStyle option must be ''normal'' or ''modal''.';
                        return;
                    end
                    Options.WindowStyle = char(val);
                case 'interpreter'
                    if numel(val)~=1 || ~any(strcmpi(val,{'latex','tex','none'}))
                        err{2} = 'Interpreter option must be ''latex'', ''tex'', or ''none''.';
                        return;
                    end
                    Options.Interpreter = char(val);
                case 'cancelbutton'
                    if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
                        err{2} = 'CancelButton option must be ''on'' or ''off''.';
                        return;
                    end
                    Options.CancelButton = char(val);
                case 'applybutton'
                    if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
                        err{2} = 'ApplyButton option must be ''on'' or ''off''.';
                        return;
                    end
                    Options.ApplyButton = char(val);
                case 'buttonnames'
                    if ~iscellstr(val)
                        err{2} = 'ButtonNames option must be of cellstr or char type.';
                        return;
                    end
                    
                    % if not all 3 button names are given, use default for unspecified
                    N = numel(val);
                    if (N>3)
                        err{2} = 'ButtonNames option takes up to 3 button names.';
                        return;
                    end
                    Options.ButtonNames(1:N) = val;
                case 'sep'
                    if numel(val)~=1 || ~isnumeric(val) || val<0
                        err{2} = 'Sep option must be non-negative scalar value.';
                        return;
                    end
                    Options.Sep = val;
                case 'aligncontrols'
                    if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
                        err{2} = 'AlignControls option must be ''on'' or ''off''.';
                        return;
                    end
                    Options.AlignControls = char(val);
                case 'unitsmargin'
                    error('UnitMargin option has been deplicated. Use Formats.margin to set individual control''s label margins.');
                case 'fontsize'
                    if ~(isscalar(val) && isnumeric(val) && val>0)
                        err{2} = 'FontSize option must be non-negative scalar value.';
                        return;
                    end
                    Options.FontSize = val;
                case 'dialogsizeminimum'
                    if ~(isnumeric(val) && numel(val)==2 && ~any(isnan(val)))
                        err{2} = 'DialogSizeMinimum option must be 2-element vector.';
                    end
                    % no minimum if not positive
                    idx = val<=0;
                    Options.DialogSizeMinimum(idx) = inf;
                case 'createfcn'
                    if ~(isempty(val) || isa(val,'function_handle'))
                        err{2} = 'CreateFcn option must be a function handle.';
                    end
                    Options.CreateFcn = val;
                case 'deletefcn'
                    if ~(isempty(val) || isa(val,'function_handle'))
                        err{2} = 'DeleteFcn option must be a function handle.';
                    end
                    Options.DeleteFcn = val;
                otherwise
                    warning('inputsdlg:InvalidOption','%s is not a valid option name.',fname);
            end
        end
        
        err = {}; % all cleared
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILDGUI :: Builds the dialog box and returns handles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Formats,handles,sinfo] = buildgui(Title,Prompt,Units,TooltipStrings,FieldNames,Formats,TileMode,Options,Step)
        % 1. Create handle graphic objects for all controls (embedded in uipanels)
        % 2. Generate sinfo to assist object positioning in doResize
        
        sep = Options.Sep;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Tile the controls
        coltile = TileMode(2)>1;
        if any(TileMode>1)
            
            Formats = repmat(Formats,TileMode);
            
            num = numel(Prompt)*prod(TileMode);
            Prompt = reshape(repmat(Prompt,TileMode),num,1);
            FieldNames = reshape(repmat(FieldNames,TileMode),num,1);
            Units = reshape(repmat(Units,TileMode),num,1);
            TooltipStrings = reshape(repmat(TooltipStrings,TileMode),num,1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% determine how to utilize 'none' space
        % Place all the elements at (0,0)
        dim = size(Formats); % display grid dimension
        free = strcmp('none',{Formats.type}); % location of empty block maybe to be occupied by neighbor entry
        num = sum(~free); % number of controls
        map = zeros(dim); % determine which control occupies which block(s)
        order = zeros(1,num); % uicontrol placement order (for tab order)
        n = 1;
        for f = 1:prod(dim)
            if coltile % traverse column-first
                [i,j] = ind2sub(dim,f);
            else % traverse row-first
                [j,i] = ind2sub(dim([2 1]),f);
            end
            m = sub2ind(dim,i,j);
            
            if free(m)
                mode = diff(Formats(m).limits);
                [i,j] = ind2sub(dim,m);
                if mode>0 && j>1, map(m) = map(sub2ind(dim,i,j-1)); % copy from left
                elseif mode<0 && i>1, map(m) = map(sub2ind(dim,i-1,j)); % copy from above
                end % other wise, 0 (nothing occupying)
            else
                map(m) = n;
                order(n) = m;
                n = n + 1;
            end
        end
        
        % remove none's from Formats and order the rest in Prompt order
        Formats = Formats(order).'; % removes all none-types
        
        FigColor=get(0,'DefaultUicontrolBackgroundcolor');
        
        fig = dialog(           ...
            'Visible'     ,'off'   , ...
            'Name'       ,Title   , ...
            'Pointer'     ,'arrow'  , ...
            'Units'      ,'pixels'  , ...
            'UserData'     ,'Cancel'  , ...
            'Tag'       ,'Inputsdlg'   , ...
            'HandleVisibility' ,'callback' , ...
            'Color'      ,FigColor  , ...
            'WindowStyle'   ,Options.WindowStyle, ...
            'DoubleBuffer'   ,'on'    , ...
            'Resize'      ,Options.Resize    ...
            );
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Create controls
        %%%%%%%%%%%%%%%%%%%%%
        CommonInfo = {'Units'  'pixels'
            'FontSize'      Options.FontSize
            'FontWeight'   get(0,'DefaultUicontrolFontWeight');
            'HandleVisibility'  'callback'}';
        
        props.edit = [CommonInfo {...
            'Style'      'edit';
            'HorizontalAlignment' 'left';
            'BackgroundColor' 'white'}'];
        
        props.checkbox = [CommonInfo {...
            'Style'      'checkbox';
            'HorizontalAlignment' 'center';
            'BackgroundColor' FigColor}'];
        
        props.popupmenu = [CommonInfo {...
            'Style'      'popupmenu';
            'HorizontalAlignment' 'left';
            'BackgroundColor' 'white'}'];
        
        props.listbox = [CommonInfo {...
            'Style'      'listbox';
            'HorizontalAlignment' 'left';
            'BackgroundColor' 'white'}'];
        
        props.slider = [CommonInfo {...
            'Style'      'slider';
            }'];
        
        props.uibuttongroup = [CommonInfo {...
            'BackgroundColor' FigColor
            }'];
        
        props.radiobutton = props.checkbox;
        props.radiobutton{2,strcmp(props.checkbox(1,:),'Style')} = 'radiobutton';
        
        props.pushbutton = [CommonInfo {...
            'Style'        'pushbutton';
            'HorizontalAlignment' 'center'}'];
        
        props.togglebutton = props.pushbutton;
        props.togglebutton{2,strcmp(props.pushbutton(1,:),'Style')} = 'togglebutton';
        
        props.table = CommonInfo;
        
        % Add VerticalAlignment here as it is not applicable to the above.
        props.text = [CommonInfo {...
            'Style'               'text'
            'Units'               'normalized'
            'Position'            [0 0 1 1]
            'Visible'             'off'}'];
        
        props.label = [CommonInfo {...
            'BackgroundColor'     FigColor;
            'HorizontalAlignment' 'left';
            'VerticalAlignment'   'bottom';
            'Color'        get(0,'FactoryUIControlForegroundColor');
            'Interpreter'     Options.Interpreter}'];
        
        % For each control, place a uicontrol and an axes (axes is to enable LaTeX)
        figsz = get(fig,'Position');
        
        props.uipanel = [CommonInfo {...
            'BackgroundColor' FigColor
            'BorderType','none'
            }.'];
        props.axes = [CommonInfo {...
            'Units' 'normalized'
            'Position' [0 0 1 1]
            'XLim' [0 figsz(3)]
            'YLim' [0 figsz(4)]
            'Visible' 'off'}.'];
        
        isbtngrp = false(num,1);
        istext = reshape(arrayfun(@(s)strcmp(s.type,'text'),Formats),num,1);
        autosize = zeros(num,2); % 0-fixed, 1-autosize, 2-resize with window
        rowspan = zeros(num,1); % # of rows control occupies
        colspan = zeros(num,1); % # of columns control occupies
        hPanels = zeros(num,2); % [uipanel|axes]
        hCtrls = zeros(num,3); % [Prompt|Edit|Unit]
        for m = 1:num % for each control
            
            % get current control's Format spec
            fmt = Formats(m);
            
            % determine the span of the control
            [i,j] = find(map==m);
            rowspan(m) = numel(unique(i));
            colspan(m) = numel(unique(j));
            
            % set autosize (if fixed, set only once)
            %    autosize(m,:) = (fmt.size<0) + (fmt.size<=0);
            autosize(m,:) = (fmt.size<=0);
            
            % if autosize, set to default size, if fixed
            if autosize(m,1) % width
                if fmt.size(1)==0
                    fmt.size(1) = 1; % temp size
                else
                    fmt.size(1) = -fmt.size(1);
                end
            end
            if autosize(m,2) % height
                if fmt.size(2)==0
                    fmt.size(2) = 1; % temp size
                else
                    fmt.size(2) = -fmt.size(2);
                end
            end
            
            % for uicontrols
            hPanels(m,1) = uipanel('Parent',fig,'Position',[0 0 fmt.size],props.uipanel{:});
            
            % for text controls & labels
            hPanels(m,2) = axes(props.axes{:},'Parent',hPanels(m,1)); %#ok
            
            % always add labels even if control is not using them
            hCtrls(m,2) = text('Parent',hPanels(m,2),props.label{:},'String',Prompt{m});
            hCtrls(m,3) = text('Parent',hPanels(m,2),props.label{:},'String',Units{m});
            
            idx = strcmp(fmt.style,{'radiobutton','togglebutton'});
            isbtngrp(m) = any(idx);
            if isbtngrp(m)
                % create the UI Button Group object
                h = uibuttongroup('Parent',hPanels(m,1),props.uibuttongroup{:},...
                    'Position',[0 0 fmt.size],'Tag',FieldNames{m});%,'Title',Prompt{m}
                hCtrls(m,1) = h;
                
                % Create button objects and record their extents
                dim_btns = size(fmt.items);
                kvalid = find(~cellfun(@isempty,fmt.items));
                Nvalid = numel(kvalid);
                hButtons = zeros(Nvalid,1);
                btn_w = zeros(dim_btns);
                btn_h = zeros(dim_btns);
                for n = 1:numel(kvalid)
                    [i,j] = ind2sub(dim_btns,kvalid(n));
                    hButtons(n) = uicontrol('Parent',h,'Style',fmt.style,props.(fmt.style){:},...
                        'String',fmt.items{i,j},'Min',0,'Max',1,'UserData',kvalid(n),...
                        'TooltipString',TooltipStrings{m});
                    pos = get(hButtons(n),'Extent');
                    btn_w(i,j) = pos(3);
                    btn_h(i,j) = pos(4);
                end
                
                % set buttons sizes and extra margins to account for the button width
                if idx(1) % radiobutton
                    % each column to have the same width
                    margin = [20 0];
                    btn_w = max(btn_w,[],1) + margin(1); % button widths with extra pixels for non-text part of control
                else % togglebutton
                    % all buttons to have the same width
                    margin = [12 2];
                    btn_w = repmat(max(btn_w(:)) + margin(1),1,dim_btns(2)); % button widths with extra pixels for non-text part of control
                end
                btn_h = max(btn_h,[],2) + margin(2); % button heights with extra pixels for non-text part of conrol
                
                % Button positions
                btn_sep = sep*3/4;
                x0 = cumsum([0 btn_w]+btn_sep)-btn_sep*3/8;
                y0 = flipud(cumsum([0;btn_h]+btn_sep)-btn_sep*3/8);
                
                % set positions of buttons
                kvalid = find(hButtons~=0);
                for n = 1:Nvalid
                    [i,j] = ind2sub(dim_btns,kvalid(n)); % i-col, j-row
                    pos = [x0(j) y0(i+1) btn_w(j) btn_h(i)];
                    set(hButtons(n),'Position',pos);
                end
                
                % set the size
                set(h,'Position',[0 0 x0(end) y0(1)],'UserData',hButtons);
                autosize(m,:) = 0; % no autosize
                
            elseif strcmp(fmt.style,'table') % uitable
                
                hCtrls(m,1) = uitable('Parent',hPanels(m,1), props.(fmt.style){:}, ...
                    'Position',[0 0 fmt.size],'ColumnName',fmt.items,...
                    'ColumnFormat',fmt.format,'ColumnWidth',fmt.limits,...
                    'ColumnEditable',true,...
                    'Tag',FieldNames{m});
                
            else % uicontrols
                
                % create a uipanel, embed a uicontrol, and 2 text labels
                hc = uicontrol('Parent',hPanels(m,1), 'Style',fmt.style, props.(fmt.style){:},...
                    'Position',[0 0 fmt.size],'Tag',FieldNames{m});
                hCtrls(m,1) = hc;
                
                % set min and max if not a numeric edit box
                if ~any(strcmp(fmt.style,'edit')) || ~any(strcmp(fmt.format,{'float','integer','date','dir'}))
                    lim = fmt.limits;
                    if any(isinf(lim)), lim = [0 2]; end % edit:vector
                    set(hc,'Min',lim(1),'Max',lim(2));
                end
                
                % style-dependent configuration
                switch fmt.style
                    case 'text' % static text (only a label)
                        
                        % display the label on the upper left hand corner of the display grid cell
                        %%%% DEFAULT :: set(hCtrls(m,2),'Units','normalized','Position',[0 1],'VerticalAlignment','top');
                        set(hCtrls(m,2),'Units','normalized','Position',[0 0.1],'VerticalAlignment','bottom');
                        
                        % create invisible dummy control
                        set(hc,'FontName',get(hCtrls(m,2),'FontName'),'String',Prompt{m});
                        
                        % if not autowidth, go ahead and map-out
                        if ~autosize(m,1)
                            set(hc,'Position',[0 0 fmt.size]);
                            msg = textwrap(h,Prompt(m));
                            str = sprintf('%s\n',msg{:});
                            set(hCtrl(m,1),'String',str(1:end-1));
                            autosize(m,2) = 0; % only auto-height if auto-width
                        end
                        
                        % Use Axes Text object (in order to render LaTeX text)
                        set(hCtrls(m,3),'String','','Visible','off');
                        
                    case 'edit' % edit type
                        
                        % check if multi-line control
                        dlim = round(diff(fmt.limits));
                        multiline = strcmp(fmt.format,'vector') || (any(strcmp(fmt.format,{'text','file'})) && dlim>1);
                        
                        % change alignment for numeric formats
                        if any(strcmp(fmt.format,{'float','integer'}))
                            set(hc,'HorizontalAlignment','center');
                        elseif strcmp(fmt.format,'vector')
                            set(hc,'HorizontalAlignment','left');
                        end
                        
                        % set auto-height (no resize allowed if single-line)
                        if autosize(m,2)>0 % auto-height adjustment
                            % set
                            if multiline % set to have dlim lines
                                if strcmp(fmt.format,'vector')
                                    nrows = max(fmt.limits(1),min(fmt.limits(2),5));
                                else
                                    nrows = dlim-1;
                                end
                                set(hc,'String',repmat(sprintf(' \n'),1,nrows-1));
                            else % single-line
                                set(hc,'String',' ');
                                autosize(m,2) = 0; % no need to adjust dynamically
                            end
                            ext = get(hc,'Extent');
                            set(hc,'String','');
                            fmt.size(2) = ext(4); % set to the font height
                            set(hc,'Position',[0 0 fmt.size]);
                        end
                        
                    case 'checkbox' % no labels
                        
                        % Show the prompt label with the control
                        set(hc,'String',Prompt{m});
                        set(hCtrls(m,2),'String','','Visible','off');
                        
                        % Set the control size (fixed, no resize)
                        pos = get(hc,'Extent'); % width reflects only label width
                        pos(3) = pos(3) + 20; % pad extra for the checkbox itself
                        set(hc,'Position',pos);
                        autosize(m,:) = 0; % no resizing
                        
                    case 'popupmenu'
                        
                        % get the width of the widest entry
                        if autosize(m,1) % auto-width
                            w = 0;
                            for n = 1:numel(fmt.items)
                                set(hc,'String',fmt.items{n});
                                ext = get(hc,'Extent');
                                w = max(w,ext(3));
                            end
                            fmt.size(1) = w + 20; % additional width for the pulldown
                        end
                        
                        if autosize(m,2) % auto-height
                            % only list 1 item and get the extent
                            set(hc,'String',fmt.items{1});
                            ext = get(hc,'Extent');
                            fmt.size(2) = ext(4);
                        end
                        
                        % re-set position
                        if any(autosize(m,:))
                            set(hc,'Position',[0 0 fmt.size]);
                            autosize(m,:) = 0; % no resizing
                        end
                        
                        % Set menu & choose the first entry
                        set(hc, 'String',fmt.items);
                        
                    case 'listbox'
                        
                        % Set menu & choose the first entry
                        set(hc,'String',fmt.items);
                        
                        if any(autosize(m,:))
                            % determine the optimal size
                            ext = get(hc,'Extent');
                            if autosize(m,1) % auto-width
                                fmt.size(1) = ext(3) + 20;
                                if autosize(m,2)~=1 && fmt.size(2)<ext(4) % with vertical scroller
                                    fmt.size(1) = fmt.size(1); % add vertical scroll bar width
                                end
                                autosize(m,1) = 0; % no resizing
                            end
                            if autosize(m,2) % auto-height -> set to the tallest
                                % Restrict the height if a maximum number of lines was specified
                                if fmt.limits(1) > 0 && fmt.limits(1) < numel(fmt.items)
                                    set(hc,'String',fmt.items(1:fmt.limits(1)));
                                    ext = get(hc,'Extent');
                                    set(hc,'String',fmt.items);
                                end
                                fmt.size(2) = ext(4);
                                if fmt.limits(1)>0
                                    autosize(m,2) = 0;
                                end
                            end
                            
                            % re-set position
                            set(hc,'Position',[0 0 fmt.size]);
                        end
                        
                    case 'slider'
                        
                        if any(autosize(m,:))
                            if autosize(m,1) && ~autosize(m,2) % vertical slider, auto-width -> fixed width
                                fmt.size(1) = 16;
                                autosize(m,1) = 0;
                            elseif autosize(m,2) % auto-height -> fixed height
                                fmt.size(2) = 16;
                                autosize(m,2) = 0;
                            end
                            set(hc,'Position',[0 0 fmt.size]);
                        end
                        
                    case 'pushbutton' % button & color types
                        
                        if strcmp(fmt.type,'color') % color type
                            set(hc,'String',' ') % space to set auto-height
                        else % button type
                            % Show the prompt label with the control
                            set(hc,'String',Prompt{m});
                            set(hCtrls(m,2),'String','','Visible','off');
                        end
                        
                        if autosize(m,2) % auto-height -> fix
                            ext = get(hc,'Extent');
                            fmt.size(2) = ext(4) + 6;
                            set(hc,'Position',[0 0 fmt.size]);
                            autosize(m,2) = 0;
                        end
                end
            end
        end % for m = 1:num
        
        % layout each control and labels within the panel (and axes)
        labelbase = zeros(num,2,2); % lower-left corner coordinates
        labelpos = zeros(num,2,2);
        labelsiz = zeros(num,2,2);
        labelalign = ones(num,4);
        ctrlpos = zeros(num,4);
        set(hCtrls(istext,2),'Units','pixels');
        alignprops = cell(2,2); % {'HorizontalAlignemt','VerticalAlignment'}
        for m = 1:num
            
            fmt = Formats(m);
            ext = cell2mat(get(hCtrls(m,[2 3]),'Extent'));
            if istext(m) % text type does not use uicontrol
                if autosize(m,1) % if autosize, minimum size to be the square
                    ext(1,3) = ext(1,4);
                end
                pos = [0 0 0 0];
            else
                pos = get(hCtrls(m,1),'Position');
                if isbtngrp(m)
                    pos(1) = pos(1) + 1;
                end
            end
            wd = [ext(1,3);pos(3);ext(2,3)];
            ht = [ext(1,4);pos(4);ext(2,4)];
            x0 = zeros(2,1); x = x0;
            y0 = zeros(2,1); y = y0;
            
            % place prompt label w.r.t. the control @ (0,0)
            tok = regexp(fmt.labelloc,'^(left|top)(.+)$','tokens','once');
            labelalign(m,1) = strcmp(tok{1},'left');
            if labelalign(m,1) % left
                % if control height is autoadjusted, make sure it is higher than label height
                if autosize(m,2)
                    ht(2) = max(ht(2),ht(1));
                end
                
                alignprops{1,1} = 'left'; alignprops(1,2) = tok(2);
                labelalign(m,2) = find(strcmp(tok{2},{'bottom','middle','top'}))-1;
                x0(1) = -wd(1) - fmt.margin(1); % move label to left
                x(1) = x0(1);
                y0(1) = (ht(2)-ht(1))*labelalign(m,2)/2;
                y(1) = ht(2)*labelalign(m,2)/2;
                
            else % above
                % if control width is autoadjusted, make sure it is wider than label width
                if autosize(m,1)
                    wd(2) = max(wd(2),wd(1));
                end
                
                alignprops{1,2} = 'bottom'; alignprops(1,1) = tok(2);
                labelalign(m,2) = find(strcmp(tok{2},{'left','center','right'}))-1;
                x0(1) = (wd(2)-wd(1))*labelalign(m,2)/2;
                x(1) = wd(2)*labelalign(m,2)/2;
                y0(1) = ht(2) + fmt.margin(1); % move label up
                y(1) = y0(1);
            end
            
            % units label w.r.t. the control @ (0,0)
            tok = regexp(fmt.unitsloc,'^(right|bottom)(.+)$','tokens','once');
            labelalign(m,3) = strcmp(tok{1},'right');
            if labelalign(m,3) % right
                % if control height is autoadjusted, make sure it is higher than unit height
                if autosize(m,2)
                    ht(2) = max(ht(2),ht(3));
                end
                
                alignprops{2,1} = 'left'; alignprops(2,2) = tok(2);
                labelalign(m,4) = find(strcmp(tok{2},{'bottom','middle','top'}))-1;
                x0(2) = wd(2) + fmt.margin(2); % move label right
                x(2) = x0(2);
                y0(2) = (ht(2)-ht(3))*labelalign(m,4)/2;
                y(2) = ht(2)*labelalign(m,4)/2;
            else % below
                % if control width is autoadjusted, make sure it is wider than unit width
                if autosize(m,1)
                    wd(2) = max(wd(2),wd(3));
                end
                
                alignprops{2,2} = 'top'; alignprops(2,1) = tok(2);
                labelalign(m,4) = find(strcmp(tok{2},{'left','center','right'}))-1;
                x0(2) = (wd(2)-wd(3))*labelalign(m,4)/2;
                x(2) = wd(2)*labelalign(m,4)/2;
                y0(2) = -ht(3)-fmt.margin(2); % move label down
                y(2) = -fmt.margin(2);
            end
            
            % translate so that all coordinates are on the first quadrant
            xmin = min(min(x0),0);
            ymin = min(min(y0),0);
            
            labelbase(m,1,:) = x0 - xmin;
            labelbase(m,2,:) = y0 - ymin;
            labelpos(m,1,:) = x - xmin;
            labelpos(m,2,:) = y - ymin;
            labelsiz(m,:,:) = ext(:,[3 4]).';
            
            ctrlpos(m,:) = [-xmin -ymin wd(2) ht(2)];
            if isbtngrp(m)
                ctrlpos(m,1) = ctrlpos(m,1) + 1;
                ctrlpos(m,3) = ctrlpos(m,3) - 1;
            end
            
            if ~istext(m)
                set(hCtrls(m,[2 3]),{'HorizontalAlignment','VerticalAlignment'},alignprops);
            end
        end
        set(hCtrls(istext,2),'Units','normalized');
        
        % minimum panel size
        pnpos = [zeros(num,2) max(max(labelbase+labelsiz,[],3),ctrlpos(:,[1 2])+ctrlpos(:,[3 4]))];
        
        % if strcmpi(Options.AutoLayout,'on')
        %    % minimum panel size
        %    pnpos = [zeros(num,2) max(max(labelbase+labelsiz,[],3),ctrlpos(:,[1 2])+ctrlpos(:,[3 4]))];
        %    map = autolayout(pnpos,Options.DesiredFigureWidth,sep);
        %    dim = size(map);
        % end
        
        % Optionally align adjacent controls in each column
        if strcmpi(Options.AlignControls, 'on')
            for n = 1:dim(2) % for each column
                if n==1
                    notspanned = true(dim(1),1);
                else
                    notspanned = diff(map(:,[n-1 n]),[],2)~=0;
                end
                idx = setdiff(unique(map(notspanned,n)),0); % must be a control & remove duplicates
                x0 = max(ctrlpos(idx,1));
                dx = x0-ctrlpos(idx,1);
                ctrlpos(idx,1) = ctrlpos(idx,1) + dx;
                labelpos(idx,1,:) = bsxfun(@plus,labelpos(idx,1,:),dx);
                labelbase(idx,1,:) = bsxfun(@plus,labelbase(idx,1,:),dx);
            end
            
            % recompute the minimum panel sizes
            pnpos = [zeros(num,2) max(max(labelbase+labelsiz,[],3),ctrlpos(:,[1 2])+ctrlpos(:,[3 4]))];
        end
        
        % position the controls w/in their respective panel, and set panel size
        nottext = ~istext;
        Nnottext = sum(nottext);
        set(hCtrls(nottext,1),{'Position'},mat2cell(ctrlpos(nottext,:),ones(Nnottext,1),4));
        set(hCtrls(nottext,2),{'Position'},mat2cell(labelpos(nottext,:,1),ones(Nnottext,1),2));
        set(hCtrls(:,3),{'Position'},mat2cell(labelpos(:,:,2),ones(num,1),2));
        set(hPanels(:,1),{'Position'},mat2cell(pnpos,ones(num,1),4));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine the maximum and minimum column widths & row heights
        
        idx = map>0;
        wctrl = (pnpos(:,3)-sep*(colspan-1))./colspan;
        wcell = zeros(dim); % -> widths of cells
        wcell(idx) = wctrl(map(idx));
        hctrl = (pnpos(:,4)-sep*(rowspan-1))./rowspan; % heights of controls
        hcell = zeros(dim);
        hcell(idx) = hctrl(map(idx)); % -> heights of cells
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine which columns & rows to autosize. Only column or row is
        % autosized for spanned controls, favoring the left-most column or
        % bottom-most row
        
        % get preferred column to autosize
        % Nauto = autosize(:,1);
        Nauto = false(dim);
        Nauto(idx) = autosize(map(idx),1);
        Nauto = sum(Nauto,1); % number of autowidth controls in each column
        [~,col] = sort(fliplr(Nauto),'descend');
        col(:) = dim(2)-col+1; % preferred order of columns to be autosized
        
        % assign which column to autosize for each control
        Ictrl = find(autosize(:,1)>0); % controls with auto-sized width
        autowidth = false(1,dim(2));
        for m = Ictrl.' % for each control
            [~,j] = find(map==m);
            [tf,I] = ismember(col,j);
            autowidth(j(I(find(tf,1)))) = true;
        end
        
        % get preferred row to autosize
        Nauto = false(dim);
        Nauto(idx) = autosize(map(idx),2);
        Nauto = sum(Nauto,2); % number of autowidth controls in each column
        [~,row] = sort(flipud(Nauto),'descend');
        row(:) = dim(1)-row+1; % preferred order of columns to be autosized
        
        % assign which row to autosize for each control
        Ictrl = find(autosize(:,2)>0&~istext);
        autoheight = false(dim(1),1);
        for m = Ictrl.'
            [i,~] = find(map==m);
            [tf,I] = ismember(row,i);
            autoheight(i(I(find(tf,1)))) = true;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Create the button panel
        
        % All buttons on a uipanel
        hBtnPanel = uipanel('Parent',fig,props.uipanel{:});
        
        % OK Button
        hBtn(1) = uicontrol(hBtnPanel, props.pushbutton{:}, ...
            'String', Options.ButtonNames{1}, 'UserData', 'OK');
        
        editable = ~all(istext);
        
        % Cancel Button
        if editable && strcmpi(Options.CancelButton,'on')
            hBtn(2) = uicontrol(hBtnPanel, props.pushbutton{:}, ...
                'String', Options.ButtonNames{2}, 'UserData', 'Cancel');
        end
        
        % Apply Button
        if editable && strcmpi(Options.ApplyButton,'on')
            hBtn(end+1) = uicontrol(hBtnPanel, props.pushbutton{:}, ...
                'String', Options.ButtonNames{3}, 'UserData', 'Apply');
        end
        
        % set size
        offset = 25;
        minwidth = 37;
        Nbtns = numel(hBtn);
        ext = cell2mat(get(hBtn,{'Extent'}));
        btnw = max(minwidth,max(ext(:,3)))+offset;
        btnh = max(ext(:,4)) + 6;
        btnsize = repmat([btnw btnh],Nbtns,1);
        btnpos = [(0:Nbtns-1).'*(btnw+sep)+sep repmat(sep,Nbtns,1)];
        
        set(hBtn,{'Position'},mat2cell([btnpos btnsize],ones(Nbtns,1),4));
        btnpanelsiz = [Nbtns*(btnw+sep)+sep btnh+2*sep];
        set(hBtnPanel,'Position',[0 0 btnpanelsiz]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % output handle struct
        handles.fig = fig;
        handles.panels = hPanels;
        handles.ctrls = hCtrls;
        handles.btnpanel = hBtnPanel;
        handles.btns = hBtn;
        
        % output positioning info
        sinfo.map = map; % Formats grid mapping
        sinfo.istext = istext; % true if static text
        sinfo.pnminsiz = pnpos(:,[3 4]); % minimum panel size
        sinfo.ctrlpos = ctrlpos; % control position w/in panel
        sinfo.labelpos = labelpos; % label positions w/in panel
        sinfo.labelalign = labelalign; % label alignments
        sinfo.btnpnsiz = btnpanelsiz; % button panel size
        
        sinfo.w_max = max(wcell,[],1);
        sinfo.w_min = min(wcell,[],1);
        sinfo.h_max = max(hcell,[],2);
        sinfo.h_min = min(hcell,[],2);
        
        sinfo.w_delta = sinfo.w_max-sinfo.w_min;
        sinfo.h_delta = sinfo.h_max-sinfo.h_min;
        
        % total non-adjustable width & height over all controls
        sinfo.w_totalnonadj = sum(sinfo.w_max(~autowidth))+sum(sinfo.w_min(autowidth));
        sinfo.h_totalnonadj = sum(sinfo.h_max(~autoheight))+sum(sinfo.h_min(autoheight));
        
        % autosize related fields
        sinfo.autosize = logical(autosize); % 1 to autosize once, 2 to resize with window
        sinfo.autowidth = find(autowidth);
        sinfo.autoheight = find(autoheight);
        
        figw = max(sum(sinfo.w_max)+sep*(dim(2)+1),sinfo.btnpnsiz(1)+2*sep);
        figh = sum(sinfo.h_max)+sep*(dim(1)+1)+sinfo.btnpnsiz(2);
        sinfo.figsizmin = min([figw figh], Options.DialogSizeMinimum);
        
        % Set default figure size
        Nauto = [numel(sinfo.autowidth) numel(sinfo.autoheight)];
        
        if strcmpi(Step,'initialize-histo') || strcmpi(Step,'initialize-immuno');
            scrsz = get(0,'ScreenSize');
            figpos = [0.21*scrsz(3) 0.3*scrsz(4) 0.1*scrsz(3) 0.1*scrsz(4)];
        elseif strcmpi(Step,'processing-options-histo')
            scrsz = get(0,'ScreenSize');
            figpos = [0.42*scrsz(3) 0.38*scrsz(4) 0.1*scrsz(3) 0.25*scrsz(4)];
        elseif strcmpi(Step,'processing-options-immuno')
            scrsz = get(0,'ScreenSize');
            figpos = [0.42*scrsz(3) 0.38*scrsz(4) 0.1*scrsz(3) 0.25*scrsz(4)];
        elseif strcmpi(Step,'threshold-manual')
            scrsz = get(0,'ScreenSize');
            figpos = [0.055*scrsz(3) 0.4*scrsz(4) 0.1*scrsz(3) 0.1*scrsz(4)];
        elseif strcmpi(Step,'selection')
            scrsz = get(0,'ScreenSize');
            figpos = [0.4*scrsz(3) 0.4*scrsz(4) 0.1*scrsz(3) 0.1*scrsz(4)];
        elseif isempty(Step);
            figpos = get(fig,'Position');
        end
        
        if Nauto(1)==0 % all fixed columns
            figpos(3) = sinfo.figsizmin(1);
        else
            %%%%% DEFAULT :: figpos(3) = max(figpos(3),sinfo.figsizmin(1)+59*Nauto(1));
            figpos(3) = max(figpos(3),sinfo.figsizmin(1)+0*Nauto(1));
        end
        if Nauto(2)==0 % all fixed rows
            figpos(4) = sinfo.figsizmin(2);
        else
            figpos(4) = max(figpos(4),sinfo.figsizmin(2)+19*Nauto(2));
        end
        set(fig,'Position',figpos);
        movegui(fig,'onscreen');
        
    end

    function resizegui(handles,sinfo,Options)
        % This function places all controls in proper place.
        % Must be called before the GUI is made visible as buildgui function
        % just creates uicontrols and do not place them in proper places.
        
        % KNOWN ISSUE: Resize event does not fire all the time
        
        sep = Options.Sep; % separation between controls
        
        % get current figure size
        figPos = get(handles.fig,'Position');
        figSize = figPos(3:4);
        
        % force figure size to be larger than the minimum size allowed
        idx = figSize<sinfo.figsizmin;
        if any(idx)
            figSize(idx) = sinfo.figsizmin(idx);
            figPos([3 4]) = figSize;
            set(handles.fig,'Position',figPos);
            movegui(handles.fig,'onscreen');
        end
        
        % Place the button panel (lower right hand corner)
        pos = [figSize(1)-sinfo.btnpnsiz(1)-sep sep sinfo.btnpnsiz];
        set(handles.btnpanel,'Position',pos);
        
        dim = size(sinfo.map);
        num = size(handles.panels,1);
        
        % Determine the column widths
        w_total = figSize(1)-sep*(dim(2)+1); % sans spacers
        autowidth = sinfo.autowidth;
        totalnonadj = sinfo.w_totalnonadj;
        w_adj = (w_total-totalnonadj)/numel(autowidth);
        widths = sinfo.w_max;
        widths(autowidth) = sinfo.w_min(autowidth) + w_adj;
        
        % make sure all the fixed width components fit in the auto-width columns
        idx = widths(autowidth)<sinfo.w_max(autowidth);
        while any(idx)
            % fix the row with maximum difference b/w max and min
            idx = autowidth(idx); % idx now contains actual row index
            [~,I] = max(sinfo.w_delta(idx)); % max diff
            idx = idx(I); % narrow down to 1
            autowidth = setdiff(autowidth,idx); % remove it from autowidth list
            widths(idx) = sinfo.w_max(idx); % fix the width to its max
            totalnonadj = totalnonadj + sinfo.w_max(idx) - sinfo.w_min(idx);
            
            % recompute the adjustable heights
            w_adj = (w_total-totalnonadj)/numel(autowidth);
            widths(autowidth) = sinfo.w_min(autowidth) + w_adj;
            idx = widths(autowidth)<sinfo.w_max(autowidth);
        end
        
        % Determine the grid's x-coordinates
        grid_x0 = cumsum([0 widths(1:end-1)]+sep);
        grid_x1 = grid_x0 + widths;
        
        % Set heights of autoadjusted text controls
        h_min = sinfo.h_min;
        heights = sinfo.h_max;
        idx = sinfo.istext & sinfo.autosize(:,1);
        for m = find(idx).'
            % get the grid cell position
            [i,j] = find(sinfo.map==m);
            x = min(grid_x0(j));
            pwidth = max(grid_x1(j))-x;
            ppos = [x 0 pwidth 1]; % panel position, only set width
            
            set(handles.panels(m,1),'Position',ppos);
            
            % use textwrap to obtain an initial wrapped candidate
            str = get(handles.ctrls(m,1),'String');
            msg = textwrap(handles.ctrls(m,1),{str});
            
            % often too wide, so reduce words on every line.
            % Known Bug: if there is an explicit newline in the string, there
            % is a possible for it to be ignored as a result of the code below
            Nmsgs = numel(msg);
            str = '';
            for k = 1:Nmsgs
                if k==1
                    str = msg{1};
                else
                    str = sprintf('%s\n%s',str,msg{k});
                end
                set(handles.ctrls(m,2),'String',str);
                ext = get(handles.ctrls(m,2),'Extent');
                while ext(3)>1
                    % send last word off to the next line
                    [idx,word] = regexp(str,'\s(\S+)$','start','tokens','once');
                    str(idx:end) = [];
                    set(handles.ctrls(m,2),'String',str);
                    ext = get(handles.ctrls(m,2),'Extent');
                    if k==Nmsgs && numel(msg)==Nmsgs
                        msg(k+1) = word;
                    else
                        msg{k+1} = sprintf('%s %s',word{1},msg{k+1});
                    end
                end
            end
            
            % set height
            i = unique(i);
            heights(i) = max(heights(i),ext(4)/numel(i));
            h_min(i) = min(h_min(i),ext(4)/numel(i));
        end
        
        % Determine the automatically adjusted row heights
        h_total = figSize(2)-sep*(dim(1)+1)-sinfo.btnpnsiz(2);
        autoheight = sinfo.autoheight;
        totalnonadj = sinfo.h_totalnonadj;
        h_adj = (h_total-totalnonadj)/numel(autoheight);
        heights(autoheight) = h_min(autoheight) + h_adj;
        
        idx = heights(autoheight)<sinfo.h_max(autoheight);
        while any(idx)
            % fix the row with maximum difference b/w max and min
            idx = autoheight(idx); % idx now contains actual row index
            [~,I] = max(sinfo.h_delta(idx));
            idx = idx(I);
            autoheight = setdiff(autoheight,idx);
            heights(idx) = sinfo.h_max(idx);
            totalnonadj = totalnonadj + sinfo.h_max(idx) - sinfo.h_min(idx);
            
            % recompute adjustable heights
            h_adj = (h_total-totalnonadj)/numel(autoheight);
            heights(autoheight) = sinfo.h_min(autoheight) + h_adj;
            idx = heights(autoheight)<sinfo.h_max(autoheight);
        end
        
        % grid y-coordinates
        grid_y0 = flipud(cumsum([0;flipud(heights(2:end))]+sep)) + sinfo.btnpnsiz(2);
        yoffset = figSize(2) - sep - grid_y0(1) - heights(1); % aligned to top edge
        grid_y0(:) = grid_y0 + yoffset;
        grid_y1 = grid_y0 + heights;
        
        % obtain position of each control panel
        for m = 1:num
            % get the grid cell position
            [i,j] = find(sinfo.map==m);
            x = min(grid_x0(j));
            y = min(grid_y0(i));
            ytop = max(grid_y1(i));
            pwidth = max(grid_x1(j))-x;
            ppos = [x y pwidth ytop-y]; % panel position
            
            % if fixed height, align the panel to the upper edge of the grid cell
            if ~sinfo.autosize(m,2)
                ppos(2) = ytop-sinfo.pnminsiz(m,2);
                ppos(4) = sinfo.pnminsiz(m,2);
            end
            
            set(handles.panels(m,1),'Position',ppos);
            
            if ~sinfo.istext(m) && any(sinfo.autosize(m,:)) % for all other control types
                
                cpos = {sinfo.ctrlpos(m,:);sinfo.labelpos(m,:,1);sinfo.labelpos(m,:,2)};
                
                % adjust the width of the control
                if sinfo.autosize(m,1) % auto width
                    
                    % determine the control width
                    dw = ppos(3)-sinfo.pnminsiz(m,1);
                    wd = cpos{1}(3) + dw; % new control width
                    cpos{1}(3) = wd;
                    
                    x0 = sinfo.ctrlpos(m,1);
                    if ~sinfo.labelalign(m,1) % prompt label at top
                        switch sinfo.labelalign(m,2)
                            case 2 % right
                                cpos{2}(1) = x0 + wd;
                            case 1 % center
                                cpos{2}(1) = x0 + wd/2;
                        end
                    end
                    
                    if sinfo.labelalign(m,3) % unit label to right
                        cpos{3}(1) = x0 + wd;
                    else % unit label at bottom
                        switch sinfo.labelalign(m,4)
                            case 2 % right
                                cpos{3}(1) = x0 + wd;
                            case 1 % center
                                cpos{3}(1) = x0 + wd/2;
                        end
                    end
                end
                
                if sinfo.autosize(m,2) % auto height
                    dh = ppos(4)-sinfo.pnminsiz(m,2);
                    ht = cpos{1}(4) + dh;
                    cpos{1}(4) = ht;
                    
                    y0 = sinfo.ctrlpos(m,2);
                    if sinfo.labelalign(m,1) % prompt label to left
                        switch sinfo.labelalign(m,2)
                            case 2 % top
                                cpos{2}(2) = y0 + ht;
                            case 1 % middle
                                cpos{2}(2) = y0 + ht/2;
                        end
                    else % prompt label at top
                        cpos{2}(2) = y0 + ht;
                    end
                    
                    if sinfo.labelalign(m,3) % units label to right
                        switch sinfo.labelalign(m,4)
                            case 2 % top
                                cpos{3}(2) = ht;
                            case 1 % middle
                                cpos{3}(2) = ht/2;
                        end
                    end
                end
                
                set(handles.ctrls(m,:),{'Position'},cpos);
            end
        end
    end


% Copyright (c) 2009-2013, Takeshi Ikuma
% Copyright (c) 2010, Luke Reisner
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%   * Neither the names of its contributors may be used to endorse or
%     promote products derived from this software without specific prior
%     written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

end