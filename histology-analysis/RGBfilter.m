function I = RGBfilter(A,Lthresh,Uthresh)

global stain SF

% -------------------------------------------------------------------------
% Threshold RGB image based on background color (white or black)
% -------------------------------------------------------------------------

% Mask pixels outside of threshold bounds and assign background color
if strcmpi(stain,'dPSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
    Im = (Lthresh <= A & A <= Uthresh) | (10 >= A);     Im = all(Im,3);     bground = 0;    
else
    Im = (Lthresh <= A & A <= Uthresh) | (245 <= A);    Im = all(Im,3);     bground = 255;    
end

% Remove background in each image channel
A1 = A(:,:,1);      A1(Im == 1) = bground;
A2 = A(:,:,2);      A2(Im == 1) = bground;
A3 = A(:,:,3);      A3(Im == 1) = bground;

% Recompile thresholded image
A = cat(3,A1,A2,A3);


% Generate binary mask of thresholded image
if strcmpi(stain,'dPSR') || strcmpi(stain,'POL') || strcmpi(stain,'EnF') || strcmpi(stain,'IF')
    Abw = im2bw(A,0.01);
else
    Abw = imcomplement(im2bw(A,0.99));
end

% Extract area of each pixel group in binary mask
bwprop = regionprops(Abw,'Area','PixelList');

% Initialize and write pixel areas from structure to array
bwarea = zeros(length(bwprop(:,1)),1);

for II = 1:length(bwprop(:,1))
    bwarea(II,1) = bwprop(II,1).Area;
end

% Remove 'small' pixel groups
for II = 1:length(bwarea);
    
    % Define area threshold...
    if bwarea(II,1) < round(15*SF(1)*SF(2));
        
        % Extract number of pixels in the current group...
        num = length(bwprop(II,1).PixelList(:,1));
        
        % Recolor to be the same as the background...
        for JJ = 1:num
            A(bwprop(II,1).PixelList(JJ,2),bwprop(II,1).PixelList(JJ,1),:) = bground;
        end
        
    end
    
end


% Redefine thresholded image
I = A;

% Show result in window
figure(1)
subplot('Position',[0.05 0.29 0.92 0.68]);
imshow(I)
