function histRGB(I,Lthresh,Uthresh,Lcolor,Ucolor)

% Isolate R-,G-, and B- channels from image
R = I(:,:,1);    G = I(:,:,2);    B = I(:,:,3);

% Define image size parameters
figure(1)
h2 = subplot('Position',[0.05 0.05 0.93 0.2]);

% Clear current axis
cla(h2)

% Create histograms of all non-white pixel values in RED channel
[countsr,x] = imhist(R(R~=255));	
area(x,countsr,'FaceColor','r','EdgeColor','r','LineWidth',1.5);
hold on

% Create histograms of all non-white pixel values in GREEN channel
[countsg,x] = imhist(G(G~=255));	
area(x,countsg,'FaceColor','g','EdgeColor','g','LineWidth',1.5);

% Create histograms of all non-white pixel values in BLUE channel
[countsb,x] = imhist(B(B~=255));
area(x,countsb,'FaceColor','b','EdgeColor','b','LineWidth',1.5);

% Histogram transparency
alpha(0.6)

% Set plot boundaries
xmin = 0;   xmax = 255;                                 xlim([xmin xmax])
ymin = 0;   ymax = max([countsr; countsg; countsb]);    ylim([ymin 1.02*ymax])

% Plot vertical lines for current threshold values
plot([Lthresh Lthresh],[ymin 1.02*ymax],'Color',Lcolor,'LineStyle','-','LineWidth',1.25)
plot([Uthresh,Uthresh],[ymin 1.02*ymax],'Color',Ucolor,'LineStyle','-','LineWidth',1.25)

% Plot threshold region as rectangle
area([Lthresh Uthresh],[1.02*ymax 1.02*ymax],'FaceColor','k','EdgeColor','none');
alpha(0.4)

