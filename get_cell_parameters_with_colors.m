close all;
clear all; 
clc;
warning('off');

scale = 0.65; % 1/3.0769; % conversion µm/pixel

%% load image, binarize it and invert it

j=1; % number corresponding to the image
I=imread('Mask1.tif');
% figure, imshow(I);

I=imbinarize(I);
I=imcomplement(I);
% figure, imshow(I);

%% find and vizualise all connected components = cells
CC = bwconncomp(I);
labeled = labelmatrix(CC);
RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
% figure;
% imshow(RGB_label);
% title('connected components');
 
%% erase cells at the borders 
dilatedLabel = imclearborder(I);
% figure, imshow(dilatedLabel);

CC = bwconncomp(dilatedLabel);
labeled = labelmatrix(CC);
RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
% figure, imshow(RGB_label);
% title('connected components');


%% erase small components 
data = regionprops(CC,'Area');

Max_size = max([data.Area]);
Min_size = 0; %(mean([data.Area]))/10; 

K = bwareafilt(dilatedLabel,[Min_size Max_size]);
CC = bwconncomp(K);

%% visualize binary image before and after filtering by size
% figure;
% subplot(1,2,1);
% imshow(dilatedLabel);
% title('original binary image');
% subplot(1,2,2);
% imshow(K);
% title('filtered binary image');

%% extract area and perimeter of cells

data = regionprops(CC,{...
    'Area',...
    'Perimeter',...
    'Centroid',...
    'MajorAxisLength',...
    'MinorAxisLength',...
    'Orientation'});


table = struct2table(data); % convert structure into table

% converting data
cell_size = [data.Area].*(scale^2);  % conversion en µm²
cell_perimeter = [data.Perimeter].*scale; % conversion en µm
major_axis_length = [data.MajorAxisLength].*scale; % conversion en µm
minor_axis_length = [data.MinorAxisLength].*scale; % conversion en µm


% data are transposed to be put in the columns of table
table.Area = cell_size.';
table.Perimeter = cell_perimeter.';
table.Orientation = [data.Orientation].';
table.MajorAxisLength = major_axis_length.';
table.MinorAxisLength = minor_axis_length.';

%% visualize remaining cells

mal = abs(table.Orientation); %Plug in Orientation here.
[~,idx] = histc(mal,0:90);  
L = zeros(CC.ImageSize); %preallocate
    for ii = 1:CC.NumObjects
      L(CC.PixelIdxList{ii}) = idx(ii);    %fill in indices
    end
% cmap = parula(90);  %a colormap
% Lrgb = label2rgb(L,cmap); %build rgb image
% figure;
% imshow(Lrgb);
% % imwrite(Lrgb, [ '00' num2str(j) '-cell-colors.bmp']);
% caxis([0, 90])
% colorbar;
% title('Cell orientation (°)');

%% Find edges and adjacency matrix

se = strel('disk',4);
dilateLabel = imdilate(labeled,se);
LBL = imLabelEdges(dilateLabel);
% Show them
cmap2 = parula(max(max(LBL)));
LBLrgb = label2rgb(LBL,cmap2); %build rgb image
% figure;
% imshow(LBLrgb);
% title('Edges');

adjacencyMatrix = zeros(max(max(dilateLabel)));
N = double(max(dilateLabel(:)));
for c=1:N
    lbl = immultiply(dilateLabel, imdilate(dilateLabel==c, se));
    lbl = lbl(lbl~=0);
    % Removing contacts between cells that share very little space
    [GC,lbl] = groupcounts(lbl);
    lbl = lbl(GC>20);
    %
    adjacencyMatrix(c,lbl) = 1;
end
adjacencyMatrix = adjacencyMatrix | adjacencyMatrix';

g = graph(adjacencyMatrix,'omitselfloops');

% remove
xydata = vertcat(data.Centroid);
xdata = xydata(:,1);
ydata = xydata(:,2);
figure
imshow(K);
hold on
plot(g,'XData',xdata,'YData',ydata) 
% remove

maxcycle = 8;
mincycle = 8;
[cycles,edgecycles] = allcycles(g,'MaxCycleLength',maxcycle,'MinCycleLength',mincycle);

xydata_cell = {data.Centroid}';
xycycle_cell = xydata_cell(cell2mat(cycles)); 
xycycle_cell = mat2cell(cell2mat(xycycle_cell),ones(size(xycycle_cell,1),1),2*maxcycle);

% remove complex polygons
complex_pol_index = @(cycle_xy) polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).NumRegions;
complex_pol_indices = cellfun(complex_pol_index,xycycle_cell);
cycles = cycles(complex_pol_indices == 1);
edgecycles = edgecycles(complex_pol_indices == 1);
xycycle_cell = xycycle_cell(complex_pol_indices == 1);

% sortedcycles
clockwise_index = @(cycle_xy) ispolycw(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)]));
clockwise_indices = cellfun(clockwise_index,xycycle_cell);
ccw_cycle = @(cycle, cycle_xy) [flip(cycle)*clockwise_index(cycle_xy) + cycle*(1-clockwise_index(cycle_xy))];
ccw_cycles = cellfun(ccw_cycle,cycles,xycycle_cell,'UniformOutput',false);
ccw_edgecycles = cellfun(ccw_cycle,edgecycles,xycycle_cell,'UniformOutput',false);

% solidity filter
solidity = @(cycle_xy)...
    polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).area /...
    polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).convhull.area;
solidities = cellfun(solidity,xycycle_cell);
% solid_ccw_cycles = ccw_cycles(solidities > 0.9);
% solid_ccw_edgecycles = ccw_edgecycles(solidities > 0.9);

% roundness filter % AND SOLIDITY
roundness = @(cycle_xy)...
    4*pi*polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).area /...
    (polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).perimeter)^2;
roundnesses = cellfun(roundness,xycycle_cell);
round_solid_ccw_cycles = ccw_cycles(logical((solidities > 0.9).*(roundnesses > 0.8)));
round_solid_ccw_edgecycles = ccw_edgecycles(logical((solidities > 0.9).*(roundnesses > 0.8)));

topo_wrapper = @(cell_cycle) topological_charge(cell_cycle, data);

topologicalCharges = cellfun(topo_wrapper,round_solid_ccw_cycles);
plusOneDefs = find(topologicalCharges==1);

cmap3 = parula(length(plusOneDefs));
cmap3 = cmap3(randperm(size(cmap3, 1)), :);

figure
% tiledlayout flow
imshow(K);
hold on

for i = 1:length(plusOneDefs)
%     nexttile
%     imshow(K)
%     hold on
    highlight(plot(g,'XData',xdata,'YData',ydata,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',round_solid_ccw_edgecycles{plusOneDefs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
%     title("Defect " + i)
%     labels = {'1','2','3','4','5','6','7','8'};
%     plot(xdata(solid_ccw_cycles{plusOneDefs(i)}),ydata(solid_ccw_cycles{plusOneDefs(i)}),'Marker','.','MarkerSize',6,'MarkerFaceColor','r')
%     text(xdata(solid_ccw_cycles{plusOneDefs(i)}),ydata(solid_ccw_cycles{plusOneDefs(i)}),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
%     xmean = mean(xdata(solid_ccw_cycles{plusOneDefs(i)}));
%     ymean = mean(ydata(solid_ccw_cycles{plusOneDefs(i)}));
%     plot(xmean,ymean,'Marker','*','MarkerSize',6,'MarkerFaceColor','r')
end
plot(g,'XData',xdata,'YData',ydata,'EdgeColor','#0072BD','NodeColor','#0072BD')

