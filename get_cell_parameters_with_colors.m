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

%% ellipse visualization thanks to its parametric equation superimposed to original image
figure;
imshow(K);

t = linspace(0,2*pi,50);
hold on
for k = 1:length(data)
    a = data (k).MajorAxisLength/2;
    b = data (k).MinorAxisLength/2;
    Xc = data (k).Centroid(1);
    Yc = data (k).Centroid(2);
    phi = deg2rad(-data(k).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'b','Linewidth',2);
end
hold off

%% Orientation vector visualization superimposed to original image
figure;
imshow(K);

hold on
for k = 1:length(data)
    vlength = major_axis_length(k);
    t = linspace(-vlength/2,vlength/2,3);
    Xc = data(k).Centroid(1);
    Yc = data(k).Centroid(2);
    phi = deg2rad(-data(k).Orientation);
    x = Xc + t*cos(phi);
    y = Yc + t*sin(phi);
    plot(x,y,'b','Linewidth',2);
end
hold off
grid on
xticks(0:10:500)
yticks(0:10:500)

% Structured orientation vector visualization superimposed to original image
figure;
imshow(K);
vlength = 14;
xstep = 18;
ystep = 18;
xcoords = 0:xstep:475;
xcoords = xcoords + vlength/2;
ycoords = 0:ystep:470; %475
ycoords = ycoords + vlength/2;
[X,Y] = meshgrid(xcoords,ycoords);
phimat = zeros(length(ycoords),length(xcoords));

hold on
for k = 1:length(xcoords)
    for h = 1:length(ycoords)
        t = linspace(-vlength/2,vlength/2,3);
        Xc = X(h,k);
        Yc = Y(h,k);
        pixelindex = ((k-1)*ystep+vlength/2)*470 + (h-1)*xstep + vlength/2;
        for i = 1:length(data)
          cellindex = i*ismember([pixelindex],[CC.PixelIdxList{i}]);
              if cellindex > 0
                  break
              end
        end
        if cellindex < 1
            checksize = 5;
            pixelindexmat = zeros(1+2*checksize);
            cellindexmat = zeros(1+2*checksize);
            pixelindexmat(1+checksize,:) = pixelindex;
            for row = 1:checksize
                pixelindexmat(row,:) = pixelindex - 470*(checksize-row+1);
                pixelindexmat(1+2*checksize-row+1,:) = pixelindex + 470*(checksize-row+1);
            end
            for col = 1:checksize
                pixelindexmat(:,col) = pixelindexmat(:,col) - (checksize-col+1);
                pixelindexmat(:,1+2*checksize-col+1) = pixelindexmat(:,1+2*checksize-col+1) + (checksize-col+1);
            end
            for i = 1:length(data)
                cellindexmat = cellindexmat + i*ismember([pixelindexmat],[CC.PixelIdxList{i}]);
            end
            if any(any(cellindexmat))
                cellindexmat(cellindexmat == 0) = NaN;
                cellindex = mode(cellindexmat, "all");
            end           
        end
        if cellindex < 1
            plot(Xc,Yc,'.','Color','g');
        else 
            phi = deg2rad(-data(cellindex).Orientation);
            phimat(h,k) = phi;
            phimat(phimat == 0) = NaN;
            x = Xc + t*cos(phi);
            y = Yc + t*sin(phi);
            plot(x,y,'r','Linewidth',2);
        end
    end
end
hold off
grid on
xticks(0:20:500)
yticks(0:20:500)

%% Saving pictures
name_picture=sprintf('00%d-Colors.tif',j);
saveas(gcf,name_picture);

%% calculations Shape Index
table.SI = 4 *pi* table.Area ./ (table.Perimeter.^2);
table.SIjamming = table.Perimeter./sqrt(table.Area);

table_name=sprintf('00%d-CY5.csv',j);
% writetable(table,table_name);

plot(g,'XData',xdata,'YData',ydata,'EdgeColor','#0072BD','NodeColor','#0072BD')