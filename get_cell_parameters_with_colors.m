close all;
clear all; 
clc;

global I analysis_foldername j

%% parameters
origin_filename = 'Mask2.tif'; % path + filename of the image/mask to be analyzed
cell_color = 'Black'; 
pixel_scale = 0.65; % conversion µm/pixel
micrometers = 'off'; % "on" = working in microns; "off" or else = working in pixels
dilate_strel = strel('diamond',1); % structuring element for dilation/erosion of binary images 
minimum_contact_length = @(perimeter) 0.05*perimeter; % minimum shared area ¿in pixels? for cells to be in contact
max_cycle_elements = 8; % maximum cycle size to look for in graphs
min_cycle_elements = 8; % minimum cycle size to look for in graphs
j=1; % number corresponding to the image (for file naming purposes)

%% Create analysis folder where data and images are saved
analysis_foldername = extractBefore(origin_filename,'.');
[status, msg, msgID] = mkdir(analysis_foldername);

%% flags --> 0 = off; 1 = on; 2 = on + plot; 3 (or >) = on + plot + save
image2binary_flag = 1; % load image, binarize it and invert it if necessary
borderoff_flag = 1; % erase cells at the borders 
smalloff_flag = 1; % erase small components
cc_flag = 1; % find all connected components = cells
celldata_flag = 2; % extract geometrical data from cells and visualize orientation
graph_flag = 1; % find adjacency matrix and graph
topocycles_flag = 1; % find cycles that have a certain topological charge
complexpoloff_flag = 1; % remove complex polygons
solidityfilter_flag = 0; % solidity filter
roundnessfilter_flag = 0; % roundness filter
naive_flag = 0; % naive on = convexhull
plusonedefs_flag = 2; % Look for +1 defects
minusonedefs_flag = 2; % Look for -1 defects
plushalfdefs_flag = 2; % Look for +1/2 defects
minushalfdefs_flag = 2; % Look for -1/2 defects

% visualize ellipses superimposed to original image
% visualize orientation vectors superimposed to original image

%% load image, binarize it and invert it if necessary
[binary_fig] = image2binary(image2binary_flag, origin_filename, cell_color);

%% erase cells at the borders 
[borderoff_fig] = borderoff(borderoff_flag);

%% erase small components and visualize binary image before and after filtering by size
[smalloff_fig] = smalloff(smalloff_flag);

%% find label and visualize all connected components = cells
[cc, labeled_cells, labeledcells_fig] = connectedcomp(cc_flag);

%% extract geometrical data from cells and applying the scaling if required
[cell_data, cell_orientation_fig] = celldata(celldata_flag, micrometers, cc, labeled_cells);

%% Find adjacency matrix and graph
[connectivity_graph, cells_xc, cells_yc, connectivity_fig] = conngraph(graph_flag, cc, labeled_cells, dilate_strel, cell_data, minimum_contact_length);

%% Find cycles that have a certain topological char ge
[contact_cycles, edgecycles, topologicalCharges, tiledfig, singlefig, ...
    plus_one_defs, plusone_x, plusone_y, ...
    minus_one_defs, minusone_x, minusone_y, ...
    plus_half_defs, plushalf_x, plushalf_y, ...
    minus_half_defs, minushalf_x, minushalf_y] = topocycles( ...
    topocycles_flag,topocycles_flag, solidityfilter_flag, roundnessfilter_flag, ...
    plusonedefs_flag, minusonedefs_flag, plushalfdefs_flag, minushalfdefs_flag, ...
    cell_data, connectivity_graph, max_cycle_elements, min_cycle_elements, cells_xc, cells_yc);


%% ellipse visualization thanks to its parametric equation superimposed to original image
figure;
imshow(K);

t = linspace(0,2*pi,50);
hold on
for k = 1:length(cell_data)
    a = cell_data (k).MajorAxisLength/2;
    b = cell_data (k).MinorAxisLength/2;
    Xc = cell_data (k).Centroid(1);
    Yc = cell_data (k).Centroid(2);
    phi = deg2rad(-cell_data(k).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'b','Linewidth',2);
end
hold off

%% Orientation vector visualization superimposed to original image
figure;
imshow(K);

hold on
for k = 1:length(cell_data)
    vlength = cell_major_axis(k);
    t = linspace(-vlength/2,vlength/2,3);
    Xc = cell_data(k).Centroid(1);
    Yc = cell_data(k).Centroid(2);
    phi = deg2rad(-cell_data(k).Orientation);
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
        for i = 1:length(cell_data)
          cellindex = i*ismember([pixelindex],[cc.PixelIdxList{i}]);
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
            for i = 1:length(cell_data)
                cellindexmat = cellindexmat + i*ismember([pixelindexmat],[cc.PixelIdxList{i}]);
            end
            if any(any(cellindexmat))
                cellindexmat(cellindexmat == 0) = NaN;
                cellindex = mode(cellindexmat, "all");
            end           
        end
        if cellindex < 1
            plot(Xc,Yc,'.','Color','g');
        else 
            phi = deg2rad(-cell_data(cellindex).Orientation);
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

clear *_flag
% warning('off')
% warning('on')
