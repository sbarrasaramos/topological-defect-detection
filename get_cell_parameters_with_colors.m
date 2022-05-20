close all;
clear all; 
clc;

%% parameters
origin_filename = 'Mask2.tif'; % path + filename of the image/mask to be analyzed
cell_color = 'Black'; 
pixel_scale = 0.65; % conversion µm/pixel
dilate_strel = strel('diamond',1); % structuring element for dilation/erosion of binary images 
minimum_contact_length = @(perimeter) 0.05*perimeter; % minimum shared area ¿in pixels? for cells to be in contact
max_cycle_elements = 8; % maximum cycle size to look for in graphs
min_cycle_elements = 8; % minimum cycle size to look for in graphs
j=1; % number corresponding to the image (for file naming purposes)

%% Create analysis folder where data and images are saved
analysis_foldername = extractBefore(origin_filename,'.');
[status, msg, msgID] = mkdir(analysis_foldername);

%% flags --> 0 = off; 1 = on; 2 = on + plot; 3 (or >) = on + plot + save
image2binary_flag = 3; % load image, binarize it and invert it if necessary
borderoff_flag = 3; % erase cells at the borders 
smalloff_flag = 3; % erase small components
cc_flag = 3; % find all connected components = cells
celldata_flag = 3; % extract geometrical data from cells and visualize orientation
orientation_plottype = 'Colormap'; % choose how to show cell orientation: 'colormap', 'ellipse', 'major_axis'
graph_flag = 3; % find adjacency matrix and graph
topocycles_flag = 3; % find cycles that have a certain topological charge
complexpoloff_flag = 1; % remove complex polygons
solidityfilter_flag = 1; % solidity filter
roundnessfilter_flag = 1; % roundness filter
naive_flag = 0; % naive on = convexhull
plusonedefs_flag = 3; % Look for +1 defects
minusonedefs_flag = 3; % Look for -1 defects
plushalfdefs_flag = 3; % Look for +1/2 defects
minushalfdefs_flag = 3; % Look for -1/2 defects

%% load image, binarize it and invert it if necessary
if image2binary_flag > 0
    I=imbinarize(imread(origin_filename));
    switch lower(cell_color)
        case 'black'
            I=imcomplement(I);
        case 'white'
        otherwise
            disp('Image is not a mask')
    end
    if image2binary_flag > 1
        figure, imshow(I);
        if image2binary_flag > 2
            saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-binary.tif',j)));
        end
    end
else
    error('You cannot perform an analysis without a binary image. Set image2binary_flag to 1 or a higher value.')
end

%% erase cells at the borders 
if borderoff_flag > 0
    I = imclearborder(I);
    if borderoff_flag > 1
        figure, imshow(I);
        if borderoff_flag > 2
            saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-binary_WOEdges.tif',j)));
        end
    end
else
    warning("Parameters obtained from cells at the borders might result in biased or incorrect results.")
end

%% erase small components and visualize binary image before and after filtering by size
if smalloff_flag > 0
    cell_sizes = regionprops(I, 'Area');
    max_cell_size = max([cell_sizes.Area]);
    min_cell_size = mean([cell_sizes.Area])/10; 
    I = bwareafilt(I,[min_cell_size max_cell_size],4);
    if smalloff_flag > 1
        figure, imshow(I);
        if smalloff_flag > 2
            saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-binary_WOSmall.tif',j)));
        end
    end
    clear cell_sizes  max_cell_size min_cell_size
else
    warning("Very small cells are usually the result of issues during the segmentation process.")
end

%% find label and visualize all connected components = cells
if cc_flag > 0
    cc = bwconncomp(I,4); % Using connectivity 4, we can have 1px borders
    labeled_cells = labelmatrix(cc);
    if cc_flag > 1
        cc_RGBlabel = label2rgb(labeled_cells, @spring, 'c', 'shuffle');
        figure, imshow(cc_RGBlabel);
        title('connected components');
        if cc_flag > 2
            saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-connected_components.tif',j)));
        end
    end
end

%% extract geometrical data from cells and applying the scaling if required
if celldata_flag > 0
    cell_data = regionprops(cc,{...
        'Area',...
        'Perimeter',...
        'Centroid',...
        'MajorAxisLength',...
        'MinorAxisLength',...
        'Orientation'});
    if celldata_flag > 1
        switch lower(orientation_plottype)
            case 'colormap'
                integer_orientation = zeros(cc.ImageSize);
                integer_orientation(labeled_cells>0) = ceil(abs([cell_data(labeled_cells(labeled_cells>0)).Orientation]));
                integer_orientationcc_RGBlabel = label2rgb(integer_orientation, parula(90));
                figure, imshow(integer_orientationcc_RGBlabel);
                caxis([0, 90]);
                colorbar;
                title('Cell orientation (°)');
            case 'ellipse'
                angle_distribution = linspace(0,2*pi,50);
                cellprops = mat2cell( [ ...
                    cell_data.Orientation; ...
                    cell_data.MajorAxisLength; ...
                    cell_data.MinorAxisLength; ...
                    vertcat(cell_data.Centroid)'; ...
                    ]',ones(cc.NumObjects,1),[5]);
                ellipse_calculator = @(cell_properties) [
                    cell_properties(4) + cell_properties(2)*cos(angle_distribution)*cosd(-cell_properties(1))/2 - cell_properties(3)*sin(angle_distribution)*sind(-cell_properties(1))/2;
                    cell_properties(5) + cell_properties(2)*cos(angle_distribution)*sind(-cell_properties(1))/2 + cell_properties(3)*sin(angle_distribution)*cosd(-cell_properties(1))/2]';
                ellipses = cellfun(ellipse_calculator,cellprops,'UniformOutput',false);
                figure, imshow(I);
                hold on 
                plot(getcolumn([ellipses{:}],1:2:2*cc.NumObjects),getcolumn([ellipses{:}],2:2:2*cc.NumObjects),'b','Linewidth',2);
            case 'major_axis'
                cellprops = mat2cell( [ ...
                    cell_data.Orientation; ...
                    cell_data.MajorAxisLength; ...
                    cell_data.MinorAxisLength; ...
                    vertcat(cell_data.Centroid)'; ...
                    ]',ones(cc.NumObjects,1),[5]);
                axis_calculator = @(cell_properties) [
                    cell_properties(4) + 0.8*linspace(-cell_properties(2)/2,cell_properties(2)/2,3)*cosd(-cell_properties(1)); 
                    cell_properties(5) + 0.8*linspace(-cell_properties(2)/2,cell_properties(2)/2,3)*sind(-cell_properties(1))]';
                axis = cellfun(axis_calculator,cellprops,'UniformOutput',false);
                figure, ax = axes; imshow(I);
                hold on
                plot(getcolumn([axis{:}],1:2:2*cc.NumObjects),getcolumn([axis{:}],2:2:2*cc.NumObjects),'b','Linewidth',2);
            otherwise
                disp('This option is not available and it will not be plotted or saved')
        end
        
        if celldata_flag > 2
            saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-cell_orientation.tif',j)));
        end
    end
end

%% Find adjacency matrix and graph
if graph_flag > 0
    dilated_labeled_cells = imdilate(labeled_cells,dilate_strel);
%     cell_edges = imLabelEdges(dilated_labeled_cells);
%     cell_edges_RGBlabel = label2rgb(cell_edges,parula(max(cell_edges,[],'all')));
%     figure, imshow(cell_edges_RGBlabel);
%     title('Edges');
    
    adjacency_matrix = zeros(max(dilated_labeled_cells,[],'all'));
    for cell=1:cc.NumObjects
        overlapping_pixels = immultiply(dilated_labeled_cells, imdilate(dilated_labeled_cells==cell, dilate_strel));
        nz_overlapping_pixels = overlapping_pixels(overlapping_pixels~=0);
        [GC,overlapping_cells] = groupcounts(nz_overlapping_pixels); % Removing contacts between cells that share very little space
        gc_min = minimum_contact_length(mean([cell_data(cell).Perimeter]));
        overlapping_cells = overlapping_cells(GC>gc_min);
        adjacency_matrix(cell,overlapping_cells) = 1;
    end
    adjacency_matrix = adjacency_matrix | adjacency_matrix';
    
    g = graph(adjacency_matrix,'omitselfloops');

    cell_centroids = vertcat(cell_data.Centroid);
    cells_xc = cell_centroids(:,1);
    cells_yc = cell_centroids(:,2);

    if graph_flag > 1
        figure, imshow(I);
        hold on
        plot(g,'XData',cells_xc,'YData',cells_yc) 
        if graph_flag > 2
            saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-adjacency_graph.tif',j)));
        end
    end
end

%% Find cycles that have a certain topological charge
if topocycles_flag > 0
    [contact_cycles,edgecycles] = allcycles(g,'MaxCycleLength',max_cycle_elements,'MinCycleLength',min_cycle_elements);

    cell_centroids_cell = {cell_data.Centroid}';
    cycle_centroids_cell = cell_centroids_cell(cell2mat(contact_cycles)); 
    cycle_centroids_cell = mat2cell(cell2mat(cycle_centroids_cell),ones(size(cycle_centroids_cell,1),1),2*max_cycle_elements);

    % remove complex polygons
    if complexpoloff_flag > 0
        complex_pol_index = @(cycle_xy) polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).NumRegions;
        complex_pol_indices = cellfun(complex_pol_index,cycle_centroids_cell);
        contact_cycles = contact_cycles(complex_pol_indices == 1);
        edgecycles = edgecycles(complex_pol_indices == 1);
        cycle_centroids_cell = cycle_centroids_cell(complex_pol_indices == 1);
    end

    % all cycles --> counter-clockwise
    clockwise_index = @(cycle_xy) ispolycw(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)]));
    clockwise_indices = cellfun(clockwise_index,cycle_centroids_cell);
    ccw_cycle = @(cycle, cycle_xy) [flip(cycle)*clockwise_index(cycle_xy) + cycle*(1-clockwise_index(cycle_xy))];
    contact_cycles = cellfun(ccw_cycle,contact_cycles,cycle_centroids_cell,'UniformOutput',false);
    contact_edgecycles = cellfun(ccw_cycle,edgecycles,cycle_centroids_cell,'UniformOutput',false);

    % Solidity and roundness filters
    if (solidityfilter_flag > 0 || roundnessfilter_flag > 0)
        % Solidity
        solidity = @(cycle_xy)...
            polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).area /...
            polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).convhull.area;
        solidities = cellfun(solidity,cycle_centroids_cell);
        % Roundness
        roundness = @(cycle_xy)...
            4*pi*polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).area /...
            (polyshape(cycle_xy([1:2:length(cycle_xy)]),cycle_xy([2:2:length(cycle_xy)])).perimeter)^2;
        roundnesses = cellfun(roundness,cycle_centroids_cell);

        contact_cycles = contact_cycles(logical((solidities > (0.9*logical(solidityfilter_flag))).*(roundnesses > (0.8*logical(roundnessfilter_flag)))));
        contact_edgecycles = contact_edgecycles(logical((solidities > (0.9*logical(solidityfilter_flag))).*(roundnesses > (0.8*logical(roundnessfilter_flag)))));
    end

    topo_wrapper = @(cell_cycle) topological_charge(cell_cycle, cell_data);

    topologicalCharges = cellfun(topo_wrapper,contact_cycles);

    if (plusonedefs_flag > 1 || minusonedefs_flag > 1 || plushalfdefs_flag > 1 || minushalfdefs_flag > 1)
        tiledfig = figure; % tiledlayout flow with cycles
        singlefig = figure; % superimposed defect centers
        imshow(I);
    end

    if plusonedefs_flag > 0
        plus_one_defs = find(topologicalCharges==1);
        plusone_x = mean(cells_xc(vertcat(contact_cycles{plus_one_defs}))');
        plusone_y = mean(cells_yc(vertcat(contact_cycles{plus_one_defs}))');
        if plusonedefs_flag > 1
            % defect cycles plot
            cmap3 = parula(length(plus_one_defs));
            cmap3 = cmap3(randperm(size(cmap3, 1)), :);
            figure(tiledfig)
            nexttile, imshow(I);
            hold on
            for i = 1:length(plus_one_defs)
                highlight(plot(g,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{plus_one_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
            end
            hold off
            title("Defects: +1")
            % defect centers plot
            figure(singlefig)
            hold on
            plot(plusone_x,plusone_y,'Marker','o','MarkerSize',10,'MarkerFaceColor','#EDB120','MarkerEdgeColor','#EDB120','LineStyle','none')
            hold off
        end
    end

    if minusonedefs_flag > 0
        minus_one_defs = find(topologicalCharges==-1);
        minusone_x = mean(cells_xc(vertcat(contact_cycles{minus_one_defs}))');
        minusone_y = mean(cells_yc(vertcat(contact_cycles{minus_one_defs}))');
        if minusonedefs_flag > 1
            % defect cycles plot
            cmap3 = parula(length(minus_one_defs));
            cmap3 = cmap3(randperm(size(cmap3, 1)), :);
            figure(tiledfig)
            nexttile, imshow(I);
            hold on
            for i = 1:length(minus_one_defs)
                highlight(plot(g,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{minus_one_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
            end
            hold off
            title("Defects: -1")
            % defect centers plot
            figure(singlefig)
            hold on
            plot(minusone_x,minusone_y,'Marker','o','MarkerSize',10,'MarkerFaceColor','#77AC30','MarkerEdgeColor','#77AC30','LineStyle','none')
            hold off
        end
    end

    if plushalfdefs_flag > 0
        plus_half_defs = find(topologicalCharges==0.5);
        plushalf_x = mean(cells_xc(vertcat(contact_cycles{plus_half_defs}))');
        plushalf_y = mean(cells_yc(vertcat(contact_cycles{plus_half_defs}))');
        if plushalfdefs_flag > 1
            % defect cycles plot
            cmap3 = parula(length(plus_half_defs));
            cmap3 = cmap3(randperm(size(cmap3, 1)), :);
            figure(tiledfig)
            nexttile, imshow(I);
            hold on
            for i = 1:length(plus_half_defs)
                highlight(plot(g,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{plus_half_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
            end
            hold off
            title("Defects: +1/2")
            % defect centers plot
            figure(singlefig)
            hold on
            plot(plushalf_x,plushalf_y,'Marker','o','MarkerSize',10,'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','LineStyle','none')
            hold off
        end
    end

    if minushalfdefs_flag > 0
        minus_half_defs = find(topologicalCharges==-0.5);
        minushalf_x = mean(cells_xc(vertcat(contact_cycles{minus_half_defs}))');
        minushalf_y = mean(cells_yc(vertcat(contact_cycles{minus_half_defs}))');
        if minushalfdefs_flag > 1
            % defect cycles plot
            cmap3 = parula(length(minus_half_defs));
            cmap3 = cmap3(randperm(size(cmap3, 1)), :);
            figure(tiledfig)
            nexttile, imshow(I);
            hold on
            for i = 1:length(minus_half_defs)
                highlight(plot(g,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{minus_half_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
            end
            hold off
            title("Defects: -1/2")
            % defect centers plot
            figure(singlefig)
            hold on
            plot(minushalf_x,minushalf_y,'Marker','o','MarkerSize',10,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E','LineStyle','none')
            hold off
        end
    end

    if (plusonedefs_flag > 2 || minusonedefs_flag > 2 || plushalfdefs_flag > 2 || minushalfdefs_flag > 2)
        saveas(tiledfig,fullfile(analysis_foldername, sprintf('00%d-defect_cycles.tif',j)));
        saveas(singlefig,fullfile(analysis_foldername, sprintf('00%d-defect_centers.tif',j)));
    end

end

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

clear *_flag
% warning('off')
% warning('on')
