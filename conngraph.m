function [connectivity_graph, cells_xc, cells_yc,fig] = conngraph(graph_flag, cc, labeled_cells, dilate_strel, cell_data, minimum_contact_length)

global I analysis_foldername j

    fig = nan;
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
        adjacency_matrix = adjacency_matrix | adjacency_matrix'; % Solving symmetry conflicts due to dilate_strel
        
        connectivity_graph = graph(adjacency_matrix,'omitselfloops');
    
        cell_centroids = vertcat(cell_data.Centroid);
        cells_xc = cell_centroids(:,1);
        cells_yc = cell_centroids(:,2);
    
        if graph_flag > 1
            fig = figure;
            imshow(I);
            hold on
            plot(connectivity_graph,'XData',cells_xc,'YData',cells_yc) 
            if graph_flag > 2
                saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-adjacency_graph.tif',j)));
            end
        end
    end
end