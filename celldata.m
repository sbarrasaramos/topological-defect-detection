function [cell_data, fig] = celldata(celldata_flag, micrometers, cc, labeled_cells)

global I analysis_foldername j

    fig = nan;
    if celldata_flag > 0
        cell_data = regionprops(cc,{...
            'Area',...
            'Perimeter',...
            'Centroid',...
            'MajorAxisLength',...
            'MinorAxisLength',...
            'Orientation'});
        switch lower(micrometers)
            case 'on' % scaling data
                scaling = mat2cell( [ ...
                    cell_data.Perimeter; ...
                    cell_data.MajorAxisLength; ...
                    cell_data.MinorAxisLength; ...
                    vertcat(cell_data.Centroid)'; ...
                    cell_data.Area ...
                    ]'*pixel_scale.*[ones(1,5) pixel_scale],ones(cc.NumObjects,1),[1 1 1 2 1]);
            [cell_data.Perimeter, cell_data.MajorAxisLength, cell_data.MinorAxisLength, cell_data.Centroid, cell_data.Area] = scaling{:};
            otherwise
                warning('Pixels are not necessarily comparable between units. This should not affect orientation or graphs')
        end
        if celldata_flag > 1
            integer_orientation = zeros(cc.ImageSize);
            integer_orientation(labeled_cells>0) = ceil(abs([cell_data(labeled_cells(labeled_cells>0)).Orientation]));
            integer_orientationcc_RGBlabel = label2rgb(integer_orientation, parula(90));
            fig = figure;
            imshow(integer_orientationcc_RGBlabel);
            caxis([0, 90]);
            colorbar;
            title('Cell orientation (°)');
            if celldata_flag > 2
                saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-cell_orientation.tif',j)));
            end
        end
    end
end