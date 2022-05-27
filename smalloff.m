function [I, fig] = smalloff(smalloff_flag) 

global I analysis_foldername j

    fig = nan;
    if smalloff_flag > 0
        cell_sizes = regionprops(I, 'Area');
        max_cell_size = max([cell_sizes.Area]);
        min_cell_size = mean([cell_sizes.Area])/10; 
        I = bwareafilt(I,[min_cell_size max_cell_size],4);
        if smalloff_flag > 1
            fig = figure;
            imshow(I);
            if smalloff_flag > 2
                saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-binary_WOSmall.tif',j)));
            end
        end
        % clear cell_sizes  max_cell_size min_cell_size
    else
        warning("Very small cells are usually the result of issues during the segmentation process.")
    end
end    