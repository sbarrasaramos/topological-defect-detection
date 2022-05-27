function [cc, labeled_cells, fig] = connectedcomp(cc_flag)

global I analysis_foldername j

fig = nan;
    if cc_flag > 0
        cc = bwconncomp(I,4); % Using connectivity 4, we can have 1px borders
        labeled_cells = labelmatrix(cc);
        if cc_flag > 1
            cc_RGBlabel = label2rgb(labeled_cells, @spring, 'c', 'shuffle');
            fig = figure;
            imshow(cc_RGBlabel);
            title('connected components');
            if cc_flag > 2
                saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-connected_components.tif',j)));
            end
        end
    end
end