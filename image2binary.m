function [I, fig] = image2binary(image2binary_flag, origin_filename, cell_color)

global I analysis_foldername j

    fig = nan;
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
            fig = figure;
            imshow(I);
            if image2binary_flag > 2
                saveas(gcf,fullfile(analysis_foldername, sprintf('00%d-binary.tif',j)));
            end
        end
    else
        error('You cannot perform an analysis without a binary image. Set image2binary_flag to 1 or a higher value.')
    end
end