function I = borderoff(borderoff_flag, I, analysis_foldername) 
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
end