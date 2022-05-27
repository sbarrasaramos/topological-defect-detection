function [ ...
    contact_cycles, edgecycles, topologicalCharges, tiledfig, singlefig, ...
    plus_one_defs, plusone_x, plusone_y, ...
    minus_one_defs, minusone_x, minusone_y, ...
    plus_half_defs, plushalf_x, plushalf_y, ...
    minus_half_defs, minushalf_x, minushalf_y ...
    ] = topocycles( ...
    topocycles_flag, complexpoloff_flag, solidityfilter_flag, roundnessfilter_flag, ...
    plusonedefs_flag, minusonedefs_flag, plushalfdefs_flag, minushalfdefs_flag, ...
    cell_data, connectivity_graph, max_cycle_elements, min_cycle_elements, cells_xc, cells_yc ...
)

global I analysis_foldername j

    if topocycles_flag > 0
        [contact_cycles,edgecycles] = allcycles(connectivity_graph,'MaxCycleLength',max_cycle_elements,'MinCycleLength',min_cycle_elements);
    
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
    
        tiledfig = nan;
        singlefig = nan;
    
        if (plusonedefs_flag > 0 || minusonedefs_flag > 0 || plushalfdefs_flag > 0 || minushalfdefs_flag > 0)
            tiledfig = figure; % tiledlayout flow with cycles
            singlefig = figure; % superimposed defect centers
            imshow(I);
        end
    
        plus_one_defs = nan;
        plusone_x = nan;
        plusone_y = nan;
    
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
                    highlight(plot(connectivity_graph,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{plus_one_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
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
    
        minus_one_defs = nan;
        minusone_x = nan;
        minusone_y = nan;

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
                    highlight(plot(connectivity_graph,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{minus_one_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
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
    
        plus_half_defs = nan;
        plushalf_x = nan;
        plushalf_y = nan;

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
                    highlight(plot(connectivity_graph,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{plus_half_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
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
    
        minus_half_defs = nan;
        minushalf_x = nan;
        minushalf_y = nan;

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
                    highlight(plot(connectivity_graph,'XData',cells_xc,'YData',cells_yc,'EdgeColor','#0072BD','NodeColor','#0072BD'),'Edges',contact_edgecycles{minus_half_defs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
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
