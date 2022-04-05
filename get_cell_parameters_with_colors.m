close all;
clear all; 
clc;

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
% 
%% erase cells at the borders 
dilatedLabel = imclearborder(I);
% figure, imshow(dilatedLabel);

CC = bwconncomp(dilatedLabel);
labeled = labelmatrix(CC);
RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
% figure, imshow(RGB_label);
% % title('connected components');


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

se = offsetstrel('ball',3,3);
dilateLabel = imdilate(labeled,se);
LBL = imLabelEdges(dilateLabel);
% Show them
cmap2 = parula(max(max(LBL)));
LBLrgb = label2rgb(LBL,cmap2); %build rgb image
% figure;
% imshow(LBLrgb);
% title('Edges');

figure
tiledlayout flow
% imshow(K);
% hold on

adjacencyMatrix = zeros(max(max(dilateLabel)));
N = double(max(dilateLabel(:)));
se = strel('diamond', 2);
for c=1:N
    lbl = intersect(1:N, unique(immultiply(dilateLabel, imdilate(dilateLabel==c, se))));
    lbl = lbl(lbl~=c);
    adjacencyMatrix(c,lbl) = 1;
end

g = graph(adjacencyMatrix);
xydata = vertcat(data.Centroid);
xdata = xydata(:,1);
ydata = xydata(:,2);
%  p = plot(g,'XData',xdata,'YData',ydata);

[cycles,edgecycles] = allcycles(g,'MaxCycleLength',6,'MinCycleLength',6);

topo_wrapper = @(cell_cycle) topological_charge(cell_cycle, data);

topologicalCharges = cellfun(topo_wrapper,cycles);
plusOneDefs = find(topologicalCharges==0.5);

cmap3 = parula(length(plusOneDefs));
cmap3 = cmap3(randperm(size(cmap3, 1)), :);
for i = 181:182 %1:length(plusOneDefs)
    nexttile
    imshow(K)
    hold on
    highlight(plot(g,'XData',xdata,'YData',ydata),'Edges',edgecycles{plusOneDefs(i)},'EdgeColor',cmap3(i,:),'LineWidth',1.5,'NodeColor',cmap3(i,:),'MarkerSize',6)
    title("Defect " + i)
    labels = {'1','2','3','4','5','6'};
    plot(xdata(cycles{plusOneDefs(i)}),ydata(cycles{plusOneDefs(i)}),'Marker','.','MarkerSize',6,'MarkerFaceColor','r')
    text(xdata(cycles{plusOneDefs(i)}),ydata(cycles{plusOneDefs(i)}),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
    xmean = mean(xdata(cycles{plusOneDefs(i)}));
    ymean = mean(ydata(cycles{plusOneDefs(i)}));
    plot(xmean,ymean,'Marker','*','MarkerSize',6,'MarkerFaceColor','r')
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

% % Orientation vector visualization superimposed to original image
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

%% Adjacent cells calculation

figure;
imshow(K);
hold on

C = vertcat(data.Centroid);
DT = delaunay(C(:,1), C(:,2))
triplot(DT,C(:,1), C(:,2));

% pixelinlinemat = zeros(length(data), length(data),100,2);
% 
% for k = 1:length(data) 
%     for h = 1:length(data)
%         cellmat(h,k,:) = data(k).Centroid - data(h).Centroid;
%         celldistmat(h,k) = norm(reshape(cellmat(h,k,:),[1,2]))*(norm(reshape(cellmat(h,k,:),[1,2]))<100) ;
%         cellanglemat(h,k) = atan2(cellmat(h,k,2),cellmat(h,k,1))*(norm(reshape(cellmat(h,k,:),[1,2]))<100);
%         dispmat(h,k,:) = [cos(cellanglemat(h,k)) sin(cellanglemat(h,k))];
%         for i = 1:ceil(celldistmat(h,k))
%             pixelinlinemat(h,k,i,:) = round(data(h).Centroid + i*reshape(dispmat(h,k,:),[1,2]));
%         end
%         if h < 203 && h > 201
%             xu = reshape(pixelinlinemat(h,k,:,1), [1,100]);
%             yu = reshape(pixelinlinemat(h,k,:,2), [1,100]);
%             isNZ=(~reshape(pixelinlinemat(h,k,:,1), [1,100])==0);
%             xu = xu(isNZ);
%             yu = yu(isNZ);
%             plot(xu,yu,'r','Linewidth',2)
%         end
%     end
% end
% 
% 
% % plot([data(1).Centroid(1) data(2).Centroid(1)], ...
% %     [data(1).Centroid(2) data(2).Centroid(2)],'g','Linewidth',2)


%% Saving pictures

name_picture=sprintf('00%d-Colors.tif',j);
saveas(gcf,name_picture);


%% calculations Shape Index
table.SI = 4 *pi* table.Area ./ (table.Perimeter.^2);
table.SIjamming = table.Perimeter./sqrt(table.Area);


table_name=sprintf('00%d-CY5.csv',j);
% writetable(table,table_name);

function topoCharge = topological_charge(cycle, data_struct)
% For clockwise cycles  
    xydata = vertcat(data_struct.Centroid);
    xdata = xydata(:,1);
    ydata = xydata(:,2);
    x = xdata(cycle);
    y = ydata(cycle);
    ch = convhull(x,y,'Simplify',true);
    ch = ch(1:end-1);
    % plot(x(ch),y(ch),'g')
    oDef = [data_struct(cycle).Orientation];
    oDef_ch = oDef(ch);
    oDef_ch_aux = [oDef_ch(end) oDef_ch(1:end-1)];
    difODef_ch = (oDef_ch_aux - oDef_ch)/180*pi;
    difODef_ch_norm =  [difODef_ch(abs(difODef_ch)<=pi/2), difODef_ch(difODef_ch<-pi/2) + pi, difODef_ch(difODef_ch>pi/2) - pi];
    topoCharge = sum(difODef_ch_norm)/(2*pi);
end

function topoCharge = topological_charge_naive(cycle, data_struct)
    % For clockwise cycles
    oDef = [data_struct(cycle).Orientation];
    oDef_aux = [oDef(end) oDef(1:end-1)];
    difODef = (oDef_aux - oDef)/180*pi;
    difODef_norm =  [difODef(abs(difODef)<=pi/2), difODef(difODef<-pi/2) + pi, difODef(difODef>pi/2) - pi];
    topoCharge = sum(difODef_norm)/(2*pi);
end

