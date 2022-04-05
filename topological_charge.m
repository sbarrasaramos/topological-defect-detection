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
