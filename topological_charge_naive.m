function topoCharge = topological_charge_naive(cycle, data_struct)
    % For clockwise cycles
    oDef = [data_struct(cycle).Orientation];
    oDef_aux = [oDef(end) oDef(1:end-1)];
    difODef = (oDef_aux - oDef)/180*pi;
    difODef_norm =  [difODef(abs(difODef)<=pi/2), difODef(difODef<-pi/2) + pi, difODef(difODef>pi/2) - pi];
    topoCharge = sum(difODef_norm)/(2*pi);
end