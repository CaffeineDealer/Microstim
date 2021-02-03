function [qs cmap] = qsign(x,y,type)

% qsign function provides color code to relate scatter plot to RF shift

% Inputs: 
% x,y: arrow coordinates for arrow or scatter plot axis for scatter
% type: scatter or arrow

% Outputs:
% qs: position relative to unity line or quadrant position
% cmap: color map

% Written by Yavar Korkian
% Created on 28.10.2020

switch type
    case 'scatter'
        if x < y
            qs = 1;
            cmap = [0 0 1];
        elseif x > y
            qs = -1;
            cmap = [1 0 0];
        elseif x == y
            qs = 0;
            cmap = [0 1 0];
        end
        
    case 'arrow'
        if x>0 && y>0
            qs = 1;
            cmap = [0 0 1];
        elseif x<0 && y>0
            qs = 2;
            cmap = [0 1 0];
        elseif x<0 && y<0
            qs = 3;
            cmap = [1 0 0];
        elseif x>0 && y<0
            qs = 4;
            cmap = [0 0 0];
        end
end