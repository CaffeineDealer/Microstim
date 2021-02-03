classdef fitVonMises < handle
    
    properties (Access = public)
        Out
        x
    end
    
    
    methods
         function y = fitVonMises(x,name)
             data = load(name,'npd');
             x = x + 2;
             y.x = x;
             y.Out = data;     
         end
         
    end
end