classdef gaussDataTriangle < handle
    %Returns the gauss points and gauss weight. The values are in natural
    % co-ordinates.
    %Input: Order of integration.
    properties
        pt
        wt
        noGp
    end
    methods
        
        
        function dat=gaussDataTriangle(i)
            
            switch(i)
                case 1
                    dat.pt=[1/3,1/3];
                    dat.wt=0.5;
                case 2
                    dat.pt=[1/6,1/6;2/3,1/6;1/6,2/3];
                    dat.wt=[1/6;1/6;1/6];
                case 3
                    dat.pt=[1/3 1/3 ;3/5 1/5;1/5 3/5; 1/5 1/5];
                    dat.wt=[-9/32;25/96;25/96;25/96];
                case 4
                    dat.pt=[0 0 ;1/2  0;1  0; 1/2 1/2;0  1 ;0  1/2 ;1/3 1/3];
                    dat.wt=[1/40;1/15;1/40;1/15;1/40;1/15;9/40];
                    
                otherwise
                    msg = 'Invalid option: Order of integration greter than 4';
                    error(msg);
            end
            dat.noGp=size(dat.wt,1);
        end
        
    end
end