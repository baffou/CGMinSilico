classdef FcropParams
% parameters of the crop in the Fourier plane
% G. Baffou, CNRS, Jan 2022

properties
        x   % x-coordinate of the center of the crop
        y   % y-coordinate of the center of the crop
        R   % radius of the circular crop
        Npx % Npx*Npx : total number of the pixels of the image
    end
    properties(Dependent)
        shiftx
        shifty
        angle
    end
    
    methods
        
        function obj = FcropParams(x,y,R,Npx)
            if nargin==4
                obj.Npx = Npx;
                obj.x = x;
                obj.y = y;
                obj.R = R;
            end
        end
        
        function obj = set.x(obj,x0)
            if isnumeric(x0)
                if x0>=0 && x0<=obj.Npx
                    obj.x = x0;
                else
                    error(['wrong value for x, which equals ' num2str(x0) 'while Nx= ' num2str(obj.Npx) '.'])
                end
            else
                error('x must be a number')
            end
        end
        
        function obj = set.y(obj,y0)
            if isnumeric(y0)
                if y0>=0 && y0<=obj.Npx
                    obj.y = y0;
                else
                    error('wrong value for y')
                end
            else
                error('y must be a number')
            end
        end
        
        function val = get.shiftx(obj)
            val = round(obj.x-(obj.Npx/2+1));
        end
        
        function val = get.shifty(obj)
            val = round(obj.y-(obj.Npx/2+1));
        end
        
        function val = get.angle(obj)
            nshiftx = obj.shiftx/obj.Npx;
            nshifty = obj.shifty/obj.Npx;
            
            nshiftr = sqrt(nshiftx^2+nshifty^2);
            
            val.cos = nshiftx/nshiftr;
            val.sin = nshifty/nshiftr;
        end
        
        function obj2 = rotate90(obj)
            % Rotate a circular crop by 90° aroud the center of the image
            x0 = obj.Npx/2+1;
            y0 = obj.Npx/2+1;
            dx = (obj.x-x0);
            dy = obj.y-y0;
            
            x2 = x0-dy;
            y2 = y0+dx;
            
            obj2 = FcropParameters(x2,y2,obj.R,obj.Npx);
        end
        
        function obj2 = rotate180(obj)
            % Rotate a circular crop by 180° aroud the center of the image
            x0 = obj.Npx/2+1;
            y0 = obj.Npx/2+1;
            dx = obj.x-x0;
            dy = obj.y-y0;
            
            x2 = x0-dx;
            y2 = y0-dy;
            
            obj2 = FcropParameters(x2,y2,obj.R,obj.Npx);
        end
        
    end
end
