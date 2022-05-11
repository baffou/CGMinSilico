classdef CGmatrix
% class of images that are either grating transmittance or E-field
% G. Baffou, CNRS, Jan 2022

    properties
        im
        pxSize %[m]
    end
    
    properties(Dependent)
        Npx    %[px]
        imSize %[m]
    end
    
    methods

        function M = CGmatrix(im0,imS)
            % im: image matrix
            % imS: image size in µm
            M.im = im0;
            M.pxSize = imS/M.Npx;
        end
        
        function M = set.im(M,im0)
            if size(im0,1)==size(im0,2)
                M.im = im0;
            else
                error('The matrix must be square')
            end
        end
        
        function val = get.Npx(M)
            val = size(M.im,1);
        end
        
        function val = get.imSize(M)
            val = M.pxSize*M.Npx;
        end
        
        function M = removeBorderLines(M,n)
            % removes pixels at the boundary of the image. Usefule when
            % propagation of light blurries the boundary of the
            % calculated interferogram
            M.im = M.im(1:end-n,1:end-n);
        end
        
        function M = redimension(M,camPxSize,pxRatio)
            % Resize the image, ie shrink Nx by a factor of pxRatio
            % and set the pixel size in [m]
            fac = pxRatio;
            M.im = imresize(M.im,fac);
            M.pxSize = camPxSize;
        end
        
        function [M,fac] = setI0(M,I0)
            % sets the maximum intensity of the image
            fac = 1/max(abs(M.im(:)))*I0;
            M.im = M.im*fac;
        end
        
        function [M,fac] = timesC(M,fac)
            % Multiply the image values by a factor
            M.im = M.im*fac;
        end
        
        function M = tile(M,N)
            % Tile a unit cell
            % N: new matrix size
            factor=ceil(N/M.Npx);
            M.im=repmat(M.im,factor,factor);
            M.im=M.im(1:N,1:N); % in case gratingSize is not an integer, it was overestimated in the previous line so we cut here the image properly à Ni
        end
        
        function M=square(M)
            % Compute the intensity matrix from the amplitude matrix
            M.im=abs(M.im).^2;
        end
        
        function M=propagation(M,lambda,L,dx,dy)
            % Propagate the electric field over a distance L
            % lambda: wavelength of propagation
            % L: distance of propagation (positive or negative)
            % dx, dy: shift in x and y, due to an titled illumination
            if nargin==3
                dx = 0;
                dy = 0;
            end
            M.im = improp(M.im,M.imSize,lambda,L,dx,dy);
        end
        
        function M = TileRot5(M,th)
            % Tile 5 unit cells into a super unit cell and rotate by an
            % angle theta.
            fac = 5;
            if nargin==1
                theta = acos(3/5);
            elseif nargin==2
                theta = th;
                warning('theta must equal acos(3/5) or 0 in this code version')
            else
                error('wrong number of inputs')
            end
            [X,Y] = meshgrid(1:M.Npx*fac,1:M.Npx*fac);
            M.im = rotation(X,Y,M.im,theta);
        end
        
        function figure(M,h)
            % Display the matrix
            % h: pre-existing figure handle, if a figure already exists.
            if nargin==2
                figure(h)
            else
                figure
            end
            ai = imagesc(1e6*(1:M.Npx)*M.pxSize,1e6*(1:M.Npx)*M.pxSize,real(M.im));
            set(ai.Parent,'DataAspectRatio',[1 1 1])
            xlabel(gca,'um')
            colorbar
            colormap(gca,CGcolorScale)
        end

        function obj=times(obj1,obj2)
            % Overload of the times methods. Enable the multiplication of
            % two matrices.
            if obj1.Npx~=obj2.Npx
                error('Different image sizes')
            elseif obj1.pxSize~=obj2.pxSize
                error('Different pixel sizes')
            end
            obj = obj1;
            obj.im = obj1.im.*obj2.im;
        end
    end
    
end












