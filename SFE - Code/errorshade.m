function h = errorshade(x,y,sigma,color,varargin)
% errorshade plots a shaded region to indicate gaussian uncertainty. This function
% works by generating an RGB image of a specified color, and setting transparency 
% of the RGB image corresponding to uncertainty values. 
% 
%% Syntax 
% 
%  errorshade(x,y,sigma,color) 
%  errorshade(...,'resolution',res) 
%  errorshade(...,'average',N)
%  h = errorshade(...)
% 
%% Description 
% 
% errorshade(x,y,sigma,color) plots a gaussian shaded region centered about the 
% line given by x,y.  The input sigma represents one standard deviation of shading
% weight, and color is a three-element vector containing rgb values of the shading color.
%
% errorshade(...,'resolution',res) specifies resolution of the underlying RGB image. res
% can be a scalar value or a two-element vector in the form [xres yres].  Default resolution
% is 2500 pixels.  Larger values may take longer to plot, smaller values may appear jittery.
%
% errorshade(...,'average',N) specifies an N point moving average to smooth out local 
% spikes in data.  Use an odd-numbered integer value because even numbers will result in a 
% slight offset on the horizontal direction.  The averaging option requires the Image Processing 
% Toolbox. 
% 
% h = errorshade(...) returns a handle h of the plotted RGB image. 

%% Error checks: 
narginchk(4,inf) 
assert(numel(color)==3,'Input error: color must be three element vector.')
assert(numel(x)==numel(y),'Input error: Dimensions of x and y must agree.') 
assert(isscalar(sigma)==1,'Input error: sigma must be a scalar.') 
%% Input parsing 
tmp = strncmpi(varargin,'resolution',3); 
if any(tmp)
   res = varargin{find(tmp)+1}; 
   if isscalar(res) 
      res = [res res]; 
   else
      assert(numel(res)==2,'Input error: resolution must be a scalar or a two-element vector.') 
   end
else
   res = 2500*[1 1]; 
end
tmp = strncmpi(varargin,'average',2); 
if any(tmp)
   avg = varargin{find(tmp)+1}; 
   assert(isscalar(avg)==1,'Input error: moving average distance must be a scalar.')
   assert(license('test','image_toolbox')==1,'Cannot find a license for the Image Processing Toolbox, which is required for the moving average feature.') 
   y = imfilter(y(:),fspecial('average',[avg 1]),'replicate');
end
buffer = 3*sigma; % This is the padding to add around all measurements in the vertical dimension. 
%% Get limits: 
% Make a grid corresponding to the dimensions of the data +/- buffer: 
[X,Y] = meshgrid(linspace(min(x),max(x),res(1)),linspace(min(y)-buffer,max(y)+buffer,res(2))); 
% Find y locations along all x points of the grid
yi = interp1(x,y,linspace(min(x),max(x),res(1)));
% normal distribution: 
P = (1/sqrt(2*pi*sigma^2)) * exp(-(bsxfun(@minus,Y,yi)).^2/(2*sigma^2));
% Distance of each point from yi will be used as a measure of transparency:
Z = P-min(P(:));
Z = Z/max(Z(:)); 
% Create an RGB image of the specified color: 
RGB = cat(3,color(1)*ones(size(Z)),color(2)*ones(size(Z)),color(3)*ones(size(Z))); 
% Plot the RGB image of color: 
h = image(RGB,'xdata',X(1,:),'ydata',Y(:,1)); 
axis xy 
% Set transparency: 
set(h,'alphadata',Z)
% Set renderer to OpenGL because transparency only works with OpenGL: 
set(gcf,'renderer','OpenGL'); 
% Send to bottom: 
uistack(h,'bottom')
%% clean up: 
if nargout==0 
   clear h
end
end