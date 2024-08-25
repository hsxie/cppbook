function [up,down] = envelope(x,y,interpMethod)

%ENVELOPE gets the data of upper and down envelope of the known input (x,y).
%
% Input parameters:
% x the abscissa of the given data
% y the ordinate of the given data
% interpMethod the interpolation method
%
% Output parameters:
% up the upper envelope, which has the same length as x.
% down the down envelope, which has the same length as x.
%
% See also DIFF INTERP1

% Designed by: Lei Wang, <WangLeiBox@hotmail.com>, 11-Mar-2003.
% Last Revision: 21-Mar-2003.
% Dept. Mechanical & Aerospace Engineering, NC State University.
% $Revision: 1.1 $ $Date: 3/21/2003 10:33 AM $

if length(x) ~= length(y)
error('Two input data should have the same length.');
end

if (nargin < 2)|(nargin > 3),
error('Please see help for INPUT DATA.');
elseif (nargin == 2)
interpMethod = 'linear';
end



% Find the extreme maxim values 
% and the corresponding indexes
%----------------------------------------------------
extrMaxValue = y(find(diff(sign(diff(y)))==-2)+1);
extrMaxIndex = find(diff(sign(diff(y)))==-2)+1;



% Find the extreme minim values 
% and the corresponding indexes
%----------------------------------------------------
extrMinValue = y(find(diff(sign(diff(y)))==+2)+1);
extrMinIndex = find(diff(sign(diff(y)))==+2)+1;



up = extrMaxValue;
up_x = x(extrMaxIndex);

down = extrMinValue;
down_x = x(extrMinIndex);



% Interpolation of the upper/down envelope data
%----------------------------------------------------
up = interp1(up_x,up,x,interpMethod); 
down = interp1(down_x,down,x,interpMethod);