function xa = xaxis(arg)
%XAXIS manipulate x-axis of plot
%
% 	V = XAXIS returns a row vector containing the scaling of the
%   x-axis of the current plot.
%
% 	XAXIS([XMIN XMAX]) sets scaling for the x-axis on the current plot.

% Author: Ren\'e Laterveer,  3 dec. '92
% Code stolen form axis.m in matlab toolbox
% 26/2/93 - added support for axes increasing from top to bottom

ax = gca;

if (nargin == 0)
    xa = get(ax,'XLim');
else
    if ((max(size(arg)) == 2))
        if (arg(1)<arg(2))
          dir = 'normal';
        else
          dir = 'reverse';
          a = arg(1);
          b = arg(2);
          arg = [b a];
        end
        set(ax,'XLim',arg,...
               'XLimMode','manual',...
               'XDir',dir);
        if ~ishold
            view(2);
        end
    else
        error('Vector must have 2 elements.')
    end
end
