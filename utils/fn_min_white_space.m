function [ax] = fn_min_white_space(orig)
% Minimize the white space on a figure
% INPUTS:
%   orig [figure axis] - origianl figure axis
% OUTPUTS:
%   ax [figure axis] - adjusted figure axis

ax = orig;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

end