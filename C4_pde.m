% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
% 1) Export the required variables from pdetool and create a MATLAB script
%    to perform operations on these.
% 2) Define the problem completely using a MATLAB script. See
%    https://www.mathworks.com/help/pde/examples.html for examples
%    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',2);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1 0.87499999999999978 1]);
set(ax,'PlotBoxAspectRatio',[1.4999999999999998 1 1.2499999999999996]);
set(ax,'XLim',[-0.20000000000000001 2.2000000000000002]);
set(ax,'YLim',[-0.20000000000000001 1.2]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');
pdetool('gridon','on');

% Geometry description:
pderect([0 2 1 0],'R1');
pdecirc(0.5,0.5,0.10000000000000009,'C1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1-C1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(8,...
'dir',...
2,...
char('1','0','0','1'),...
char('1','0'))
pdesetbd(7,...
'dir',...
2,...
char('1','0','0','1'),...
char('1','0'))
pdesetbd(6,...
'dir',...
2,...
char('1','0','0','1'),...
char('1','0'))
pdesetbd(5,...
'dir',...
2,...
char('1','0','0','1'),...
char('1','0'))
pdesetbd(4,...
'dir',...
2,...
char('1','0','0','1'),...
char('1','0'))
pdesetbd(3,...
'dir',...
2,...
char('1','0','0','1'),...
char('1','0'))
pdesetbd(2,...
'dir',...
2,...
char('1','0','0','1'),...
char('0','0'))
pdesetbd(1,...
'neu',...
2,...
char('1','0','0','1'),...
char('0','0'))

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('refine')

% PDE coefficients:
pdeseteq(2,...
'0',...
char('uy(1,:)','0.0','0.0','uy(2,:)'),...
char('-ux(1,:)','-ux(2,:)'),...
char('1.0','0.0','0.0','1.0'),...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['0       ';...
'0.0     ';...
'uy(1,:) ';...
'0.0     ';...
'-ux(1,:)';...
'1.0     ';...
'0.0     ';...
'0.0     ';...
'0       ';...
'0.0     ';...
'uy(2,:) ';...
'-ux(2,:)';...
'0.0     ';...
'1.0     '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1500','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 7 1 0 0 0 11 1 1 0 1 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')
