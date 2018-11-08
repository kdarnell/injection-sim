function [Grid] = build_grid(Grid,varargin)
   % author: Kristopher Darnell
   % date: 02/03/15
   % description:
   % This function computes takes in minimal definition of the computational
   % domain and grid and computes all containing all pertinent information
   % about the grid.
   % Input:
%    Grid.xmin = left boundary of the domain
%    Grid.xmax = right bondary of the domain
%    Grid.Nx   = number of grid cells
%
%
%    Output:
   coordoptions = {'cartesian','polar','spherical'};
   if isempty(varargin)
       Grid.Geom = 1; Grid.Geominfo = coordoptions{1};
   else
       Grid.Geominfo = varargin; Grid.Geom = find(strcmp(varargin,coordoptions));
   end
   Grid.Lx = Grid.xmax - Grid.xmin;
   Grid.dx = Grid.Lx / (Grid.Nx-1);
   Grid.xc = linspace(Grid.xmin,Grid.xmax,Grid.Nx)';
   Grid.xf = [Grid.xc - Grid.dx/2; Grid.xmax + Grid.dx/2];
   Grid.Nfx = length(Grid.xf);
   Grid.dof = ones(Grid.Nx,1);
   Grid.dof_xmin = [1;zeros(Grid.Nx-1,1)];  
   Grid.dof_xmax = [zeros(Grid.Nx-1,1);1];
  