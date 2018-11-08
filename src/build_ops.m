function [D,G,I]=build_ops(Grid)
   % author: Kristopher Darnell
   % date: 02/03/15
   % Description:
   % This function computes the discrete divergence and gradient matrices on a
   % regular staggered grid using central difference approximations. The
   % discrete gradient assumes homogeneous boundary conditions.
   % Input:
   % Grid = structure containing all pertinent information about the grid.
   % Output:
   vec = ones(Grid.Nx+1,1);
   D = (1/Grid.dx)*spdiags([-vec vec],0:1,Grid.Nx,Grid.Nx+1);
   Xf = spdiags(Grid.xf.^(Grid.Geom-1),0,Grid.Nx+1,Grid.Nx+1);
   Xcinv = spdiags(Grid.xc.^(1-Grid.Geom),0,Grid.Nx,Grid.Nx);
   G = -D'; 
   G(1,:) = zeros(1,Grid.Nx); 
   G(end,:) = G(1,:);
   D = Xcinv*D*Xf;
   I = speye(Grid.Nx);
   % Example call:
   % >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 5;
   % >> Grid = build_grid(Grid);
   % >> [D,G,I]=build_ops(Grid);