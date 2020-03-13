function plotsolution( grid, u )
%-----------------------------------------------------------
% Plot the temperature distribution in the fin
%
% INPUTS:   grid : one of the grids {coarse, medium, fine}
%           u    : solution vector
%
%-----------------------------------------------------------

clf;
axis off;
hold on;

caxis([min(u) max(u)]);
for i=1:5
   for j=grid.theta{i}'
      fill(grid.coor(j,1),grid.coor(j,2),u(j));
   end
end

title('Temperature distribution')
colorbar;

shading interp;

hold off;

