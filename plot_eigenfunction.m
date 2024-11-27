function [out] = plot_eigenfunction(V, F, evecs, evnum)
clf
figure
p=patch('Vertices', V, 'Faces', F', 'FaceVertexCData', evecs(:,evnum), 'LineStyle', 'none');
p.FaceColor='interp';
axis equal
colorbar
end