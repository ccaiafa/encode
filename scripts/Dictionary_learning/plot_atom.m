function [  ] = plot_atom( v, bvecs, color )

%v = v - min(v)*ones(size(v));
v = v/norm(v);
points = bvecs.*repmat(v,1,3);
scatter3(points(:,1),points(:,2),points(:,3),'.','MarkerFaceColor',color);
hold on

% [COEFF, SCORE, LATENT]=pca(points');
% plot3([0 SCORE(1,1)],[0 SCORE(2,1)],[0 SCORE(3,1)],'Color',color)
% plot3([0 SCORE(1,2)],[0 SCORE(2,2)],[0 SCORE(3,2)],'Color',color)

end

