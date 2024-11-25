%% Read STL Geometry 

clear;
TR = stlread("PlateSquareHolePlanar.stl");
nodes = (TR.Points(:,1:2))';
elements = (TR.ConnectivityList)';
edges = freeBoundary(TR)';

% TR = triangulation(elements',nodes');
% stlwrite(TR,'PlateSquareHolePlanar.stl');
geometrywrite(nodes, elements, edges, 'geometry.dat');

%% Figure Trianglated Geometry

figure (1); clf; 
patch('Faces',elements','Vertices',nodes','FaceColor',[0.85, 0.85, 0.85], ...
    'EdgeColor','blue');

% hold on
% plot([nodes(1, edges(1,:)); nodes(1, edges(2,:))] , ...
%      [nodes(2, edges(1,:)); nodes(2, edges(2,:))], ...
%      Color = 'r', Marker='.',MarkerEdgeColor='b', MarkerSize = 8);
% hold on
% text(nodes(1,:),nodes(2,:), 'n' + string(1:size(nodes,2)));
% hold on
% text((nodes(1, edges(1,:)) + nodes(1, edges(2,:)))/2, ...
%      (nodes(2, edges(1,:)) + nodes(2, edges(2,:)))/2, ...
%     'e' + string(1:size(edges,2)));

axis equal off

exportgraphics(gca,'PlateSquareHolePlanar-Geo.pdf','ContentType','vector');
savefig(gcf,'PlateSquareHolePlanar-Geo.fig','compact');
