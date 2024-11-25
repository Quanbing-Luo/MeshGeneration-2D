%% Geometry generation with delaunayTriangulation and isInterior 

clear;
n=50;
theta=2*pi*((n-1):-1:0)/n;
nodes = [ cos(theta) -2  2  2 -2;
          sin(theta) -2 -2  2  2 ];
edges =[1:54; 2:50 1 52:54 51];

DT = delaunayTriangulation(nodes', edges');
TF = isInterior(DT);
elements=DT(TF,:)';

TR = triangulation(elements',nodes');
stlwrite(TR,'SquareWithHole.stl');
geometrywrite(nodes, elements, edges, 'geometry.dat');

%% Figure Trianglated Geometry

% [nodes, elements, edges] = geometryread('geometry.dat');

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

exportgraphics(gca,'SquareWithHole-Geo.pdf','ContentType','vector');
savefig(gcf,'SquareWithHole-Geo.fig','compact');
