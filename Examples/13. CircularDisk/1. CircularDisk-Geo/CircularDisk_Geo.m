%% Geometry

clear;
n=80;
theta=2*pi*(0:(n-1))/n;
nodes = [cos(theta);sin(theta)];
edges =[1:n; 2:n,1];

DT = delaunayTriangulation(nodes', edges');
TF = isInterior(DT);
elements=DT(TF,:)';


TR = triangulation(elements',nodes');
stlwrite(TR,'CircularDisk.stl');
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

exportgraphics(gca,'CircularDisk-Geo.pdf','ContentType','vector');
savefig(gcf,'CircularDisk-Geo.fig','compact');
