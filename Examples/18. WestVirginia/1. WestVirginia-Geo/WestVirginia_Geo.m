%% Geometry

clear;

load usastates.mat
lat = usastates(46).Lat(1:end-3);
lon = usastates(46).Lon(1:end-3);


theta = mean(lat);    % Central latitude
x = lon * cosd(theta);% Scale longitude
y = lat;              % Latitude remains unchanged

nodes = [x;y];
n= size(nodes,2);
edges =[1:n; 2:n,1];

DT = delaunayTriangulation(nodes', edges');
TF = isInterior(DT);
elements=DT(TF,:)';



TR = triangulation(elements',nodes');
stlwrite(TR,'Tennessee.stl');
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

exportgraphics(gca,'WestVirginia-Geo.pdf','ContentType','vector');
savefig(gcf,'WestVirginia-Geo.fig');
