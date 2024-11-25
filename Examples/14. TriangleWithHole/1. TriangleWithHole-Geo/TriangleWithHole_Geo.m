%% Geometry

clear;
x1=[-3    1   1 ];
y1=[0   -1   1  ];

N=50;r=0.65;
theta=-2*pi*(0:(N-1))/N;
x2=r*cos(theta);
y2=r*sin(theta);

nodes = [x1 x2; y1 y2];

edges =[1:3 4:53; 2:3,1 5:53, 4];

DT = delaunayTriangulation(nodes', edges');
TF = isInterior(DT);
elements=DT(TF,:)';

TR = triangulation(elements',nodes');
stlwrite(TR,'TriangleWithHole.stl');
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

exportgraphics(gca,'TriangleWithHole-Geo.pdf','ContentType','vector');
savefig(gcf,'TriangleWithHole-Geo.fig','compact');
