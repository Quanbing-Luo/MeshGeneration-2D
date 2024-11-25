%% Geometry generation with delaunayTriangulation and isInterior 

clear;
x1=[0.5, -0.5,  -0.5,  0.5 ];
y1=[0.5,  0.5,  -0.5, -0.5 ];
n1=size(x1,2);

x2=[0.4,  0.05,  0.05,  0.4];
y2=[0.4,  0.4 ,  0.05,  0.05];
n2=size(x2,2)+n1;

x3=[-0.4,  -0.05, -0.05,  -0.4];
y3=[0.4,   0.4,   0.05,   0.05];
n3=size(x3,2)+n2;

NN=25;r=0.15/2; 
dx=(0.5-2*r)*2/3+r;  dy=(0.55-2*r)/2+r;
theta=2*pi*(0:(NN-1))/NN;

x4= 0.5-dx+r*cos(theta);
y4=-0.5+dy+r*sin(theta);
n4=size(x4,2)+n3;

x5=-0.5+dx+r*cos(theta);
y5=-0.5+dy+r*sin(theta);
n5=size(x5,2)+n4;


nodes=[x1,x2,x3,x4,x5;y1,y2,y3,y4,y5];
edges=[1:n5; 2:n1,1, (n1+2):n2,(n1+1), (n2+2):n3, ...
    (n2+1), (n3+2):n4,(n3+1), (n4+2):n5,(n4+1)];

DT = delaunayTriangulation(nodes', edges');
TF = isInterior(DT);
elements=DT(TF,:)';

% TR = triangulation(elements',nodes');
% stlwrite(TR,'FourHoles.stl');
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

exportgraphics(gca,'FourHoles-Geo.pdf','ContentType','vector');
savefig(gcf,'FourHoles-Geo.fig','compact');
