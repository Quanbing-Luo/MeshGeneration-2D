%% Geometry generation with delaunayTriangulation and isInterior

clear;
r=0.45;
theta=pi:-pi/10:-pi/2;
%theta=pi:-pi/20:-pi/2;
%theta=pi:-pi/30:-pi/2;

nodes = [-0.5 -0.5 0.5 0.5  r*cos(theta),-r;  
         0.5 -0.5 -0.5 0.5  r*sin(theta),-r];

sz=length(theta); szt=4+sz+1;
edges =[1:szt;
  2:4,1,6:szt,5];

DT = delaunayTriangulation(nodes', edges');
TF = isInterior(DT);
elements=DT(TF,:)';

% TR = triangulation(elements',nodes');
% stlwrite(TR,'SquareLargeHole.stl');

%% Figure Trianglated Geometry

figure (1); clf;
patch('Faces',elements','Vertices',nodes','FaceColor',[0.85, 0.85, 0.85], ...
    'EdgeColor','b');
% hold on
% plot([nodes(1, edges(1,:));nodes(1, edges(2,:))] , ...
%      [nodes(2, edges(1,:)); nodes(2, edges(2,:))], ...
%      Color = 'r', Marker='.',MarkerEdgeColor='b', MarkerSize = 8);

hold on
plot(nodes(1, :), nodes(2, :), 'r.', MarkerSize = 10);

% hold on
% text(nodes(1,:),nodes(2,:), 'n' + string(1:size(nodes,2)));
% hold on
% text((nodes(1, edges(1,:)) + nodes(1, edges(2,:)))/2, ...
%      (nodes(2, edges(1,:)) + nodes(2, edges(2,:)))/2, ...
%     'e' + string(1:size(edges,2)));

axis equal off
%exportgraphics(gca,'SquareLargeHole-Geo.pdf','ContentType','vector');
savefig(gcf,'SquareLargeHole-Geo.fig','compact');
