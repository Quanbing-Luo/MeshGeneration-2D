%% Mesh generation
clear;

% Create a PDE model.
model = createpde;

%Geometry from Triangulated Mesh
importGeometry(model,"SquareLargeHole.stl");

figure (1); clf;
pdegplot(model,VertexLabels = "on",EdgeLabels = "on");
axis equal off

mesh = generateMesh(model, GeometricOrder = "linear", ...
    Hedge = {3,0.01}, Hgrad=1.1);
nodes = mesh.Nodes; elements =mesh.Elements;
save("mesh.mat","nodes","elements");

%% Import geometry and mesh data

TR = stlread('SquareLargeHole.stl');
nodes_geo = TR.Points(:,1:2)';
elements_geo=TR.ConnectivityList'; 
load("mesh.mat","nodes","elements");


%% Figure mesh

figure (2); clf;
patch('Faces',elements','Vertices',nodes','FaceColor',[0.85, 0.85, 0.85], ...
    'EdgeColor','b');

hold on
plot(nodes_geo(1, :), nodes_geo(2, :), 'r.', MarkerSize = 10);

axis equal off
% exportgraphics(gca,'SquareLargeHole-MeshToolboxHedge.pdf','ContentType','vector');
% savefig(gcf,'SquareLargeHole-MeshToolboxHedge.fig','compact');


%% Figure Mesh Quality Analysis

xs=nodes(1,:);ys=nodes(2,:);
t1=elements(1,:);t2=elements(2,:);t3=elements(3,:);
la=hypot(xs(t2)-xs(t1), ys(t2)-ys(t1));
lb=hypot(xs(t3)-xs(t1), ys(t3)-ys(t1));
lc=hypot(xs(t3)-xs(t2), ys(t3)-ys(t2));
rho=(lc+lb-la).*(lc+la-lb).*(la+lb-lc)./ (la.*lb.*lc);

figure (3); clf;
histogram(rho,0:0.1:1,'Normalization','probability');
set(gcf,'Position',[230 250 400 400]) ;
set(gca,'Position',[0.12 0.14 0.85 0.8]);
xlim([0 1]);
xlabel('Radius ratio','interpreter','latex');
ylabel('Percentage','interpreter','latex');
grid on
% exportgraphics(gca,'SquareLargeHole-MeshToolboxQualityHedge.pdf','ContentType','vector');
% savefig(gcf,'SquareLargeHole-MeshToolboxQualityHedge.fig','compact');