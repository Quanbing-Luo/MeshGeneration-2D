%% Mesh generation

clear;
meshgeneration('geometry.dat','meshradii.dat', Hmax = 0.12);


%% Import geometry and mesh data
[nodes_geo, elements_geo, edges_geo] = geometryread('geometry.dat');
[nodes, elements, radii] = meshradiiread('meshradii.dat');

%% Figure mesh

figure (1); clf;
patch('Faces',elements','Vertices',nodes','FaceColor',[0.85, 0.85, 0.85], ...
    'EdgeColor','b');

% hold on
% plot(nodes_geo(1, :), nodes_geo(2, :), 'r.', MarkerSize = 10);

axis equal off
exportgraphics(gca,'CircularDisk-MeshHmax.pdf','ContentType','vector');
savefig(gcf,'CircularDisk-MeshHmax.fig','compact');


%% Figure Mesh Quality Analysis

xs=nodes(1,:);ys=nodes(2,:);
t1=elements(1,:);t2=elements(2,:);t3=elements(3,:);
la=hypot(xs(t2)-xs(t1), ys(t2)-ys(t1));
lb=hypot(xs(t3)-xs(t1), ys(t3)-ys(t1));
lc=hypot(xs(t3)-xs(t2), ys(t3)-ys(t2));
rho=(lc+lb-la).*(lc+la-lb).*(la+lb-lc)./ (la.*lb.*lc);

figure (2); clf;
histogram(rho,0:0.1:1,'Normalization','probability');
set(gcf,'Position',[230 250 400 400]) ;
set(gca,'Position',[0.12 0.14 0.85 0.8]);
xlim([0 1]);
xlabel('Radius ratio','interpreter','latex');
ylabel('Percentage','interpreter','latex');
grid on
exportgraphics(gca,'CircularDisk-MeshQualityHmax.pdf','ContentType','vector');
savefig(gcf,'CircularDisk-MeshQualityHmax.fig','compact');


%% Figure Mesh with Bubbles

figure (3); clf;
patch('Faces',elements','Vertices',nodes','FaceColor',[0.85, 0.85, 0.85], ...
    'EdgeColor','b');

hold on  
for i=1:length(radii)
    r=radii(i);
    if r < 1e10
        rectangle('Position',[nodes(1,i)-r, nodes(2,i)-r, 2*r, 2*r], ...
    'Curvature', [1 1], 'EdgeColor','r');
    end
end
axis equal off
exportgraphics(gca,'CircularDisk-MeshRadiiHmax.pdf','ContentType','vector');
savefig(gcf,'CircularDisk-MeshRadiiHmax.fig','compact');