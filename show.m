clear 
clc
addpath 'D:\matcode\helpers'
% source_dir = 'org/';
% source_dir = [source_dir,'/',name];
name = 'dragon_o.off';
% source_dir = 'org/';
% cover_mesh = [source_dir,'/',name];
%[vertex_o, face_o] = read_obj(name);
[vertex_o, face_o] = read_off(name);
vertex_o = vertex_o*10;
% vertex_o = vertex_o';
% vertex_o = vertex_o-repmat(mean(vertex_o),[size(vertex_o,1), 1]);
figure;
%figure('Position', [10, 10, 80, 60]); % 设置窗口大小
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]); 
set(gca, 'XDir', 'reverse');
%trisurf(F,x,y,z); 
plot_mesh(vertex_o, face_o); 
camzoom(0.5);
%shading interp; colormap jet(256);

Mesh.v=vertex_o';
Mesh.f=face_o';

name = 'dragon_w.off';
% source_dir = 'enc/';
% enc_mesh = [source_dir,'/',name];

%[vertex_e, face_e] = read_obj(enc_mesh);
[vertex_e, face_e] = read_off(name);
vertex_e = vertex_e*10;
figure;
%figure('Position', [10, 10, 80, 60]); % 设置窗口大小
axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]); 
set(gca, 'XDir', 'reverse');
%trisurf(F,x,y,z); 
plot_mesh(vertex_e, face_o); 
camzoom(0.5);
%vertex_e = vertex_e';
%vertex_e = vertex_e-repmat(mean(vertex_e),[size(vertex_e,1), 1]);

% 计算差的平方
differences = (vertex_o - vertex_e).^2;

% 对每列求和
sumOfSquares = sum(differences, 1);

% 对每个顶点计算欧几里得距离
euclideanDistances = sqrt(sumOfSquares)';
euclideanDistances = euclideanDistances * 120;


options.face_vertex_color = perform_saturation(euclideanDistances,1.2);


figure;
%camorbit(180, 0, 'z');
renderMesh(Mesh, euclideanDistances,15,180);%不显示边
colorbar