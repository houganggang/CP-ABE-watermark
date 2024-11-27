
%%
%[V,F] = readOBJ('VASE.obj');
clear;
addpath(genpath('D:/matcode/toolbox_graph-master'));
%[V, F] = read_obj('test.obj');  %3*n
[V, F] = read_off('bunny.off');
% vertex = V';
% [max_x, idax] = max(vertex(:,1)); [min_x, idix] = min(vertex(:,1)); %求出坐标的最大和最小
% [max_y, iday] = max(vertex(:,2)); [min_y, idiy] = min(vertex(:,2));
% [max_z, idaz] = max(vertex(:,3)); [min_z, idiz] = min(vertex(:,3));
%  xiao = [min_x min_y min_z]; x_in = max_x -min_x; y_in = max_y - min_y; z_in = max_z - min_z;
%  da = [max_x max_y max_z]; dda = max(x_in,y_in,z_in);
%  x_in = 1/x_in; y_in = 1/y_in; z_in = 1/z_in;
%  %matr = diag([x_in y_in z_in]);  
% matr = diag([dda dda dda]);
%  inv_matr = matr^-1;
% ver_01 = (vertex - repmat(xiao,[size(vertex,1), 1]))*matr;%转化为0到1之间的比例
% V = ver_01;
% [V1, F1] = read_obj('tes11.obj');  %3*n
% [hd] = HausdorffDist(V',V1');

%[vertex_new] = preprocess( V );
%V = vertex_new;
V = V'; F = F';  %n*3
%%
 vertex = V-repmat(mean(V),[size(V,1), 1]);
 [azimuth,elevation,ra] = cart2sph(vertex(:,1),vertex(:,2),vertex(:,3));
 raa = max(ra);
 vertex = vertex./raa;%将模型缩放到单位球 a^2+b^2+c^2=1
 write_mesh('bunny_o.off', vertex, F);
% 
% write_mesh('Arm2.off', vertex, F);
figure;
% trisurf(F,V(:,1),V(:,2),V(:,3)); 
% colormap(gray); axis equal;
% title('orignal');
plot_mesh(vertex',F'); shading interp; colormap jet(256);
title('original');
%vertex = V-repmat(mean(V),[size(V,1), 1]);
ryi = mean(vertex);
message_length = 32 ;  %水印长度
message = round(rand(message_length,1));  %水印
%encode = round(rand(message_length,1));   %加密序列
%encode = -0.5 + rand(1, message_length);
lengt  = size(vertex,1);
[azimuth,elevation,r] = cart2sph(vertex(:,1),vertex(:,2),vertex(:,3));%坐标转换
maxr = max(r); minr = min(r);
stand = (maxr - minr)/message_length; %每个水印的间隔
index_w = cell(message_length,1);  %水印序列存储
index_c = cell(message_length,1);   %归一化之后的顶点值 
qu_st = minr;  %开始的位置
%求出每个点所属的水印序列号
for i = 1:message_length
    ind_num = [];
    for j = 1:size(r,1)
        if (r(j)>=qu_st) &&(r(j)<qu_st+stand)
            ind_num = [ind_num j];
        end
    end
    index_w{i,1} = ind_num;
    qu_st = qu_st + stand;
end
for i = 1:message_length  %归一化
    ind_sub = index_w{i,1};ff = [];
    pp = min(r(ind_sub)); qq = max(r(ind_sub));
    for j = 1:size(ind_sub,2)
        aa = (r(ind_sub(j))-pp)/(qq-pp);
        ff(j) = aa;
    end
    index_c{i,1} = ff;
end
%%
%水印嵌入
kn = 0.001;%水印变化强度
a = 0.05;%水印强度参数
num = 0;qu_v = [];index_mean = 0;index_meann = [];
for i = 1:message_length
    if message(i)==0
        index_mean = mean(index_c{i,1});
        k = 1;meann = 0.5;k = k + kn;
        while(index_mean > meann - a)   %嵌入0 ，但是平均值大于 1/2，若小于1/2，则直接跳过
            for j=1:size(index_c{i,1},2)
                index_c{i,1}(j) = index_c{i,1}(j)^(k);
            end
            index_mean = mean(index_c{i,1});
            k = k + kn;
        end
    else   %嵌入l,但平均值小于1/2
        index_mean = mean(index_c{i,1});
        k = 1;meann = 0.5;k = k - kn;
        while(index_mean < meann + a)   %嵌入1 ，但是平均值小于 1/2，
            for j=1:size(index_c{i,1},2)
                index_c{i,1}(j) = index_c{i,1}(j)^(k);
            end
            index_mean = mean(index_c{i,1});
            k = k - kn;
        end
    end
    index_meann(i) = index_mean;
end

%%
%恢复模型
for i = 1:message_length
    ind_sub = [];
    for j = 1:size(index_c{i,1},2)
        ind_sub = index_w{i,1}; pp = min(r(ind_sub)); qq = max(r(ind_sub));
        %pp = min(index_c{i,1}); qq = max(index_c{i,1});
        %index_c{i,1}(j) = index_c{i,1}(j)*stand + (i-1)*stand + minr;%求出每个坐标的原始值
        index_c{i,1}(j) = index_c{i,1}(j)*(qq - pp) + pp;%求出每个坐标的原始值
    end
end
for i = 1: message_length
    for j = 1:size(index_c{i,1},2)
        r(index_w{i,1}(j)) = index_c{i,1}(j);%求出每个坐标的原始值
    end
end
%% wirte mesh
[x,y,z] = sph2cart(azimuth,elevation,r);
write_mesh('bunny_w.off', [x y z], F);
vec_rec = [x y z];
%hd值
[hd] = HausdorffDist(vertex,vec_rec);
disp('hd=');disp(hd);
figure;
%trisurf(F,x,y,z); 
plot_mesh(vec_rec', F'); shading interp; colormap jet(256);
%colormap(gray); axis equal;
%title('watermark');
%% 攻击类型
%旋转
% theta = 45; % 旋转角度 (degree)
% rotationMatrix = [cosd(theta), -sind(theta), 0;
%                   sind(theta), cosd(theta),  0;
%                   0,           0,           1];
% rotatedVertices =  vec_rec*rotationMatrix;
% vec_rec = rotatedVertices;
% figure;
% %trisurf(F,x,y,z); 
% plot_mesh(vec_rec', F'); shading interp; colormap jet(256);
% title('rotate');
%缩放
% scaleMatrix = [1.5, 0,   0;
%                0,   1.5, 0;
%                0,   0,   1.5];
% scaledVertices =  vec_rec * scaleMatrix;
% vec_rec = scaledVertices;
% figure;
% %trisurf(F,x,y,z); 
% plot_mesh(vec_rec', F'); shading interp; colormap jet(256);
% title('scale');
% %平移
% 5. 平移变换 (向右移动2，向上移动3，向前移动1)
% translationVector = [2; 3; 1];
% translationVector = repmat(translationVector', size(vec_rec,1), 1);
% translatedVertices = translationVector + vec_rec;
% vec_rec = translatedVertices;
%%
%高斯噪声，简化，平滑
% sigma = 0.01*0.01;  %标准差，方差为标准差的平方，均值为0 ，方差为多少？的噪声
% Y = awgn(vec_rec,45,'measured');%;Y2 = awgn(y,50);Y3 = awgn(z,50);
% vec_rec = Y;
% %Y = [Y1 Y2 Y3];
% ver_snr1 = meshSNR(vec_rec,Y);%直接利用原始坐标的精度更高。
%noise = normrnd(0,sigma,size(vec_rec));
%Y=imnoise(vec_rec,'gaussian',0.01,0);
%平滑
% fv1.vertices = vec_rec; fv1.faces = F;
% fv_smooth = smoothpatch(fv1, 0, 20, 0.2);
% vertex_s = fv_smooth.vertices;    %n*3
% face_s = fv_smooth.faces;
% vec_rec = vertex_s;
% %简化
vertex = vec_rec; face = F;
fe= size(vertex,1)*(0.9);
vertexflag = floor(fe);
%vertexflag = 2000;
[SimpV, SimpF] = simplification(vertex, face, vertexflag);
SimpF = SimpF+1;
vec_rec = SimpV; F = SimpF;
% figure;
% plot_mesh(SimpV', SimpF'); shading interp; axis tight;
% %  %trisurf(F, SimpV(:,1), SimpV(:,2), SimpV(:,3));
% % trimesh(face', vertex(:, 1), vertex(:, 2), vertex(:, 3));
% % colormap(gray); axis equal;
% title('simplifaction');
% %均匀量化
% depth = 7;
% [vertex_uq] = uniform_quantization(vertex, depth);

%%
%提取

%write_obj('tes11.obj', vec_rec, F);
%vec_rec = SimpV; F = SimpF; 
%vec_rec = Y;
%vec_rec = vertex_uq;
%vec_rec = vertex_s;
% figure;
% %trisurf(F,x,y,z); 
% plot_mesh(vec_rec', F'); shading interp; colormap jet(256);
cfr = mean(vec_rec);
gg = repmat(mean(vec_rec),[size(V,1), 1]);
vec_rec = vec_rec-repmat(mean(vec_rec),[size(vec_rec,1), 1]);
%diff = vertex - vec_rec;
[azimuth1,elevation1,rr] = cart2sph(vec_rec(:,1),vec_rec(:,2),vec_rec(:,3));
%difr = rr-r;
maxrr = max(rr); minrr = min(rr);
qu_str = minrr; standr = (maxrr - minrr)/message_length;
index_wr = cell(message_length,1);
%求出每个点所属的水印序列号vertex
for i = 1:message_length
    ind_num = [];
    for j = 1:size(rr,1)
        if (rr(j)>=qu_str) &&(rr(j)<=qu_str+standr)
            ind_num = [ind_num j];
        end
    end
    index_wr{i,1} = ind_num;
    qu_str = qu_str + standr;
end
index_meana = [];
for i = 1:message_length
    index = i-1;diff = [];
    index_v = index_wr{i,1};%find(index_w == index);  %每个分区索引的点坐标
    min_r = min(rr(index_v)); max_r = max(rr(index_v));
    qu_v = (rr(index_v)-min_r)/(max_r - min_r);  %归一化距离
    mean_v = mean(qu_v);
    index_meana(i) = mean_v;
    if (mean_v>0.5)
        water(i) = 1;
    else
        water(i) = 0;
    end 
end

err_dist = logical(message) - logical(water');
err_len = length(find(err_dist(:)~=0));
data_err_percent = 1-err_len/length(water);
disp('data_err_percent=');disp(data_err_percent);
% corr = corrcoef(message, water');
% disp('corr=');disp(corr);
[hd] = HausdorffDist(vertex,vec_rec);
disp('hd=');disp(hd);






