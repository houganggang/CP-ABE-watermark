%[V,F] = readOBJ('VASE.obj');
clear;
addpath(genpath('D:/matcode/toolbox_graph-master'));
[V, F] = read_off('dania.off');  %3*n 读取off模型
vertex = V';
[max_x, idax] = max(vertex(:,1)); [min_x, idix] = min(vertex(:,1)); %求出坐标的最大和最小
[max_y, iday] = max(vertex(:,2)); [min_y, idiy] = min(vertex(:,2));
[max_z, idaz] = max(vertex(:,3)); [min_z, idiz] = min(vertex(:,3));
 xiao = [min_x min_y min_z]; x_in = max_x -min_x; y_in = max_y - min_y; z_in = max_z - min_z;
 da = [max_x max_y max_z]; dda = max(x_in,y_in,z_in);
 x_in = 1/x_in; y_in = 1/y_in; z_in = 1/z_in;
 matr = diag([x_in y_in z_in]);  
 inv_matr = matr^-1;
ver_01 = (vertex - repmat(xiao,[size(vertex,1), 1]))*matr;%转化为0到1之间的比例
%[vertex_new] = preprocess( V );
%V = vertex_new;
V = V'; F = F';  %n*3
vertex = V-repmat(mean(V),[size(V,1), 1]); %归一化到以0为中心的圆上
%k = size(V,1); %特征值的个数
k = 84;ccc = true; V1=V; %k为需要的前k位特征值
a =1;   %水印强度
message_length = 16 ;  %水印长度，
message = round(rand(message_length,1));  %水印
%encode = round(rand(message_length,1));   %加密序列
encode = -0.5 + rand(1, message_length);
lengt  = size(vertex,1);
% [azimuth,elevation,r] = cart2sph(vertex(:,1),vertex(:,2),vertex(:,3));
% maxr = max(r); minr = min(r);
% stand = (maxr - minr)/message_length;
index_w = [];
%求出每个点所属的水印序列号
kn = 0.001;%水印变化强度
a = 0.03;num = 0;qu_v = [];index_mean = [];
%%
%%
%加密序列,-0.5到0。5的均匀分布
itear = 0;
while(ccc) %ccc,循环停止标记
    
%[evecs evals error] = laplace_beltrami_spectrum('bunny.off', k);
L = cotmatrix(V,F); % vertex*vertex Q some implementation may put L = 0.5*cotmatrix(V,F); this just scales ED by 0.5 and does not really matter.
%L = -L;
M = massmatrix(V,F,'barycentric');  %B,D vertex*vertex
%[EV,ED] = eigs(L,M,k,'sm','IsCholesky');  %eV是特征向量HK，largestabs  smallestabs   EV = V*k
[EV,ED] = eigs(L,M,84,'sm');%看函数说明
Minv = sqrt(diag(1./diag(M)));    % get M^{-1/2} for symmetry / faster than inv()
% Minv1 = sqrt(inv(M));
%inv_diff = Minv - Minv1;
Laplace_Beltrami = Minv * L * Minv; % get positive semi-definite discrete Laplacian
Laplace_Beltrami = Laplace_Beltrami * -1; % for positive eigenvalues
Laplace_Beltrami = (Laplace_Beltrami + Laplace_Beltrami.') * 0.5 ; % handle numerical precision issue: http://stackoverflow.com/a/33259074 

[~, eigen_val, eigen_vect] = svds(Laplace_Beltrami,84,'smallest');
%[evecs evals error] = laplace_beltrami_spectrum('dania.off', k)
%[evecs, evals, W, A] = main_mshlp('cotangent', shape, num_evecs);
%%
%特征向量，以矩阵形式返回。V 中的各列对应于沿 D 对角线方向的特征值。V 的形式和归一化取决于输入参数的组合：
%[V,D] = eigs(A) 返回矩阵 V，其各列是 A 的右特征向量，这样 A*V = V*D。V 中的特征向量已归一化，因此每个向量的 2-范数为 1。
%如果 A 为对称矩阵，则特征向量 V 为正交矩阵。
%[V,D] = eigs(A,B) 返回矩阵 V，其中各列为满足 A*V = B*V*D 的广义右特征向量。每个特征向量的 2-范数不一定为 1。
%如果 B 是对称正定矩阵，则 V 中的特征向量是归一化向量，因此每个特征向量的 B-范数均为 1。如果 A 也是对称矩阵，则特征向量与 B 正交。
%%
ED = - ED; % minus sign is needed since cotmatrix() gives a negative definite matrix.
sym = [];
for i=1:k  %看mainfod的求解过程
    Hk(:,i) = EV(:,i)/sqrt(EV(:,i)'*M*EV(:,i));  
    %Hk(:,i) = (Hk(:,i) - min(Hk(:,i)))/(max(Hk(:,i)) - min(Hk(:,i)));
    %Hk(:,i) = EV(:,i)/sqrt(sum(EV(:,i).^2));
    %Hk(:,i) = Hk(:,i)/sqrt(sum(Hk(:,i).^2));
    sym(i) = Hk(:,i)'*Hk(:,i);
    
end

x = V(:,1)'; y = V(:,2)'; z = V(:,3)';
xx = x*M*Hk; yy = y*M*Hk; zz = z*M*Hk;    %求出频谱系数,嵌入水印后的频谱系数会不同，尽量减少系数的不同
x_w = [];cc = [];y_w = [];z_w = [];
for i=1:k
    %x_w(i) = xx(i)+0.001; y_w(i) = yy(i)+0.001; z_w(i) = zz(i)+0.001;
    cc(i) = sqrt(xx(i).^2 + yy(i).^2 + zz(i).^2);      %频谱系数的修改原语，这就是嵌入水印的载体
end

vec_rec = [];
x_sp = Hk*xx'; 
y_sp = Hk*yy';
z_sp = Hk*zz';
x_d = x'-x_sp;
y_d = y'-y_sp;
z_d = z'-z_sp;

%% todo:计算步长和嵌入强度 ，通过论文的同步编码嵌入水印
beta = 0.0001;%beta = 0.00001;
S = beta*cc(2);    %量化步长= 
a =1;   %水印强度
message_length = 16;             %水印长度
message = round(rand(message_length,1));  %水印
%encode = round(rand(message_length,1));   %加密序列
encode = -0.5 + rand(1, message_length);   %加密序列,-0.5到0。5的均匀分布
Sn = [];yn = [];water = [];diff_ck = [];
num =21;j=1;
while(num <= k)   %84-20
    if mod(num,4)~=0
        i = mod(j,16);
        if i==0
           i = 16;
        end
        an = S*(message(i)/2);%+encode(i));   %抖动  -0.5  --- 1  错了，均匀量化不是向下取整
        xn = cc(num);
        water(j) = xn;
        Q = (xn)/S;
        if ((Q - floor(Q))>=0.5)
           if message(i)==1
               uk = (floor(Q) + 0.5 )*S;
           else
               uk = (floor(Q) + 1 )*S;
           end
        else
           if message(i)==1
               uk = (floor(Q) + 0.5 )*S;
           else
               uk = (floor(Q))*S;
           end
        end
%     err = floor((xn-an)/S)*S;
         %uk = Q*S + an - xn;
%     %num = floor(num);
         Sn(j) = xn + a*(uk - xn);
         j=j+1;
         %xn = 106.7*s , an  = (0.5+0.1)*s,量化后为106.1取整为106,106+0.6-106.7 =-0.1*@ = 106.6 
%     %disp('Sn='),disp(Sn(j));
%      diff_ck(j) = Sn(j) - xn;
        %Q = (xn-an)/S;
        %disp('xn='),disp(xn);
        %z = xn/S;
        %uck1 = ceil(z)*S+an;
        %uck2 = floor(z)*S+an;
%         if (abs(xn-uck1) > (abs(xn-uck2)))
%             uck = uck2;
%         else
%             uck = uck1;
%         end
        %Sn(j) = xn + a*(uck - xn);
       % j=j+1;
    end
    num=num+1;
end
% for j = 1:64
%     an = S*(message(j)/2+encode(j));   %抖动  -0.5  --- 1  错了，均匀量化不是向下取整
%     xn = cc(j+20);
%     water(j) = xn;
%     disp('xn='),disp(xn);
%     z = xn/S;
%     uck1 = ceil(z)*S+an;
%     uck2 = floor(z)*S+an;
%     if (abs(xn-uck1) > (abs(xn-uck2)))
%         uck = uck2;
%     else
%         uck = uck1;
%     end
%     Sn(j) = xn + a*(uck - xn);
%     num = (xn-an)/S;
%     if num - floor(num)>=0.5
%         num = floor(num) + 1;
%     else
%         num = floor(num);
%     end
%     err = floor((xn-an)/S)*S;
%     uk = num*S + an - xn;
%     %num = floor(num);
%     Sn(j) = xn + a*(num*S + an - xn);   %xn = 106.7*s , an  = (0.5+0.1)*s,量化后为106.1取整为106,106+0.6-106.7 =-0.1*@ = 106.6 
%     %disp('Sn='),disp(Sn(j));
%      diff_ck(j) = Sn(j) - xn;
% end
pre = [];
for j = 1:48
        i = mod(j,16);
        if i==0
           i = 16;
        end
        num = (Sn(j))/S;%-encode(i)*S
       % num = floor(num);
        if num - floor(num)>=0.5
            num = floor(num) + 1;
        else
            num = floor(num);
        end
        yn(j) = num*S - (Sn(j));%水印提取  （106.6-0.1） 106.5取整为106，减去106.5 所以106.5-106 =0.5     
        pre(j) = yn(j)/S;
end
for i = 1:message_length*3
    if abs(yn(i))>=(S/4)/a && abs(yn(i)) <=(S*0.75)/a
        yn(i) = 1;
    else
        yn(i) = 0;
    end
end
err2 = 0;
for i = 1:message_length
    if (yn(i)+yn(i+16)+yn(i+32)>=2)
        zn(i)=1;
        if (yn(i)+yn(i+16)+yn(i+32)==2)
            err2 = err2 +1;
        end
    else
        zn(i)=0;
        if (yn(i)+yn(i+16)+yn(i+32)==1)
           err2 = err2 +1;
        end 
    end
end
%%
%嵌入完水印计算参数
%  计算提取误码率
err_dist = logical(message) - logical(zn');
err_len = length(find(err_dist(:)~=0));
data_err_percent = err_len/length(zn);
disp('data_err_percent='),disp(data_err_percent);
disp('err2='),disp(err2);
j = 1;
wwww = [];
%计算嵌入水印之后的频谱系数，前20个不嵌入水印，所以是84-20
for i = 1:64
    if mod(i,4)~=0
        wwww(j) = Sn(j)/water(j);
        xx(i+20) =  Sn(j)/water(j)*xx(i+20);
        yy(i+20) =  Sn(j)/water(j)*yy(i+20);
        zz(i+20) =  Sn(j)/water(j)*zz(i+20);
        j = j+1;
    end
end

%%
%xw_rec = x_w*Hk'; yw_rec = y_w*Hk'; zw_rec = z_w*Hk';
%恢复嵌入水印的坐标值，
x_rec = xx*Hk'; y_rec = yy*Hk'; z_rec = zz*Hk';   %类似于稀疏表示
x_rec = x_rec'+ x_d; y_rec = y_rec'+ y_d; z_rec = z_rec'+ z_d;
vec_rec = [x_rec y_rec z_rec]; %带有水印的坐标
%%
%重构频谱系数 接下来是计算从带水印的坐标中提取水印，因为有因果效应，所以嵌入的水印提取会有不正确的情况
%需要迭代控制，得到最后的水印提取满足要求
%k = 300; %特征值的个数
L1 = cotmatrix(vec_rec,F); % vertex*vertex Q some implementation may put L = 0.5*cotmatrix(V,F); this just scales ED by 0.5 and does not really matter.
M1 = massmatrix(vec_rec,F,'barycentric');  %B,D vertex*vertex
[EV1,ED1] = eigs(L1,M1,k,'sm');  %eV是特征向量HK，largestabs  smallestabs   EV = V*k
%diff_ED = ED1 - ED;

ED1 = - ED1; % minus sign is needed since cotmatrix() gives a negative definite matrix.
diff_ED = ED1 - ED;
sym1 = [];
for i=1:k
    Hk1(:,i) = EV1(:,i)/sqrt(EV1(:,i)'*M1*EV1(:,i));
    %Hk1(:,i) = EV1(:,i)/sqrt(sum(EV1(:,i).^2));
   % Hk1(:,i) = Hk1(:,i)/sqrt(sum(Hk1(:,i).^2));
    sym1(i) = Hk(:,i)'*Hk1(:,i);
    %Hk1(:,i) = EV1(:,i);
end
x = vec_rec(:,1)'; y = vec_rec(:,2)'; z = vec_rec(:,3)';%重构出来的顶点坐标
%求出重构顶点和原始顶点的误差
d_x = (V(:,1) - x'); d_y = (V(:,2) - y'); d_z = (V(:,3) - z');

xx1 = x*M1*Hk1; yy1 = y*M1*Hk1; zz1 = z*M1*Hk1;  %求出频谱系数
diff_x = [];diff_y = [];diff_z = [];diff_c = [];
for i=1:84
    %x_w(i) = xx(i)+0.001; y_w(i) = yy(i)+0.001; z_w(i) = zz(i)+0.001;
    diff_x(i) = xx1(i) - xx(i); diff_y(i) = yy1(i) - yy(i); diff_z(i) = zz1(i) - zz(i); 
    cc1(i) = sqrt(xx1(i).^2 + yy1(i).^2 + zz1(i).^2);      %频谱系数的修改原语
    diff_c(i) = (cc1(i) - cc(i))/S;
end
S1 = beta*cc1(2); 
Sn1 = [];   %量化步长
%%
kk=0;pre1 = [];
for j = 1:64
    if mod(j,4) ~=0
        kk=kk+1;
        Sn1(kk) = cc1(j+20);
        %kk=kk+1;
        i = mod(kk,16);
        if i==0
           i = 16;
        end
        num = (Sn1(kk))/S1;
         if num - floor(num)>=0.5
             num = floor(num) + 1;
         else
             num = floor(num);
         end 
        yn1(kk) = num*S1 - (Sn1(kk));   %水印提取  （106.6-0.1） 106.5取整为106，减去106.5 所以106.5-106 =0.5
        pre1(kk) = yn1(kk)/S1;
    end
end

for i = 1:message_length*3
    if abs(yn1(i))>=(S1/4) && abs(yn1(i)) <=(S1*0.75)
        yn1(i) = 1;
    else
        yn1(i) = 0;
    end
end
err21 = 0;
for i = 1:message_length
    if (yn1(i)+yn1(i+16)+yn1(i+32)>=2)
        zn1(i)=1;
        if (yn1(i)+yn1(i+16)+yn1(i+32)==2)
            err21 = err21 +1;
        end
    else
        zn1(i)=0;
        if (yn1(i)+yn1(i+16)+yn1(i+32)==1)
           err21 = err21 +1;
        end 
    end
end
%%
% for j = 1:64
%     if mod(j,4) ~=0
%         Sn1(j) = cc1(j+20);
%         i = mod(j,16);
%         if i==0
%            i = 16;
%         end
%         num = (Sn1(j)-encode(j)*S1)/S1;
%         if num - floor(num)>=0.5
%             num = floor(num) + 1;
%         else
%             num = floor(num);
%         end 
%         yn1(j) = num*S1 - (Sn1(j) - encode(j)*S1);   %水印提取  （106.6-0.1） 106.5取整为106，减去106.5 所以106.5-106 =0.5
%     end
% end
% %message= 64;
% yn1 = yn1/S1;
% for i = 1:message_length
%     if abs(yn1(i))>(0.25) 
%         yn1(i) = 1;
%     else
%         yn1(i) = 0;
%     end
% end
%%
%  计算提取误码率
err_dist = logical(message) - logical(zn1');
err_len = length(find(err_dist(:)~=0));
data_err_percent = err_len/length(zn1);
disp('data_err_percent='),disp(data_err_percent);
disp('err21='),disp(err21);
if (err_len ==0 && err21 <6)
    ccc = false;
end
itear=itear+1;
disp('itear='),disp(itear);
V = vec_rec;
end
%%
snr = meshSNR(vec_rec,V1);
disp('SNR='),disp(snr);

[hd] = HausdorffDist(V1,vec_rec);
disp('Hausd='),disp(hd);
qq = sum((V1-vec_rec).^2,2);
MRSE1 = sqrt(sum(sum((V1-vec_rec).^2,2))/size(V,1));
MRSE2 = sqrt(sum(sum((vec_rec-V1).^2,2))/size(V,1));
disp('MRSE1='),disp(MRSE1);disp('MRSE2='),disp(MRSE2)
%%

%h=plot_mesh(vertex,face_o');
clf;
subplot(1,2,1);
trisurf(F, V(:,1), V(:,2), V(:,3));
%trimesh(face_o', vertex(:, 1), vertex(:, 2), vertex(:, 3));
colormap(gray); axis equal;
title('original');
subplot(1,2,2);
trisurf(F, x_rec, y_rec, z_rec);
%trimesh(face_o', vertex(:, 1), vertex(:, 2), vertex(:, 3));
 colormap(gray); axis equal;
title('recover11111');

% subplot(2,2,3);
% trisurf(F, xw_rec', yw_rec', zw_rec');
% %trimesh(face_o', vertex(:, 1), vertex(:, 2), vertex(:, 3));
%  colormap(gray); axis equal;
% title('watermark recover');

%%
%计算MSDM 
%%
%  求出模型的体积
v_num = size(F,1);m000_q = [];          %n*3   
m100_q = [];m010_q = [];m001_q = [];m200_q = [];m020_q = [];m002_q = [];m110_q = [];m101_q = [];m011_q = [];
for i = 1:v_num
    x1 = V(F(i,1),1);x2 = V(F(i,2),1);x3 = V(F(i,3),1);  %n*3
    y1 = V(F(i,1),2);y2 = V(F(i,2),2);y3 = V(F(i,3),2);
    z1 = V(F(i,1),3);z2 = V(F(i,2),3);z3 = V(F(i,3),3);
    m000 = (1/6)*abs(x1*y2*z3 - x1*y3*z2 - y1*x2*z3 + y1*x3*z2 + z1*x2*y3 - z1*x3*y2);
    m100 = (m000/4)*(x1 + x2 + x3);
    m010 = (m000/4)*(y1 + y2 + y3);
    m001 = (m000/4)*(z1 + z2 + z3);
    m200 = (m000/10)*(x1^2 + x2^2 + x3^2 +x1*x2 + x1*x3 + x2*x3);
    m020 = (m000/10)*(y1^2 + y2^2 + y3^2 +y1*y2 + y1*y3 + y2*y3);
    m002 = (m000/10)*(z1^2 + z2^2 + z3^2 +z1*z2 + z1*z3 + z2*z3);
    m110 = (m000/10)*(x1*y1 + x2*y2 + x3*y3 + (x1*y2 + x1*y3 + x2*y1 + x2*y3 + x3*y1 + x3*y2)/2);
    m101 = (m000/10)*(x1*z1 + x2*z2 + x3*z3 + (x1*z2 + x1*z3 + x2*z1 + x2*z3 + x3*z1 + x3*z2)/2);
    m011 = (m000/10)*(z1*y1 + z2*y2 + z3*y3 + (z1*y2 + z1*y3 + z2*y1 + z2*y3 + z3*y1 + z3*y2)/2);
    
    m000_q(i) = m000;m100_q(i) = m100;m010_q(i) = m010;m001_q(i) = m001;m200_q(i) = m200;m020_q(i) = m020;
    m002_q(i) = m002;m110_q(i) = m110;m101_q(i) = m101;m011_q(i) = m011;
end
%全局体积矩计算质心
c_x = sum(m100_q)/sum(m000_q); c_y = sum(m010_q)/sum(m000_q); c_z = sum(m001_q)/sum(m000_q);
C = [c_x c_y c_z];
M = [sum(m200_q) sum(m110_q) sum(m101_q);
     sum(m110_q) sum(m020_q) sum(m011_q);
     sum(m101_q) sum(m011_q) sum(m002_q)];
 
 C_dis = mean(V);%离散质心
 V2C = V - repmat(C_dis,size(V,1),1);
%%
 %将模型缩放到单位球
 %噪声攻击？
 %均匀量化攻击
 %单位球化
 %平滑
 %cc = awgn();dd = wgn();
% 1、R = normrnd(MU,SIGMA)   2、R = normrnd(MU,SIGMA,m)   
% 3、R = normrnd(MU,SIGMA,m,n) 4、假设输入信号为X,则给X加上一个均值为0,方差为1的高斯白噪声信号的方法为：
sigma = 0.001*0.001;  %标准差，方差为标准差的平方，均值为0 ，方差为多少？的噪声
noise = normrnd(0,sigma,size(V));
Y=V+normrnd(0,sigma,size(V));   %构造出噪声
%平滑
fv1.vertices = V; fv1.faces = F;
fv_smooth = smoothpatch(fv1, 1, 5);
vertex_s = fv_smooth.vertices;    %n*3
face_s = fv_smooth.faces;
%简化
vertex = V; face = F;
vertexflag = 2000;
[SimpV, SimpF] = simplification(vertex, face, vertexflag);
figure;
plot_mesh(SimpV', SimpF'); shading interp; axis tight;
title('simplifaction');
%均匀量化
depth = 7;
[vertex_uq] = uniform_quantization(vertex, depth);
%旋转，缩放，平移
aa = 2; %硕放因子
rm.vertices = V;
rm.faces = F;axiss = [0 0 1]; angle = pi/4;
N = mesh_rotate(rm, axiss, angle);

%%
snr = meshSNR(vec_rec,V);
disp('SNR='),disp(snr);

[hd] = HausdorffDist(V,vec_rec);
disp('Hausd='),disp(hd);

%%
%mrms  max rms
qq = sum((V-vec_rec).^2,2);
MRSE1 = sqrt(sum(sum((V-vec_rec).^2,2))/size(V,1));
MRSE2 = sqrt(sum(sum((vec_rec-V).^2,2))/size(V,1));
disp('MRSE1='),disp(MRSE1);disp('MRSE2='),disp(MRSE2);
%%
%xk = sum(dd);
% [V,F] = read_vtk('Callosum_Forceps_Major_surf_decpt5.vtk');
% k = 800;
% L = cotmatrix(V',F'); % some implementation may put L = 0.5*cotmatrix(V,F); this just scales ED by 0.5 and does not really matter.
% M = massmatrix(V',F','barycentric');
% [EV,ED] = eigs(L,M,k,'sm');
% ED = - ED; % minus sign is needed since cotmatrix() gives a negative definite matrix.

