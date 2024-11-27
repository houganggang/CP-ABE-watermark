function [] = mainfold_emb()

L = cotmatrix(V,F); % vertex*vertex Q some implementation may put L = 0.5*cotmatrix(V,F); this just scales ED by 0.5 and does not really matter.
%L = -L;
M = massmatrix(V,F,'barycentric');  %B,D vertex*vertex
%[EV,ED] = eigs(L,M,k,'sm','IsCholesky');  %eV是特征向量HK，largestabs  smallestabs   EV = V*k
[EV,ED] = eigs(L,M,84,'sm');
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
for i=1:k
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
    cc(i) = sqrt(xx(i).^2 + yy(i).^2 + zz(i).^2);      %频谱系数的修改原语
end
vec_rec = [];
x_sp = Hk*xx'; 
y_sp = Hk*yy';
z_sp = Hk*zz';
x_d = x'-x_sp;
y_d = y'-y_sp;
z_d = z'-z_sp;
clf;
subplot(1,2,1);
trisurf(F, V(:,1), V(:,2), V(:,3));
%trimesh(face_o', vertex(:, 1), vertex(:, 2), vertex(:, 3));
colormap(gray); axis equal;
title('original');
subplot(1,2,2);
trisurf(F, x_sp, y_sp, z_sp);
%trimesh(face_o', vertex(:, 1), vertex(:, 2), vertex(:, 3));
 colormap(gray); axis equal;
title('recover');
vec_rec = [x_sp y_sp z_sp];
% snr = meshSNR(vec_rec,V1);
% disp('SNR='),disp(snr);

%% 利用量化来形式化误差



%% todo:计算步长和嵌入强度
beta = 0.0001;%beta = 0.00001;
S = beta*cc(2);    %量化步长= 
a =1;   %水印强度
message_length = 16;             %水印长度
message = round(rand(message_length,1));  %水印
%encode = round(rand(message_length,1));   %加密序列
encode = -0.5 + rand(1, message_length);   %加密序列,-0.5到0。5的均匀分布
Sn = [];yn = [];water = [];diff_ck = [];
num =21;j=1;
while(num <= k)   %84-20,前84个低频系数，前20个低频系数不嵌入，后64个低频系数嵌入，每四个一组。
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
%  计算提取误码率
err_dist = logical(message) - logical(zn');
err_len = length(find(err_dist(:)~=0));
data_err_percent = err_len/length(zn);
disp('data_err_percent='),disp(data_err_percent);
disp('err2='),disp(err2);
j = 1;
wwww = [];
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

x_rec = xx*Hk'; y_rec = yy*Hk'; z_rec = zz*Hk';   %类似于稀疏表示
x_rec = x_rec'+ x_d; y_rec = y_rec'+ y_d; z_rec = z_rec'+ z_d;
vec_rec = [x_rec y_rec z_rec];


end