function [vertex_uq] = uniform_quantization(vertex, depth)

% 找出波形最大最小值   N*3
%scale the model in to unit cube 
[range_up]=max(vertex);   %N*3
[range_low]=min(vertex);
range=range_up-range_low;
[range_max,axis]=max(range);
for j = 1:3
    Max = range_up(j);
    Min = range_low(j);
    delv = (Max - Min)/2^depth;
    for i = 1 :2^depth+1
        m(i,j) = Min + delv*(i - 1);
    end
end

for j= 1:3
    for i =1:size(vertex,1)
        for k = 1:2^depth
            if vertex(i, j) >= m(k, j) && vertex(i, j) <= m(k+1,j)
                vertex_uq(i, j) = (m(k,j) + m(k+1,j))/2;
                break
            end
        end
    end
end

% 画出原声波以及量化之后的声波
% figure;
% plot(wav,'r'); hold on
% plot(wav_tread,'g'):title('下取整');hold off
end
