function [SimpV, SimpF] = simplification(vertices, faces, vertexflag)
% 简化封闭三角网格 / simplify closed tri-mesh
%
% 输入参数(vargin):
%          vertices: nv x 3, nv:顶点数
%             faces: nf x 3, nf:面数
%        vertexflag: ng    , ng:期望顶点数
%
% 输出参数(vargout)： 
%           SimpV：ng x 3 ,简化后的点
%           SimpF：?? x 3 ,简化后的面
%
% example：[SimpV, SimpF] = MySimp4(vertices, faces, 200)
%

% 统计模型的点数和面数
num_vertex = size(vertices, 1);

% 计算面法向量和面常量d
%faces = faces + 1; % 原始面索引从0开始，因此加1
normal_face = zeros(size(faces,1), 3);
d = zeros(size(faces,1), 1);
updata_normal_d(1:size(faces,1));

%% 计算每个点的Q
Q = zeros(num_vertex, 2); % 初始化
updataVertexQ(1:num_vertex);

%% 计算每个边的cost
temp = [faces(:,[1,2]); faces(:,[1,3]); faces(:,[2,3]);faces(:,[2,1]); faces(:,[3,1]); faces(:,[3,2])];
temp = triu(sparse(temp(:,1), temp(:,2), ones(size(temp,1), 1)));
[tempr, tempc] = find(temp~=0);
edges = [tempr, tempc];

%% 开始简化
num_newvertex = num_vertex;
while(num_newvertex > vertexflag)
    
    [~, token] = min(Q(:,1));
    ind1_coll_ver = Q(token, 2); % 保留点
    ind2_coll_ver = token; % 收缩点
    
    % 查找包含收缩点的面索引
    [neibor_ver1, ~] = neiborVertexFace(ind1_coll_ver);
    [neibor_ver2, neibor_face2] = neiborVertexFace(ind2_coll_ver);
    if( isempty(neibor_face2) )
        Q(token,:) = [intmax, nan];
        num_newvertex = sum(Q(:,1)~=intmax);
        continue;
    end
    neibor_edge2 = find(sum(edges==ind2_coll_ver, 2));
    common_vertex = intersect(neibor_ver1, neibor_ver2);
    
    if(~isempty(setdiff(neibor_ver2, [common_vertex, ind1_coll_ver])) || (numel(neibor_ver2) == numel(neibor_face2)))
        
        temp = faces(neibor_face2,:);
        temp(temp==ind2_coll_ver) = ind1_coll_ver;
        faces(neibor_face2,:) = temp;
        move_face = neibor_face2(find(sum(temp==ind1_coll_ver,2)==2));
        remain_face = setdiff(neibor_face2, move_face);
        updata_normal_d(remain_face);
        faces(move_face,:) = [];
        normal_face(move_face,:) = [];
        d(move_face) = [];
        updataVertexQ(neibor_ver2);
        temp = zeros(size(neibor_edge2, 1), 1);
        for i = [common_vertex, ind1_coll_ver]
            temp = temp + sum(edges(neibor_edge2,:)==i, 2);
        end
        move_edge = neibor_edge2(find(temp));
        remain_edge = setdiff(neibor_edge2, move_edge);
        temp = edges(remain_edge,:);
        temp(temp==ind2_coll_ver) = ind1_coll_ver;
        edges(remain_edge,:) = temp;
        
    else
        move_edge = neibor_edge2;
        faces(neibor_face2,:) = [];
        normal_face(neibor_face2,:) = [];
        d(neibor_face2) = [];
    end
    edges(move_edge,:) = [];
    Q(Q==token) = ind1_coll_ver;
    Q(token,:) = [intmax, nan];
    
    num_newvertex = sum(Q(:,1)~=intmax);
end
% 输出简化顶点和面
disp('Merging results ...')
verInd_remain = unique(faces(:));
SimpV = vertices(verInd_remain,:);
SimpF = faces;
for i = 1:length(verInd_remain)
    SimpF(SimpF == verInd_remain(i)) = i;
end
SimpF = SimpF - 1;
disp('Merging done.');

% ===============================================================================
    function updata_normal_d(index_face)
        vec1 = vertices(faces(index_face,1),:) - vertices(faces(index_face,2),:);
        vec2 = vertices(faces(index_face,1),:) - vertices(faces(index_face,3),:);
        raw_res = cross(vec1, vec2, 2);
        temp = sqrt(sum(raw_res.^2, 2));
        temp(temp < eps) = 1;
        normal_face(index_face,:) = raw_res ./ repmat(temp, 1, 3);
        d(index_face) = -sum(normal_face(index_face,:) .* vertices(faces(index_face, 1), :), 2);
    end

    function updataVertexQ(ind)
        if( size(ind,1) > size(ind,2) )
            ind = ind';
        end
        for iter = ind
            % 计算每个点的邻接点和邻接面
            [neibor_vertex, neibor_face] = neiborVertexFace(iter);
            if( isempty(neibor_face) )
                Q(iter,:) = [intmax, nan];
                continue;
            end
            % 计算每个点的邻接点的Q
            Q_neibor_vertex = zeros(length(neibor_vertex), 1);
            for j = 1:length(neibor_vertex)
                Q_neibor_vertex(j) = sum((sum(normal_face(neibor_face,:) .* repmat(vertices(neibor_vertex(j),:),length(neibor_face),1), 2) + d(neibor_face) ).^2);
            end
            % 每个点的Q=该点邻接点的Q的最大值
            [Q(iter,1), tind] = min(Q_neibor_vertex);
            Q(iter,2) = neibor_vertex(tind);
        end
    end

    function [ind_ver, ind_face] = neiborVertexFace(ind)
        ind_face = find(sum(faces==ind, 2));
        if( isempty(ind_face) )
            ind_ver = [];
            return;
        end
        ind_ver = unique(reshape(faces(ind_face, :), 1, 3*length(ind_face)));
        ind_ver(ind_ver==ind) = [];
    end
% ===============================================================================
end
