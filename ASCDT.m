% ----------------------------------------------------------------------------
% ASCDT��
%   ����Delaunay�������Ŀռ�����㷨
% Author��
%   Qiliang Liu, Jianbo Tang, Huang Jincai, China.
% Reference:
%   ������.����Delaunay������������Ӧ�ռ�����㷨[D].���ϣ����ϴ�ѧ˶ʿ����,2011
% ----------------------------------------------------------------------------
% ****************************************************************************
% **********************        ASCDT()             **************************
% ****************************************************************************
function IDX = ASCDT(geodata, MVD, A, B)
% ������
% Geodata : spatial data: <[X,Y]>
% MVD     : visual distance, 0 by default��<float: 0>
% A       : parameter for global cut��<float: 1>
% B       : parameter for local cut��<float: 1.1>
% ���أ�
% IDX     : index for each point
% 
if nargin<2||isempty(MVD)
    MVD = 0; % ��С�ĳ��ȶ�����ʶ��
end
if nargin<3||isempty(A)
    A = 1;
end
if nargin<4||isempty(B)
    B = 1.1;
end

% ����������
Graph = triedge(geodata,'n'); % �ɿռ��ʵ�幹�ɵ�ͼ(��,Delanuay������),<matrix��[StartPointID,EndPointID,Length]>
% ��Graph����ͨ�Խ��в��ԣ��ֽ�Ϊ��ͼ
G = subgraph(Graph);
% �Ը���ͼ���з���
Cluster = [];
for i=1:size(G,1)
    C = ASCDT_CLUST(G{i}, MVD, A, B);
    Cluster = [Cluster; C];
    clear C;
end
IDX = ClusterID(geodata,Cluster); % ����IDX
if nargout<1
    scattx(geodata,IDX);
end
end % // ASCDT()




% ****************************************************************************
% **********************       ASCDT_CLUST()        **************************
% ****************************************************************************
function [Cluster,SubGraph,Active_Edge] = ASCDT_CLUST(Graph, MVD, A, B)
% ������
% Graph : �ɿռ��ʵ�幹�ɵ�ͼ(��,Delanuay������),<matrix��[StartPointID,EndPointID,Length]>
% MVD   : ��С���Ӿ��룬<float: 0>
% A     : ���峤���޼�ʱ�ĵ���ϵ����<float: 1>
% B     : �ֲ������޼�ʱ�ĵ���ϵ����<float: 1>
% ���أ�
% IDX   : ͼGraph�е�ÿ����'point'�����Ĵ����.
% 
if nargin<2||isempty(MVD)
    MVD = 0; % ��С�ĳ��ȶ�����ʶ��
end
if nargin<3||isempty(A)
    A = 1;
end
if nargin<4||isempty(B)
    B = 1.1;
end
% ���峤��
[SubGraph, Active_Edge] = global_cut(Graph,MVD,A);
% �ֲ�����
SubGraph_Nums = size(SubGraph,1);
Cluster = [];
for i=1:SubGraph_Nums
    Gi = SubGraph{i};
    % ����ÿһ��Gi���оֲ������޼��������ͼ
    Sub_Gi = local_cut(Gi,MVD,B);
    Cluster = [Cluster; Sub_Gi];
end
end % ACSDT_CLUST()



% ****************************************************************************
% **********************         triedge()          **************************
% ****************************************************************************
%%
% -------------------------------------------------------------------------
% triedge()
% ��ȡ������Delaunay�������ı�Edges�ͱߵĳ���Lengths
% ������
%     [Edges,Lengths]=triedge(geodata,showfig,color,newfigure)
% ���룺
%     geodata - <double npts*nattri>: �ռ����ݼ�(��[]);
%     showfig - <char 1*n>: �Ƿ���ʾDelaunay������('no');
%     color - <char 1*1>|<double 1*3>: ����������ʱ�߳���ɫ('b').
%     newfigure - �Ƿ����´����л�ͼ
% �����
%     Edges:   ��<[�����X�е��������յ���X�е�����]>
%     Lengths: �߳�(��Edges���Ӧ).
% -------------------------------------------------------------------------
function Output=triedge(geodata,showfig,color,newfigure)
%
if nargin<2||isempty(showfig)
    showfig='y';
end
if nargin<3||isempty(color)
    color='k';
end
if nargin<4||isempty(newfigure)
    newfigure='y';
end
[ptsnums,attnums]=size(geodata);
if ptsnums<1||attnums<2
    error('   The first arguments is not suit.');
end
TRI=delaunay(geodata(:,1),geodata(:,2));
trinums=size(TRI,1);
Edges=zeros(trinums*3,2);
for i=1:trinums
    Edges(3*i-2,:)=sort([TRI(i,1),TRI(i,2)]);
    Edges(3*i-1,:)=sort([TRI(i,1),TRI(i,3)]);
    Edges(3*i,  :)=sort([TRI(i,2),TRI(i,3)]);
end
Edges=unique(Edges,'rows');
enums=size(Edges,1);
Lengths=zeros(enums,1);
switch lower(showfig)
    case {'no','n','off',0}
        for i=1:enums
            Lengths(i)=sqrt((geodata(Edges(i,1),1)-geodata(Edges(i,2),1))^2+(geodata(Edges(i,1),2)-geodata(Edges(i,2),2))^2);
        end
    case {'yes','y','on',1}
        if isequal(lower(newfigure),'yes')||isequal(lower(newfigure),'y')||isequal(newfigure,1)
            figure;
        end
        for i=1:enums
            Lengths(i)=sqrt((geodata(Edges(i,1),1)-geodata(Edges(i,2),1))^2+(geodata(Edges(i,1),2)-geodata(Edges(i,2),2))^2);
            hold on;
            plot([geodata(Edges(i,1),1),geodata(Edges(i,2),1)],[geodata(Edges(i,1),2),geodata(Edges(i,2),2)],color);
        end
        box on;
    otherwise
end
Output = [Edges,Lengths];
end % // triedge()



% ****************************************************************************
% **********************         subgraph()         **************************
% ****************************************************************************
%%
% �������ͨ��ͼ
function Sub_Graph = subgraph(Graph)
% Graph -[StartPointID, EndPointID, Length]
PointID = unique([Graph(:,1);Graph(:,2)]);  % ����Graph�ĵ��ID�ţ���ԭʼ�ռ����ݼ��е������ţ�
PointNums = length(PointID);
% �������ͨ��ͼ
IDX = zeros(PointNums,1);
IsVisited = zeros(PointNums,1);
Index = (1:PointNums);
ClusterID = 0;

for k=1:PointNums
    V = PointID(k);
    if IsVisited(k)==0
        ClusterID = ClusterID+1; % �صĸ�����1
        % �����߻�û�д����
        IsVisited(k) = 1;
        IDX(k) = ClusterID;        
        % ��ͨ��
        ConnectEdge = adjacent_edge(Graph,V,1);
        Neighbor = ConnectEdge(:,2);
        ind = Index(ismember(PointID,Neighbor));
        IsVisited(ind) = 1;
        IDX(ind) = ClusterID;
        % ��չ
        while ~isempty(Neighbor)
            V_Next = Neighbor(1);
            Neighbor(1) = [];
            ind = Index(PointID==V_Next);
            IsVisited(ind,:) = 1;
            IDX(ind,:) = ClusterID;
            ConnectEdge_next = adjacent_edge(Graph,V_Next,1);
            Neighbor_next = ConnectEdge_next(:,2);
            for u=1:length(Neighbor_next)
                ind = Index(PointID==Neighbor_next(u));
                if IsVisited(ind)==0
                    IsVisited(ind) = 1;
                    IDX(ind) = ClusterID;
                    Neighbor = [Neighbor;Neighbor_next(u)];
                end
            end
        end
    end    
end
% ɾ����С�Ĵأ���Ϊ��������
MINPTS = 5;  % ��С���ڰ�����ʵ����Ŀ����ֵ��һ������ʵ����С��MinPts����Ϊ����
if ClusterID>0
    clearnIDX = zeros(ClusterID,1);
    for i=1:(ClusterID)
        clusID=(IDX==i);
        if(sum(clusID)<=MINPTS)
            clearnIDX(i)=i;
        end
    end
    clearnIDX = clearnIDX(clearnIDX~=0);
    for i=1:length(clearnIDX)
        clus_id=clearnIDX(i);
        RemoveVertex = PointID(IDX==clus_id);
        for j=1:length(RemoveVertex)
            V=RemoveVertex(j);
            INX=(Graph(:,1)==V)|(Graph(:,2)==V);
            Graph(INX,:)=[];
        end
        IDX(IDX==clus_id)=0;
    end
    class = unique(IDX(IDX~=0),'rows');
    classnums = length(class);
    % [class,classnums] = remsame(IDX(IDX~=0));
    for c=1:classnums
        IDX(IDX==class(c))=c;
    end
end
Sub_Graph = cell(max(IDX),1);
for i=1:max(IDX)
    Sub_Graph{i}=[];
    for j=1:size(Graph,1)
        if sum(ismember(Graph(j,:),PointID(IDX==i)))>=1
            Sub_Graph{i} = [Sub_Graph{i}; Graph(j,:)];
        end
    end
end
end % // subgraph()



% ****************************************************************************
% **********************        global_cut()        **************************
% ****************************************************************************
%% global_cut()
function [SubGraph, Active_Edge] = global_cut(Graph,MVD,A)
% ͼ�и��߼��߳���ͼ�е��ID
PointID = unique([Graph(:,1);Graph(:,2)]);  % ����Graph�ĵ��ID�ţ���ԭʼ�ռ����ݼ��е������ţ�
PointNums = length(PointID);
% ���峤��
Global_Mean = mean(Graph(Graph(:,3)>MVD,3));
Global_Variation = std(Graph(Graph(:,3)>MVD,3));
Is_Remove = zeros(size(Graph,1),1);

Graph = [Graph, (1:size(Graph,1))'];
for i=1:PointNums
    Knn_Edge = adjacent_edge(Graph,PointID(i),1);
    Knn_Edge_Nums = size(Knn_Edge,1);
    if Knn_Edge_Nums>=1
        Local_Mean = mean(Knn_Edge(:,3));
        Global_Cut_Value = Global_Mean+A*(Global_Mean/Local_Mean)*Global_Variation;
        for j=1:Knn_Edge_Nums
            if (Knn_Edge(j,3) > MVD)&&(Knn_Edge(j,3) >= Global_Cut_Value)
                Is_Remove(Knn_Edge(j,4)) = 1;
            end
        end
    else
        continue;
    end
end
Active_Edge = Graph(Is_Remove==0,1:3);

clear Global_Mean ;
clear Global_Variation ;
clear Is_Remove ;
clear Knn_Edge ;
clear Knn_Edge_Nums;
clear Local_Mean;
clear Global_Cut_Value;

% �������ͨ��ͼ
IDX = zeros(PointNums,1);
IsVisited = zeros(PointNums,1);
Index = (1:PointNums);
% Vertex = PointID;
% VertexNums = length(Vertex);
ClusterID = 0;

for k=1:PointNums
    V = PointID(k);
    if IsVisited(k)==0
        ClusterID = ClusterID+1; % �صĸ�����1
        % �����߻�û�д����
        IsVisited(k) = 1;
        IDX(k) = ClusterID;        
        % ��ͨ��
        ConnectEdge = adjacent_edge(Active_Edge,V,1);
        Neighbor = ConnectEdge(:,2);
        ind = Index(ismember(PointID,Neighbor));
        IsVisited(ind) = 1;
        IDX(ind) = ClusterID;
        % ��չ
        while ~isempty(Neighbor)
            V_Next = Neighbor(1);
            Neighbor(1) = [];
            ind = Index(PointID==V_Next);
            IsVisited(ind,:) = 1;
            IDX(ind,:) = ClusterID;
            ConnectEdge_next = adjacent_edge(Active_Edge,V_Next,1);
            Neighbor_next = ConnectEdge_next(:,2);
            for u=1:length(Neighbor_next)
                ind = Index(PointID==Neighbor_next(u));
                if IsVisited(ind)==0
                    IsVisited(ind) = 1;
                    IDX(ind) = ClusterID;
                    Neighbor = [Neighbor;Neighbor_next(u)];
                end
            end
        end
    end    
end
% ɾ����С�Ĵأ���Ϊ��������
MINPTS = 5;  % ��С���ڰ�����ʵ����Ŀ����ֵ��һ������ʵ����С��MinPts����Ϊ����
if ClusterID>0
    clearnIDX = zeros(ClusterID,1);
    for i=1:(ClusterID)
        clusID=(IDX==i);
        if(sum(clusID)<=MINPTS)
            clearnIDX(i)=i;
        end
    end
    clearnIDX = clearnIDX(clearnIDX~=0);
    for i=1:length(clearnIDX)
        clus_id=clearnIDX(i);
        RemoveVertex = PointID(IDX==clus_id);
        for j=1:length(RemoveVertex)
            V=RemoveVertex(j);
            INX=(Active_Edge(:,1)==V)|(Active_Edge(:,2)==V);
            Active_Edge(INX,:)=[];
        end
        IDX(IDX==clus_id)=0;
    end
    class = unique(IDX(IDX~=0),'rows');
    classnums = length(class);
    % [class,classnums] = remsame(IDX(IDX~=0));
    for c=1:classnums
        IDX(IDX==class(c))=c;
    end
end
SubGraph = cell(max(IDX),1);
for i=1:max(IDX)
    SubGraph{i}=[];
    for j=1:size(Active_Edge,1)
        if sum(ismember(Active_Edge(j,:),PointID(IDX==i)))>=1
            SubGraph{i} = [SubGraph{i};Active_Edge(j,:)];
        end
    end
end
end % // global_cut()


% ****************************************************************************
% **********************        local_cut()         **************************
% ****************************************************************************
%% local_cut()
function [SubGraph, Active_Edge] = local_cut(Graph,MVD,B)
% ͼ�и��߼��߳���ͼ�е��ID
PointID = unique([Graph(:,1);Graph(:,2)]);  % ����Graph�ĵ��ID�ţ���ԭʼ�ռ����ݼ��е������ţ�
PointNums = length(PointID);
Graph = [Graph, (1:size(Graph,1))'];

% ������ͼGi�и�����ľֲ��߳�������
Local_Variation = zeros(PointNums,1);
Edges = cell(PointNums,1);
for i=1:PointNums
    K1_Edge = adjacent_edge(Graph,PointID(i),1); % i��1�������ڵı�
    Local_Variation(i) = std(K1_Edge(:,3));
    K2_Edge = adjacent_edge(Graph,PointID(i),2); % i��2�������ڵı�
%     Mean_GK2(i) = mean(K2_Edge(:,3));
    Edges{i} = K2_Edge; % Edges��Ŷ��ڵ��2�������ڵ�����������
    clear K1_Edge;
    clear K2_Edge;
end

% �ֲ�����
Mean_Variation = mean(Local_Variation);
clear Local_Variation;
Is_Remove = zeros(size(Graph,1),1);

for i=1:PointNums
    Gk2_Edge = Edges{i};
    Gk2_Edge_Nums = size(Gk2_Edge,1);
    if Gk2_Edge_Nums>=1
        Local_Mean = mean(Gk2_Edge(:,3)); % i��2�������ڵı߳���ֵ-Mean_G_K2
        Local_Cut_Value = Local_Mean+B*Mean_Variation;
        for j=1:Gk2_Edge_Nums
            if (Gk2_Edge(j,3) > MVD)&&(Gk2_Edge(j,3) >= Local_Cut_Value)
                Is_Remove(Gk2_Edge(j,4)) = 1;
            end
        end
    else
        continue;
    end
    clear Gk2_Edge;
end
Active_Edge = Graph(Is_Remove==0,1:3);
clear Edges;
clear Local_Mean ;
clear Local_Cut_Value ;
clear Is_Remove ;
clear Gk2_Edge_Nums;

% �������ͨ��ͼ
IDX = zeros(PointNums,1);
IsVisited = zeros(PointNums,1);
Index = (1:PointNums);
% Vertex = PointID;
% VertexNums = length(Vertex);
ClusterID = 0;

for k=1:PointNums
    V = PointID(k);
    if IsVisited(k)==0
        ClusterID = ClusterID+1; % �صĸ�����1
        % �����߻�û�д����
        IsVisited(k) = 1;
        IDX(k) = ClusterID;        
        % ��ͨ��
        ConnectEdge = adjacent_edge(Active_Edge,V,1);
        Neighbor = ConnectEdge(:,2);
        ind = Index(ismember(PointID,Neighbor));
        IsVisited(ind) = 1;
        IDX(ind) = ClusterID;
        % ��չ
        while ~isempty(Neighbor)
            V_Next = Neighbor(1);
            Neighbor(1) = [];
            ind = Index(PointID==V_Next);
            IsVisited(ind,:) = 1;
            IDX(ind,:) = ClusterID;
            ConnectEdge_next = adjacent_edge(Active_Edge,V_Next,1);
            Neighbor_next = ConnectEdge_next(:,2);
            for u=1:length(Neighbor_next)
                ind = Index(PointID==Neighbor_next(u));
                if IsVisited(ind)==0
                    IsVisited(ind) = 1;
                    IDX(ind) = ClusterID;
                    Neighbor = [Neighbor;Neighbor_next(u)];
                end
            end
        end
    end    
end
% ɾ����С�Ĵأ���Ϊ��������
MINPTS = 5;  % ��С���ڰ�����ʵ����Ŀ����ֵ��һ������ʵ����С��MinPts����Ϊ����
if ClusterID>0
    clearnIDX = zeros(ClusterID,1);
    for i=1:(ClusterID)
        clusID=(IDX==i);
        if(sum(clusID)<=MINPTS)
            clearnIDX(i)=i;
        end
    end
    clearnIDX = clearnIDX(clearnIDX~=0);
    for i=1:length(clearnIDX)
        clus_id=clearnIDX(i);
        RemoveVertex = PointID(IDX==clus_id);
        for j=1:length(RemoveVertex)
            V=RemoveVertex(j);
            INX=(Active_Edge(:,1)==V)|(Active_Edge(:,2)==V);
            Active_Edge(INX,:)=[];
        end
        IDX(IDX==clus_id)=0;
    end
    class = unique(IDX(IDX~=0),'rows');
    classnums = length(class);
    % [class,classnums] = remsame(IDX(IDX~=0));
    for c=1:classnums
        IDX(IDX==class(c))=c;
    end
end
SubGraph = cell(max(IDX),1);
for i=1:max(IDX)
    SubGraph{i}=[];
    for j=1:size(Active_Edge,1)
        if sum(ismember(Active_Edge(j,:),PointID(IDX==i)))>=1
            SubGraph{i} = [SubGraph{i};Active_Edge(j,:)];
        end
    end
end
end % // local_cut()


% ****************************************************************************
% **********************      adjacent_edge()       **************************
% ****************************************************************************
%% adjacent_edge()
function Edge = adjacent_edge(Graph,VertexID,K)
% ������һ����K��������Delaunay�������ı�
if nargin<3
    K=1;
end
% When K=1
% ��VertexIDΪ�ߵ����
Edge = Graph(Graph(:,1)==VertexID,:);
% ��VertexIDΪ�ߵ��յ�
CIndex = [2,1,(3:size(Graph,2))];
Edge = [Edge; Graph(Graph(:,2)==VertexID,CIndex)];
Edge = unique(Edge,'rows');
% When K>1
for i=1:(K-1)
    PointID = unique([Edge(:,1);Edge(:,2)]);
    PointNums = length(PointID);
    Edge = [];
    for j=1:PointNums   % ������PointID(j)ֱ�������ĵ��ID��
        % ��PointID(j)Ϊ�ߵ����
        Edge = [Edge;Graph(Graph(:,1)==PointID(j),:)];
        % ��PointID(j)Ϊ�ߵ��յ�
        Edge = [Edge;Graph(Graph(:,2)==PointID(j),CIndex)];   
    end 
    Edge = unique(Edge,'rows');
end
end % // adjacent_edge()
% ****************************************************************************



% ****************************************************************************
% **********************        ClusterID()         **************************
% ****************************************************************************
%%
function IDX = ClusterID(data,Clust)
IDX = zeros(size(data,1),1);
for i=1:size(Clust,1)
    PointID = Clust{i}(:,1:2);
    PointID = unique(PointID(:));
    IDX(PointID)=i;
    clear PointID;
end
end % // ClusterID()




