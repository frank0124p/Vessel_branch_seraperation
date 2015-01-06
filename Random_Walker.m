clear;clc;close all;
path=pwd;
path=strcat(pwd,'\Data');
addpath(path);
%--------Random Walker ---------
load no_zero.mat  %load all node in the result---- parameter: no_zero
load smooth_Ori.mat %The Dicom image
load the_radius_matrix.mat %response  each node radius  from big branch
load center_pts.mat % parameter-result the 4&5 column are the connected pair from original data
load the_radius_matrix.mat %response  each node radius  from big branch
load seed.mat %seed_node for sedd and idx for label
load dist_matrix.mat %trace the path about the real route
%Collect the Spatial ,Radius ,intensity infomation to construct the

%-----------------------Seed point processing------------------------------
seed_index=zeros(size(idx));
for i=1:length(root_node)
 index=ismember(the_radius_matrix(:,1:3),root_node(i,:),'rows');
 index=find(index==1);
 seed_index(i,1)=index;
end

seed_matrix=zeros(length(root_node),2);
seed_matrix(:,1)=seed_index;
seed_matrix(:,2)=idx;

%Reduce the seed pts 
number_seed=20;
[label1,~]=find(seed_matrix(:,2)==1);
[label2,~]=find(seed_matrix(:,2)==2);
seed=zeros(number_seed,2);
for i=1:number_seed
    if mod(i,2)==0
        seed(i,1)=seed_matrix(label2(i),1);
        seed(i,2)=seed_matrix(label2(i),2);
    else
        seed(i,1)=seed_matrix(label1(i),1);
        seed(i,2)=seed_matrix(label1(i),2);
    end
end

%adjacency matrix 
the_data_matrix=zeros(size(center_pts));
the_data_matrix(:,1:4)=the_radius_matrix;%Data member (x,y,z,radius,intensity)
for i=1:length(no_zero)
 no_zero(i,4)=smooth_Ori(no_zero(i,2),no_zero(i,1),no_zero(i,3));
end

for i=1:length(the_radius_matrix)
 index=ismember(no_zero(:,1:3),the_radius_matrix(i,1:3),'rows');
 index=find(index==1);
 the_data_matrix(i,5)=no_zero(index,4);
end

the_new_data_matrix=zeros(length(center_pts),7);
the_new_data_matrix(:,1:5)=center_pts;
the_new_data_matrix(:,6)=the_data_matrix(:,4);%for radius
the_new_data_matrix(:,7)=the_data_matrix(:,5);%for intensity value

%---------Calculate the real route for the node distance----------------

% dist_matrix=distance_trace(the_new_data_matrix,seed(:,1));
% save dist_matrix.mat dist_matrix
%--------------------------Construct the adjacency matrix -----------------
%weight matrix construct
weight_matrix=construct_adjacency_matrix(the_new_data_matrix);%get weight_matrix
%Degree matrix construct
% weight_sum=sum(weight_matrix);
% degree_matrix=zeros(size(weight_matrix));
% for i=1:length(weight_sum)
%     degree_matrix(i,i)=weight_sum(i);
% end

%Laplacian_matrix construct
% Laplacian_matrix=degree_matrix-weight_matrix;
% Laplacian_matrix=abs(Laplacian_matrix);

% Laplacian_matrix(Laplacian_matrix<0)=0;

%-------------------------Calculating the probability-------------------
% %Set up Dirichlet problem
% number_labels=2;
% boundary=zeros(length(seed_matrix),number_labels);
% for k=1:length(seed_matrix)
%     label=seed_matrix(k,2);
%     if label==1
%            boundary(k,1)=1;
%     else
%             boundary(k,2)=1;
%     end
% end
% 
% Laplacian_matrix = sparse(Laplacian_matrix);
% probabilities=dirichletboundary(Laplacian_matrix,seed_matrix(:,1),boundary);
% 
% final_label=zeros(length(probabilities),1);
% for i=1:length(probabilities)
%     temp1=probabilities(i,1);
%     temp2=probabilities(i,2);
%     
%     if temp1>temp2
%         final_label(i)=1;
%     else
%         final_label(i)=2;
%     end
% end

%-------the new version to calculate the probability about the reachability

% probability=zeros(length(the_data_matrix),length(seed));
% for i=1:length(the_data_matrix)
%     for j=1:length(seed)
%         probability(i,j)=weight_matrix(i,seed(j,1));
%     end
% end
% 
% final_label=zeros(length(probability),1);
% for i=1:length(probability)
%     [~,index] = max(probability(i,:),[],2);
%     final_label(i)=seed(index,2);
%    
% end


%-----------Distance only seperation from dist_matrix------------
seed_distance=zeros(length(dist_matrix),number_seed);
for i=1:length(dist_matrix) 
    for   j=1:number_seed
            seed_distance(i,j)=dist_matrix{i,j}.dist;
   end
end

final_label=zeros(length(seed_distance),1);
for i=1:length(seed_distance)
    [~,index] = min(seed_distance(i,:),[],2);
    final_label(i)=seed(index,2);   
end

%---------------write to the STL File----------------

result=zeros(512,512,150);
result2=zeros(512,512,150);

for i=1:size(the_data_matrix,1)
    tag=final_label(i);
    if tag==1
       result(the_data_matrix(i,2),the_data_matrix(i,1),the_data_matrix(i,3))=1;
    else
       result2(the_data_matrix(i,2),the_data_matrix(i,1),the_data_matrix(i,3))=1;
    end
end



fv = isosurface(result,0.1);  
stlwrite('result_back.stl', fv);

fv2 = isosurface(result2,0.1);  
stlwrite('result_forward.stl', fv2);

