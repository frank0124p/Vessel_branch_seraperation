%Vessel_Separation
clc;close all;clear
%Add big branch radius and pair to the node
path=pwd;
path=strcat(pwd,'\Data');
addpath(path);

load center_pts.mat % parameter-result the 4&5 column are the connected pair 
load the_radius_matrix.mat %response  each node radius  from big branch
% 
load BranchNode.mat %key point detection
% load the_radius_matrix_from_origin_Dicom.mat  %response to the each node radius  from Dicom image
keypoint=BranchNode;

%Keypoint construct the radius matrix
[n,~]=size(keypoint);
the_node_radius_matrix=zeros(n,5);
for i=1:n
    node=keypoint(i,1:3);
    index=ismember(the_radius_matrix(:,1:3),node,'rows');
    index=find(index==1);
    the_node_radius_matrix(i,1:4)=keypoint(i,1:4); %4th column indicate the connect branch
    
    the_node_radius_matrix(i,5)=the_radius_matrix(index,4);% 5th column imply the corresponding radius
end

%  [C, L, U] = SpectralClustering(connect_matrix, 2, 2);


%Write the stl file for the keypoint
[n,~]=size(keypoint);
result=zeros(512,512,150);
for i=1:n
        if the_node_radius_matrix(i,5)>1
            result(the_node_radius_matrix(i,2),the_node_radius_matrix(i,1),the_node_radius_matrix(i,3))=1;
        end
end


figure,imshow(max(result,[],3));


fv = isosurface(result,0.1);  
% fv=smoothpatch(fv,1,5);

stlwrite('KeyPoint.stl', fv);



% % Plot the Branch Point
% [n,~]=size(keypoint);
% 
% for i=1:n
%     if the_node_radius_matrix(i,5)>4.5
%     plot3(the_node_radius_matrix(i,2),the_node_radius_matrix(i,1),the_node_radius_matrix(i,3),'.','Markersize',20,'MarkerEdgeColor','r');
% 
%     else
%         
%      plot3(the_node_radius_matrix(i,2),the_node_radius_matrix(i,1),the_node_radius_matrix(i,3),'.','Markersize',10,'MarkerEdgeColor','g');
%     end
%     rotate3d on
% 
%     hold on
%     axis([ 150 350 50 350 1 150])
%     view(150,30)
%     grid on
%     h=gca; 
% 
%     title('3D Branch Visualization','fontsize',14);
%     fh = figure(1);
%     set(fh, 'color', 'white'); 
%      F=getframe;
%      disp(i)                                                    
% end
% 
% movie(F)


