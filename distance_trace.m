function dist_matrix=distance_trace(the_matrix,seed_index)
%Set a graph with all node and connect edge
nodes=zeros(length(the_matrix),4);
nodes(:,1)=the_matrix(:,4);
nodes(:,2:4)=the_matrix(:,1:3);
edges=zeros(length(the_matrix),3);
edges(:,1)=1:length(the_matrix);
edges(:,2:3)=the_matrix(:,4:5);

%Calculate the real distance matrix
dist_matrix=cell(length(the_matrix),length(seed_index));
    for i=1:length(the_matrix)
        start_idx=the_matrix(i,4);
                for j=1:length(seed_index)
                    finish_index=seed_index(j); 
                     [dist,path] = dijkstra(nodes,edges,start_idx,finish_index);   
                      field = 'dist';
                      field2 = 'path';
                      data= struct(field,dist,field2,path);
                      dist_matrix{i,j}=data;               
                end
                
        disp(i);       
 
    end
end
