% Xiaoyan Ma, Department of Genetics, Cambridge, UK

% This function generates the long table of simulated times to reach 
% specific binding sites from each starting point
% along with the corresponding binding site label  

% 1. The input matrix should have the same number of binding site in each
% cluster and the binding site in each cluster should follow each other as
% is shown in the example. 
% 2. the 3D postion should be in the scale of nanometer.
% 3. There is no need to label each binding site, because the
% jump_event_list function will do it. However, instead, the cluster number
% label is required.
% 4. input binding_site matrix should have the order of increasing x, y or 
% z in each cluster in order to make the calculation of cluster radius 
% correct;
% 5. both save_path and file_name are all character strings with quote.
% 6. between_site is the 1D distance measured in nm ( not in bp) 
% between binding sites next to each other in each cluster, 
%and it should be a fixed value for input. If you want
% non-fixed values, you also need to change the line of 97 following the
% instruction there. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jump_event_list(binding_site, between_site, save_path,file_name)

% label each binding site with a specific number
binding_site=[binding_site,(1:size(binding_site,1))'];

%resolution for calculating probability density function values for the
%diffusing particle being absorbed by the sphere over time. I set the default
%value to be 1e4, which means the time resolution is 1/1e4 second=0.0001
%second. Higher resolution can make the program running slower. 

resolution=1e4;

%Default time limit for simulation is 10sec. See detailed reason in the
%Supplementary material.
limit=10;

%cluster label
cluster=binding_site(:,4);

site_number=size(binding_site,1);

% The input matrix should have the same number of binding site in each
% cluster.
site_num_each_cluster=site_number/length(unique(cluster));

% binding_site_brief lists only the center position in every cluster with
% the cluster number and cluster absorbing radius.

binding_site_brief=[];

% for each cluster, calculate the cluster absorbing radius and the center of
% each cluster
for cluster_num=unique(cluster)',
    select_posi_in_cluster=find(binding_site(:,4)==cluster_num,1,'first');
    last_posi_in_cluster=find(binding_site(:,4)==cluster_num,1,'last');
    
    % distance between the first and the last binding site in the cluster
    d_between=sqrt((binding_site(select_posi_in_cluster,1)-binding_site(last_posi_in_cluster,1))^2+(binding_site(select_posi_in_cluster,2)-binding_site(last_posi_in_cluster,2))^2+(binding_site(select_posi_in_cluster,3)-binding_site(last_posi_in_cluster,3))^2);
    
    % s is the cluster absorbing radius; for only one binding site, the
    % default absorbing radius is 30nm, which corresponds to the definition
    % of sliding window length in FastGrip.
    s=d_between/2+30;
   
    center_in_cluster=(binding_site(select_posi_in_cluster,1:3)+binding_site(last_posi_in_cluster,1:3))/2;
    binding_site_brief=[binding_site_brief;[center_in_cluster,cluster_num,s]];
end

% calculate the distance between every pair of binding site clusters.
R_list=[];
R_list=R_list_generator(binding_site_brief,'/home/xm227/crm/matlab_script/data','R_list_try.mat');
R_list=unique(R_list,'rows');

% generate the probility density distribution for the diffusion particle
% being absorbed by a sphere given the specific radius of absorbing sphere around
% each binding site cluster.
pdf_R_list=pdf_R_list_generator(R_list,'/home/xm227/crm/matlab_script/data','pdf_try.mat','R_list_with_inte_try.mat',resolution,limit);

% R_list_with_inte is the abreveation for
% R_list_with_integral_of_probability_density. It stores the distances
% between each pair of binding site clusters along with the cluster
% absorbing radius and the integral for the absorbing probability density
% within the time limit.
load_path=strcat(save_path,'R_list_with_inte_try.mat');
to_load=load(load_path);
R_list_with_inte=to_load.R_list_with_inte;

%prob=max(R_list_with_inte(:,3));
%disp(prob);

% R_little_list stores distances between different binding sites
% within the same cluster. Our default setting assumes equal distance
% between two binding site next to each other. If you want to try unequal
% distances, please type the specific distances below instead of using 
%[1:(site_num_each_cluster-1)]'*between_site . 
R_little_list=[[1:(site_num_each_cluster-1)]'*between_site,repmat(2.0,(site_num_each_cluster-1),1)];

% generate the probility density distribution for the diffusion particle
% being absorbed by a nearby binding site in the same cluster
pdf_nearby_site=pdf_R_list_generator(R_little_list,save_path,'pdf_nearby.mat','R_list_with_inte_nearby.mat',1e6,0.1);

load_path=strcat(save_path,'R_list_with_inte_nearby.mat');
to_load=load(load_path);
R_list_with_inte_nearby=to_load.R_list_with_inte;

jump_struct=struct('start_site',[],'site_reached_and_time',[]); 
posi_from_list= 1:size(binding_site_brief,1);


% simulate the potential target one TF can reach starting from each different cluster
% and different binding sites in each cluster.
for posi_from=posi_from_list,
  for j=1:(length(cluster)/length(unique(cluster))) 
      
      %in cluster stores all binding site information within a specific cluster
      %which is the starting point of diffusion for a TF
      in_cluster=binding_site(binding_site(:,4)==posi_from,:);
      
      % binding_site_num_detail is the specific binding site number
      binding_site_num_detail=in_cluster(j,5);
      t_next_list=[];
      
      % 1e4 is the default row number in the output simulation result table
      for i=[1:1e4],
      
         % pick up the time to reach a specific target and the
         % corresponding target label according to the series of probability distribution
         % for a TF being absorbed by different targets.
         t_next=jump_distribution_fun(binding_site_brief,posi_from,R_list_with_inte,pdf_R_list,resolution,limit,pdf_nearby_site,R_list_with_inte_nearby,binding_site_num_detail,binding_site);
       
         t_next_list=[t_next_list;t_next];
      end
      
   % Store the simulated diffusion time and the target it reaches along
   % with the starting point information in t_next_struct.
   t_next_struct=struct('start_site', binding_site_brief(posi_from,:), 'site_reached_and_time', t_next_list);
   jump_struct=[jump_struct;t_next_struct];
  end
end

% delete the first empty element in jump_struct
jump_struct(1)=[];
%jump_struct_2=jump_struct;

% Since in previous steps, jump_struct only stores the time to reach each
% cluster for the diffusion between different TF clusters, we are going to
% specify the particular binding site TF reaches, given there is equal
% possibility for the TF to reach every binding site within the same
% cluster. As for little jumps within the same TF cluster, we already
% simulate the specific site it reaches when calling the function of
% jump_distribution_fun.
for posi_from=1:length(jump_struct) 
    
    t_next=jump_struct(posi_from).site_reached_and_time;
    t_next=[t_next,repmat(-1,size(t_next,1),1)];
   
    unique_cluster=unique(cluster)';
    
    %find the cluster number for a specific binding site.  
    cluster_number=floor((posi_from-1)/(length(cluster)/length(unique(cluster))))+1;
    
    posi_to_list=unique_cluster(unique_cluster~=jump_struct(posi_from).start_site(:,4));
    
    % for each potential target cluster
    for posi_to=posi_to_list,
        
        % dividing_factor is 1 over the binding site number in each cluster
        dividing_factor=1/length(find(binding_site(:,4)==posi_to));
        
        % find the binding sites in a specific cluster
        single_binding_site_list=binding_site(binding_site(:,4)==posi_to,5);
        
        % if there is only one binding site in each cluster
        if (dividing_factor==1)
            t_next(t_next(:,1)==posi_to,3)=single_binding_site_list;
            
        else 
            
          % for more than 1 site in each cluster case,  specify the 
          % particular binding site TF reaches, given there is equal
          % possibility for the TF to reach every binding 
          % site within the same cluster.
          assigning_mat=[[0:dividing_factor:(1- dividing_factor)]',(dividing_factor:dividing_factor:1)', single_binding_site_list];
      
          for row_number=find(t_next(:,1)==posi_to)',
             rand_num=rand;
             single_site_reach=assigning_mat(rand_num>assigning_mat(:,1)& rand_num<assigning_mat(:,2),3);
             t_next(row_number,3)= single_site_reach;
             %disp(single_site_reach);
          end
        end
    end
    
    % In order to avoid confusing the label between cluster and specific
    % binding site nearby, in previous steps for simulation, the labels for
    % the nearby binding sites have been turned to their original label
    % plus 100, so now we need to change them back.
    to_nearby=find(t_next(:,1)>=100);
    t_next(to_nearby,3)=t_next(to_nearby,1)-100;
    
    jump_struct(posi_from).site_reached_and_time=t_next;
    
end

% format the jump_struct to a printable version.
jump_struct_print=[];
for i=1:length(jump_struct)
   jump_struct_print=[jump_struct_print,jump_struct(i).site_reached_and_time(:,2:3)];
end

% format the printed version and specify the number of floating points to 
% print.
rep_print=repmat(['%12.7f %12.f'],1,length(jump_struct));
rep_print=[rep_print,['\n']];

full_path=strcat(save_path,file_name);
fileID = fopen(full_path,'w');
fprintf(fileID,rep_print,jump_struct_print');
fclose(fileID);

% The output long table have the colunm number equals to 2*binding site
% number and the row number of 10000

end

