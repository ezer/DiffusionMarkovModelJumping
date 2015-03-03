% Xiaoyan Ma, Department of Genetics, Cambridge, UK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%and it should be a fixed value for input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are some examples for the input binding site:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Double Cluster with 3D distance of 200 nm

%  between_site=6.7;
%  posi_mat=[0 0 0;between_site 0 0;0 200 0;between_site 200 0];      
%%%%%%%%%6.7 is for center 20bp away;

%  cluster=[1 1 2 2]';
%  cluster=[1,1:(size(posi_mat,1)-1)];
%  binding_site=[posi_mat,cluster];

% Then the binding_site should look like that:

%binding_site =

  %       0         0         0    1.0000
  %  6.7000         0         0    1.0000
  %       0  200.0000         0    2.0000
  %  6.7000  200.0000         0    2.0000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Tetrahedron with clusters in each of its corner. The edge length is
% 100nm

% between_site=8.5   
% posi_mat=[0,2/3*sqrt(3),0;1,-sqrt(3)/3,0;-1,-sqrt(3)/3,0;0,0,2/3*sqrt(6)];
% posi_mat=posi_mat .* 50;
% posi_mat_2=[posi_mat(:,1)+between_site,posi_mat(:,2:3)];
% posi_mat=[posi_mat;posi_mat_2];
% cluster=[1 2 3 4 1 2 3 4]';
% binding_site=[posi_mat,cluster];
% binding_site=[binding_site(1,:);binding_site(5,:);binding_site(2,:);...
% binding_site(6,:);binding_site(3,:);binding_site(7,:);binding_site(4,:);...
% binding_site(8,:)];

% Then
%binding_site =

 %        0   57.7350         0    1.0000
 %   8.5000   57.7350         0    1.0000
 %  50.0000  -28.8675         0    2.0000
 %  58.5000  -28.8675         0    2.0000
 % -50.0000  -28.8675         0    3.0000
 % -41.5000  -28.8675         0    3.0000
 %        0         0   81.6497    4.0000
 %   8.5000         0   81.6497    4.0000


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. cube with a binding site cluster in each of its corner. Each cluster
% contains 3 binding sites.

%between_site=6.7;
%coordi=[0 0 0;0 1 0 ;1 0 0;1 1 0; 0 0 1;1 0 1;0 1 1;1 1 1];
%edge=200;
%coo_2=[(edge.*coordi(:,1)+6.7),edge.*coordi(:,2:3)];
%coo_3=[(edge.*coordi(:,1)+13.4),edge.*coordi(:,2:3)];
%posi_mat=(reshape([edge.*coordi';coo_2';coo_3'],[3,24]))';
%cluster=[1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8]';
%binding_site=[posi_mat,cluster];

%binding_site =

%         0         0         0    1.0000
%    6.7000         0         0    1.0000
%   13.4000         0         0    1.0000
%         0  200.0000         0    2.0000
%    6.7000  200.0000         0    2.0000
%   13.4000  200.0000         0    2.0000
%  200.0000         0         0    3.0000
%  206.7000         0         0    3.0000
%  213.4000         0         0    3.0000
%  200.0000  200.0000         0    4.0000
%  206.7000  200.0000         0    4.0000
%  213.4000  200.0000         0    4.0000
%         0         0  200.0000    5.0000
%    6.7000         0  200.0000    5.0000
%   13.4000         0  200.0000    5.0000
%  200.0000         0  200.0000    6.0000
%  206.7000         0  200.0000    6.0000
%  213.4000         0  200.0000    6.0000
%         0  200.0000  200.0000    7.0000
%    6.7000  200.0000  200.0000    7.0000
%   13.4000  200.0000  200.0000    7.0000
%  200.0000  200.0000  200.0000    8.0000
%  206.7000  200.0000  200.0000    8.0000
%  213.4000  200.0000  200.0000    8.0000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An example below:

posi_mat=[0 0 0;6.7 0 0;0 200 0;6.7 200 0];
cluster=[1 1 2 2]';
binding_site=[posi_mat,cluster];

% between_site value must be specified in nanometer which is the 1D
% distance between 2 binding site in the cluster. You can use 1bp=0.34nm to
% convert bp to nm.
between_site=6.7;

% jump_event_list function generates the long table of simulated times to reach 
% specific binding sites from each starting point
% along with the corresponding binding site label 

% Please specify the storage path for the out put file as the third arguement ( the PATH should be followed with a '/') and file name as
% the forth argument.
jump_event_list(binding_site,between_site,'/home/xm227/crm/matlab_script/data/','huge_table_cube_triple_sites_5bp_d200.txt');
