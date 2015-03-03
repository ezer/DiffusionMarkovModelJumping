% Xiaoyan Ma, Department of Genetics, Cambridge, UK

% This function calculate the 3D distance between a specific point and each
% of the other points stored in the input matrix.
% posi_from is the row number for the starting point;
% binding_site is an input matrix having the same format and column content
% as the binding_site_brief matrix in jump_event_list function.
function R_list_each=R_list_each_generator(binding_site,posi_from)
r_list_each=[];   
i=posi_from;

%calculate the 3D distance between posi_from and other input points stored in the input matrix.
for j=1:size(binding_site,1)
       if(binding_site(j,4) ~= binding_site(i,4)),
          R=sqrt((binding_site(i,1)-binding_site(j,1))^2+(binding_site(i,2)-binding_site(j,2))^2+(binding_site(i,3)-binding_site(j,3))^2);
          r_list_each=[r_list_each; [binding_site(j,4),R,binding_site(j,5)]];
       end 
 end 
R_list_each=r_list_each;
end

