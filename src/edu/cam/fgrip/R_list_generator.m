% Xiaoyan Ma, Department of Genetics, Cambridge, UK

% output of this function is the matrix storing the possible distances 
%between each binding site cluster in the first colunm and the
%corresponding cluster radius in the second column.

function R_list=R_list_generator(binding_site,save_path,file_name)
r_list=[];   

%calculate the 3D distance between each input point stored in the input matrix.
for i=1:size(binding_site,1)
      for j=1:size(binding_site,1)
       if(binding_site(j,4) ~= binding_site(i,4)),
          R=sqrt((binding_site(i,1)-binding_site(j,1))^2+(binding_site(i,2)-binding_site(j,2))^2+(binding_site(i,3)-binding_site(j,3))^2);
          r_list=[r_list; [R,binding_site(j,5)]];
       end 
      end 
end

%delete the duplicate rows
r_list=unique(r_list,'rows');
R_list=r_list;
save(fullfile(save_path,file_name),'R_list');
end
