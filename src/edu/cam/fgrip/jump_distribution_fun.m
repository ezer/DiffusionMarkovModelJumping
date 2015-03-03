% Xiaoyan Ma, Department of Genetics, Cambridge, UK

%%%jump_distribution_fun simulates the jumping event from a specific
%%%starting point to different target sites. It is a function called by
%%%the jump_event_list function. This function returns an 1*2 array of 
%%%the time and the next target a TF reaches when it starts from a specific
%%%position according to diffusion time distribution. If the output is [-1
%%%-1], it means the TF will not reach any specific target before it
%%%returns to the TF pool. The time threshod for a TF to return to the pool
%%%is set to be 10 second. Please refer to the Supplementary material for
%%%detailed reason of choosing that threshod. 

function time_nextEvent=jump_distribution_fun(binding_site,posi_from,R_list_with_inte,pdf_R_list,resolution,limit,pdf_nearby_site,R_list_with_inte_nearby,binding_site_num_detail,binding_site_detailed)

site_number=size(binding_site,1);

%t_reach is the matrix storing the time for reaching each site and the site
%number.
t_reach=repmat(-1,(site_number-1),2);

%R_list_each is the list of distance between the starting point and each of
%the other clusters. In the following lines(line 21 to 46), it simulates for
% the potential diffusion between different clusters. 
R_list_each=R_list_each_generator(binding_site,posi_from);

for i=1:size(R_list_each,1)
    row_number_in_pdf_R_list = find(abs(R_list_with_inte(:,1)-R_list_each(i,2)) < eps & R_list_with_inte(:,2)==R_list_each(i,3));
    pdf_t=pdf_R_list(row_number_in_pdf_R_list,:) ;   
    
    %For each potential target cluster, generating a random number and comparing it 
    %with the cumulative probability for the TF to reach this cluster to
    %decide if this event can happen or not.
    if(rand< R_list_with_inte(row_number_in_pdf_R_list,3))
      
      % Generate the random time point according to the probility density
      % function of absorbing by a specific target sphere. 
      t_reach(i,1)=randpdf(pdf_t(1:(resolution*limit)),[(1/resolution):(1/resolution):(1*limit)],[1,1]);
      t_reach(i,2)=R_list_each(i,1);
    
    else
      t_reach(i,1)=100;
      
      %If it is set to be 100, it means in this single simulation, the TF 
      %doesn't reach this target cluster before returning to the pool.
      t_reach(i,2)=R_list_each(i,1);
    
    end 
end

%From line 50 to line 75 are the simulation for little jumps between
%different binding sites within the same homotypic cluster.

little_jump_site=binding_site_detailed(binding_site_detailed(:,4)==posi_from,[1:3,5]);
little_jump_site=[little_jump_site,repmat(2,size(little_jump_site,1),1)];

%R_list_each_little is the list of distance between the starting point and each of
%the other binding sites in the same homotypic cluster.
R_list_each_little=R_list_each_generator(little_jump_site,find(little_jump_site(:,4)==binding_site_num_detail));
t_reach_near=repmat(-1,size(R_list_each_little,1),2);

for i=1:size(R_list_each_little,1)
    nearby_site_num=R_list_each_little(i,1);
    
    row_number_in_pdf_nearby_site=find(abs(R_list_each_little(i,2)-R_list_with_inte_nearby(:,1))<0.000001);
    pdf_t_nearby=pdf_nearby_site(row_number_in_pdf_nearby_site,:) ;  
    
    if(rand<R_list_with_inte_nearby(row_number_in_pdf_nearby_site,3))
     t_reach_near(i,:)=[randpdf(pdf_t_nearby(1:1e5),[1e-6:1e-6:0.1],[1,1]),nearby_site_num+100];
   else
     t_reach_near(i,:)=[100,nearby_site_num+100];
     % using nearby_site_num+100 as the symbol for the spefic binding site
     % within the cluster is to distinguish that from the cluster number.
     % It will be restored to the real binding site number later on in the 
     %jump_event_list function.
     
    end
end
t_reach=[t_reach;t_reach_near];
  
% normally -1 will not appear, unless randpdf generate NA;
  have_value=t_reach(t_reach(:,1)~= -1 & t_reach(:,1)~=100 ,:);
  if(~isempty(have_value))
      
     % To decide which binding site the TF goes to, we choose the binding 
     %site with the minimum simlulated time, because the TF can only go to
     %one binding site at each time. The time distribution of the TF to go
     %to a specific site is actually conditioned on it doesn't bump into
     %others.
     j=find(t_reach(:,1)==min(have_value(:,1)));
     site_next_reaction=t_reach(j,2);
     
     % if hope to change the time limit for keeping the TF in the box, 'limit' can be replaced by other values
     if (t_reach(j,1)<limit)  
        time_nextEvent= [site_next_reaction t_reach(j,1)];
     
     else 
        time_nextEvent= [-1 -1];
        %disp(t_reach(j,1));
        %disp(j);
        %disp(t_reach);
     end
         
  else
      time_nextEvent=[-1 -1];
      %-1 -1 means the TF will not reach any specific target before it
      %returns to the TF pool; 
  end
end