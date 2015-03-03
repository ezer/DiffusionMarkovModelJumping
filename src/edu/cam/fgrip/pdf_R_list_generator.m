%%% Xiaoyan Ma, Department of Genetics, Cambridge, UK

%%% This function calculates the probability density for a diffusin molecule 
%%% to reach certain absorbing sphere at specific distances at each time point 
%%% given the time rsolution and time limit

function pdf_l=pdf_R_list_generator(R_list,save_path,file_name1,file_name2,resolution,limit)

time=(1/resolution):(1/resolution):limit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  
% normally the following lines is of no use, unless
% you want to set different resolution values for different time period in
% order to make the computation faster

%time_highReso= [0.0001:0.0001:1];
%time_lowReso=[1.01:0.01:100];

%%%% for time in time_lowReso, it will repeat the last point in the 0.0001reso,
%%%% and repeat it 100 times later. (See line 50-59)

%time=[time_highReso time_lowReso];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


% Default effective diffusion coefficient
D=3e5;

pdf_R_list=[];
integr=[];

% calculate the probability density function for each distance values
% stored in R_list
for k=1:size(R_list,1)
  
  pdf=[];
  R=R_list(k,1);
  
  % s is the capture distance for the target which equals to the sliding 
  % window length
  s=R_list(k,2);
  for t=time
    
    % below is the probability density function for a diffusin molecule 
    % to reach certain absorbing sphere at specific distances over time
    fun=@(t,D,R,s) s.*(R-s)./(2.*R.*sqrt(pi().*D.*t^3)).*exp(-((R-s).^2)./(4.*D.*t));
    pdf=[pdf fun(t,D,R,s)];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normally the following lines is of no use, unless
% you want to set different resolution values for different time period in
% order to make the computation faster
  
  % to_rep=pdf((1e4+1):length(pdf));
  %to_rep_new=repmat(to_rep,100,1);
  %to_rep_new=to_rep_new(1:end);
  %pdf=[pdf(1:1e4), to_rep_new];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  pdf_R_list=[pdf_R_list;pdf];
  
  %integr=[integr;[sum(pdf(1:resolution*limit))/resolution,...
  %sum(pdf(1:resolution*2))/resolution,sum(pdf(1:resolution*3))/resolution]];
  
  % calculate the integral of this probability density function over the
  % time limit
  integr=[integr;sum(pdf(1:resolution*limit))/resolution];
  
end

pdf_l=pdf_R_list;
R_list_with_inte=[R_list,integr];

save(fullfile(save_path,file_name1),'pdf_R_list');
save(fullfile(save_path,file_name2),'R_list_with_inte');

end

