clear all;
clc;
% %% 6 user channel assignment
match_c_all_new=[1,4,2,5,3,6;1,4,2,6,3,5;1,5,2,4,3,6;1,5,2,6,3,4;1,6,2,4,3,5;1,6,2,5,3,4;1,2,4,5,3,6;1,2,4,6,3,5;1,3,4,5,2,6;1,3,4,6,2,5;2,3,4,5,1,6;2,3,4,6,1,5];
% match_c_all_new=[1,2,4,3;4,3,2,1;4,2,3,1;3,1,4,2;1,4,3,2;3,2,1,4];
channel=3;
gamma=0.5;
user=6;
max_c=2; %maximum users a channel contains
delta=1e-11; %delta^2=-110dB
Pmax=1.5;  %revise  !!!
delta=1e-11; %delta^2=-110dB

%channel=3;
M=20;   %%revise!!!cov1
mode=0; %1！！CR NOMA  0！！STAR-RIS NOMA 


[h_origin,opt_h]=channel_generate(M,channel,mode,user);
%load('channel_60_6');
for j=1:1
p=ones(channel,max_c)*Pmax/user;  %%revise !!!
%alpha=0.0036*ones(channel);
alpha=1*ones(channel);
match_c=match_c_all_new(j,:);
match_c=reshape(match_c,max_c,channel)';

match_c=[1,5;2,6;3,4];
maxitr2=10;
swap=zeros(1,maxitr2);
randomswap=zeros(1,maxitr2);
 
[do_old,tr,h,h_ini,Hc_ini2,a,b]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode); 
[do,~,R,T]=do_max_combined(tr,h,channel,max_c,M);
opt_H=h_ini;
% 
%  for i=1:channel
%     for j=i+1:channel
%         for m=1:max_c
%             for n=1:max_c 
%                 [swap_or_not,match_c,p,Hc,alpha,do,weak1,weak2,randomornot]=swap_many_2(user,match_c,p,Hc,do,opt_H,[i,m],[j,n],channel,alpha);
%             end
%         end
%     end
% end

% % %     do_ini=do;
% if z<1e-11&&(~isnan(z))

%% method 1
% [alpha,z]=feasibility_initial_point(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode,gamma);  
%  for itr1=1:3
%  [Hc,opt_H,alpha,R,T,opt_rate]=opt_matrix(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode,gamma); 
%  end
%  [p,alpha,opt_rate2]=opt_p_gp(user,match_c,do,Hc,Pmax,channel,gamma);


 

%% method 2
% for itr2=1:20
%     
% for itr1=1:3
% [alpha,opt_rate2,p,x]=opt_p_cub(do,alpha,Hc,gamma,channel,Pmax);
% end
% [opt_rate2,alpha,a,b,Hc]=opt_matrix_sca(M,h,channel,tr,do,a,b,x,p);
% end

%% method oma
% w=0.5*ones(channel,max_c);
% 
% for itr=1:3
% [opt_rate2,p,w]=opt_bandpower_oma(Hc,channel,Pmax,gamma);      
% [opt_rate2,a,b,Hc]=opt_matrix_sca_oma(M,h,channel,tr,a,b,p,w,gamma);  
% opt_value(itr)=sum(sum(opt_rate2));
% end
% 
% for i=1:channel
%     for j=1:max_c
%         opt_rate_ini(i,j)=0.5*log2(2*Hc_ini2(i,j)*0.25*1e11);
%     end
% end



%% 
%  for i=1:channel
%     for j=i+1:channel
%         for m=1:max_c
%             for n=1:max_c 
%                 [swap_or_not,match_c,p,Hc,alpha,do,weak1,weak2,randomornot]=swap_many_2(user,match_c,p,Hc,do,opt_H,[i,m],[j,n],channel,alpha);
%                 if swap_or_not==1              
%                 [do_old,tr,h]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode);             
%                  for itr1=1:2
%                 [Hc,opt_H,alpha,R,T,opt_rate]=opt_matrix(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode); 
%                  end
%                 [p,alpha,opt_rate2]=opt_p_gp(user,match_c,do,Hc,Pmax,channel);
%                 end
%             end
%         end
%     end
% end





% else
%       opt_rate2=zeros(channel,max_c);
%       opt_value2(j,i)=sum(sum(opt_rate2));
% end
end









