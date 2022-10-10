clear all;
clc;
% % %% 6 user channel assignment
match_c_all_new=[1,4,2,5,3,6;1,4,2,6,3,5;1,5,2,4,3,6;1,5,2,6,3,4;1,6,2,4,3,5;1,6,2,5,3,4;1,2,4,5,3,6;1,2,4,6,3,5;1,3,4,5,2,6;1,3,4,6,2,5;2,3,4,5,1,6;2,3,4,6,1,5];
%match_c_all_new=[1,2,4,3;4,3,2,1;4,2,3,1;3,1,4,2;1,4,3,2;3,2,1,4];
channel=3;
user=6;
max_c=2; % maximum users a channel contains
delta=1e-11; %delta^2=-110dB
Pmax=2;  %revise  !!!
p=ones(channel,max_c)*Pmax/user;  %%revise !!!
alpha=0.0036*ones(channel); % feasiblity initialization
alpha=0.05*ones(channel); % feasiblity initialization
gamma=0.1;
%channel=3;
M=60;   %% revise!!!cov1
mode=0; % 1！！CR  0！！STAR
for i=1:500
[h_origin,opt_h]=channel_generate(M,channel,mode,user);
for j=1:12
%% STAR NOMA three steps
% p=ones(channel,max_c)*Pmax/user;  %%revise !!!
% match_c=match_c_all_new(j,:);
% match_c=reshape(match_c,max_c,channel)';
%[do,tr,h]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode); 
% %[swap_or_not,match_c,p,Hc]=swap_many_3(user,match_c,p,h,do,opt_h,[1,1],[1,1],channel);
% %[alpha,z]=feasibility_initial_point(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode,gamma)
%  for itr2=1:3
% [Hc,opt_H,alpha,R,T,opt_rate2]=opt_matrix(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode,gamma); 
%  end
% [p,alpha,opt_rate2]=opt_p_gp(user,match_c,do,Hc,Pmax,channel,gamma);
% cdf_6user(j,i)=sum(sum(opt_rate2));


%% STAR OMA SCA
match_c=match_c_all_new(j,:);
match_c=reshape(match_c,max_c,channel)';
[do_old,tr,h,h_ini,Hc,a,b]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode); 
p=ones(channel,max_c)*Pmax/user;  %%revise !!!
[do,tr,h]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode);
[swap_or_not,match_c,p,Hc,alpha,do,weak1,weak2,randomornot]=swap_many_2(user,match_c,p,Hc,do,opt_h,[1,1],[1,1],channel,alpha);
for itr2=1:2
[Hc,opt_H,R,T,opt_rate2]=opt_matrix_oma(match_c,tr,h,h_origin,p,M,channel,gamma); 
[p,opt_rate2,w]=opt_p_oma(user,match_c,Hc,Pmax,channel,mode,gamma);
end
oma2_cdf2(j,i)=sum(sum(opt_rate2));
end
end
