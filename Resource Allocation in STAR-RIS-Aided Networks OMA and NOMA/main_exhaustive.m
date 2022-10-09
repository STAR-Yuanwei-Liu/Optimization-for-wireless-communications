clear all;
clc;
max_c=2; %maximum users a channel contains
channel=4;
%channel=3;
M=20;   %%revise!!!cov1
mode=0; %1！！CR  0！！STAR
user=8;
%user=6
Pmax=user*0.1;
delta=1e-11; %delta^2=-110dB
%load('channel_40_8','h_origin','opt_h');
load('channel_20_8_tr','h_origin','opt_h');    %%revise!!!cov1
%match_c=[1,7;2,3;4,6;5,8;9,10]; %4 steps converge
%match_c=[10,8;2,3;4,6;5,7;9,1];
match_c=[3,1;4,2;6,5;7,8];
match_c=[1,3;4,5;7,2;6,8];
match_c=[1,8;4,5;2,7;6,3];
match_c=[1,2;3,4;5,6;7,8];
match_c=[2,5;8,1;4,6;7,3]; %%%%maximum
%match_c=[6,7;10,8;4,5;2,1;3,9];
p=ones(channel,max_c)*0.1;
alpha=0.0036*ones(channel,1);
itr1=0;
itr2=0;
swap_or_not=0;
num_total=0;
total_swap=0;
maxitr2=15;
swap=zeros(1,maxitr2);
swap(1)=1; 
[do,tr,h]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode);
for itr2=2:maxitr2
if swap(itr2-1)>=1||(itr2==1)
    p=ones(channel,max_c)*0.1;
    alpha=0.0036*ones(channel,1);
    for itr1=1:3
    %do_temp=do;
    [~,tr,h]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode);
    %do=do_temp;
    [Hc,opt_H,alpha,R,T,opt_rate,do]=opt_matrix(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode); %do?
     end
    [p,alpha,opt_rate2]=opt_p_gp(user,match_c,do,Hc,Pmax,channel);
    str=num2str(match_c);  %%if swap, display current match
    disp(str);
    sum_rate(itr2)=sum(sum(opt_rate2));
    swap(itr2)=0;
    %% swaping
flag=0;
for i=1:channel
    for j=i+1:channel
        for m=1:max_c
            for n=1:max_c 
                %if match_c(i,mod(m,2)+1)~=match_c(j,n)&&match_c(j,mod(m,n)+1)~=match_c(i,m)
                    %[swap_or_not,match_c]=swap_many(user,match_c,p_,Hc,do,opt_H,[i,m],[j,n],channel);
                    [swap_or_not,match_c,p,Hc]=swap_many_3(user,match_c,p,Hc,do,opt_H,[i,m],[j,n],channel);
                    if(swap_or_not)==1
%                         str=num2str(match_c);  %%if swap, display current match
%                         disp(str);
                        total_swap=total_swap+1;
                        swap(itr2)=swap(itr2)+1;
                        flag=1;
                    else
                        swap(itr2)=swap(itr2)+0;
                    end
                    num_total=num_total+1;
                %end 
            end
        end
    end

end
%     temp=match_c';
%     match_state(itr2,:)=temp(1:end);
else
    sum_rate(itr2)=sum_rate(itr2-1);
end
end
%% results


[V1,D1]=eig(R);
[V2,D2]=eig(T);