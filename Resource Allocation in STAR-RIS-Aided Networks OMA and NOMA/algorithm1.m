function [Hc,R,T,opt_rate2,match_c,swap,optvalue]=algorithm1(M,channel,mode,h_origin,opt_h)
max_c=2; %maximum users a channel contains
%channel=3;
user=6;
%user=6
Pmax=user*0.1;  %revise  !!!
delta=1e-11; %delta^2=-110dB

%load('channel_40_8_tr','h_origin','opt_h');    %%revise!!!cov1
%match_c=[1,2;3,4;5,6;7,8];
% match_c=[1,8;4,5;2,7];
match_c=[1,6;2,5;3,4];
%match_c=[1,6;3,2;5,4];
%match_c=[2,5;8,1;4,6;7,3]; %%%%maximum
%match_c=[8,5;2,1;4,6;7,3]; %%%%
p=ones(channel,max_c)*0.1;  %%revise !!!
alpha=0.0036*ones(channel,1);
itr1=0;
itr2=0;
swap_or_not=0;
num_total=0;
total_swap=0;
maxitr2=30;
swap=zeros(1,maxitr2);
randomswap=zeros(1,maxitr2);
swap(1)=1; 
[do,tr,h]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode); 
if mode==0||mode==1
[alpha]=feasibility_initial_point(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode);   
end

if mode==2  %%STAR-OMA
[Hc,opt_H,R,T,opt_rate]=opt_matrix_oma(match_c,tr,h,h_origin,p,M,channel); %do?
[p_,opt_rate2]=opt_p_oma(user,match_c,Hc,Pmax,channel);  
optvalue=sum(sum(opt_rate2));
end


for itr2=2:3
if swap(itr2-1)>=1||(itr2==1)||mode==2
    if mode==0||mode==1
    for itr1=1:3
    [~,tr,h]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode);
    %do=do_temp;  
    [Hc,opt_H,alpha,R,T,opt_rate,do]=opt_matrix(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode); 
     %corv1_60(itr1)=sum(sum(opt_rate));    %%revise!!!cov1
     %corv1_40(itr1)=sum(sum(opt_rate));
     end
    [p,alpha,opt_rate2]=opt_p_gp(user,match_c,do,Hc,Pmax,channel);
    end
    sum_rate(itr2)=sum(sum(opt_rate2));
    swap(itr2)=0;
    %% swaping
    flag=0;
for i=1:channel
    for j=i+1:channel
        for m=1:max_c
            for n=1:max_c 
                    %[swap_or_not,match_c]=swap_many(user,match_c,p_,Hc,do,opt_H,[i,m],[j,n],channel);
                    if mode==1
                    swap_or_not=0;
                    %[swap_or_not,match_c,p,Hc]=swap_many_3(user,match_c,p,Hc,do,opt_H,[i,m],[j,n],channel);   
                    else
                    swap_or_not=0;
                    %[swap_or_not,match_c,p,Hc]=swap_many_3(user,match_c,p,Hc,do,opt_H,[i,m],[j,n],channel);   
                    %[swap_or_not,match_c,p,Hc,alpha,do,weak1,weak2,randomornot]=swap_many_2(user,match_c,p,Hc,do,opt_H,[i,m],[j,n],channel,alpha);   
                    end
                   
                    if(swap_or_not)==1                     
                        randomswap(itr2)=randomswap(itr2)+1;
%                         str=num2str(match_c);  %%if swap, display current match
%                         disp(str);
                        total_swap=total_swap+1;
                        swap(itr2)=swap(itr2)+1;
                        flag=1;                   
                    else
                        swap(itr2)=swap(itr2)+0;
                    end
                    num_total=num_total+1;
            end
        end
    end
end

else
    sum_rate(itr2)=sum_rate(itr2-1);
    optvalue=max(sum_rate);
end
end