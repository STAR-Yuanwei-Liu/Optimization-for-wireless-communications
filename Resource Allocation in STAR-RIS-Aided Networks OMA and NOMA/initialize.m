function[do_ini,tr,h,h_ini,Hc_ini2,a,b]=initialize(user,channel,M,match_c,max_c,opt_h,h_origin,mode)
%user=6;
if mode==1
M=M/2;
end
% clear all;
% clc;
% M=10;
% user=6;
% channel=4;
%match_u=[1,1,2,2,3,3]; %current channels that users are matched with
%match_c=[1,6;2,5;2,4;3,5]; %user that channel matches  (initialized manually)
%match_c=[1,2;3,4;5,6;1,6];
%match_c=[4,5;1,6;3,2;2,4];
% max_c=2; %maximum users a channel contains
tr=zeros(channel,max_c);  % T or R user T=0 R=1

%%whether T or R according to match state
for i=1:channel
    for j=1:max_c
        if match_c(i,j)<=user/2 %%need to revise
            tr(i,j)=1;
        else
            tr(i,j)=0;
        end
    end
end
% opt_H=opt_h;  %user-channel matrix after optimizing
Hc=zeros(channel,max_c); %equal channel after matching
%h=zeros(M,channel,max_c); %original channel after matching
for i=1:channel
    for j=1:max_c
        user_index=match_c(i,j);
        Hc(i,j)=opt_h(user_index,i);
        h(:,i,j)=h_origin(:,user_index,i);
    end
end
%decoding order according to equal channel 
do=zeros(channel,max_c); %decoding order

%save('initial','do','tr','h_origin','h');


%% R user initialize 
m=max(max(opt_h(1:user/2,:)));
%m=min(min(opt_h(1:user/2,:)));
[roll1,cell1]=find(opt_h==m);
angle1=-angle(h_origin(:,roll1,cell1));
vector1=sqrt(0.5)*(cos(angle1)+1i*sin(angle1));

% vector1=randn(M,1)+1i*randn(M,1);
% vector2=randn(M,1)+1i*randn(M,1);
% for m=1:M
% vector1(m)=sqrt(0.5)*vector1(m)/norm(vector1(m));
% vector2(m)=sqrt(0.5)*vector2(m)/norm(vector2(m));
% end


for i=1:user/2 
for j=1:channel
for m=1:M
%h_ini(i,j)=h_origin(:,i,j).'*vector1; %注意为普通转置
h_ini_1(m,i,j)=h_origin(m,i,j)*vector1(m);
end
end
end
for i=1:user/2 
for j=1:channel
h_ab(i,j)=sum(h_ini_1(:,i,j));
h_ini(i,j)=norm(sum(h_ini_1(:,i,j)))^2;
end
end
%h_ini=(abs(h_ini))^2;



%%T user initialize
n=max(max(opt_h(user/2+1:end,:)));
%n=min(min(opt_h(user/2+1:end,:)));
[roll2,cell2]=find(opt_h==n);
angle2=-angle(h_origin(:,roll2,cell2));
vector2=sqrt(0.5)*(cos(angle2)+1i*sin(angle2));

for i=user/2+1:user 
for j=1:channel
for m=1:M
%h_ini(i,j)=h_origin(:,i,j).'*vector1; %注意为普通转置
h_ini_1(m,i,j)=h_origin(m,i,j)*vector2(m);
end
end
end
for i=user/2+1:user 
for j=1:channel
h_ab(i,j)=sum(h_ini_1(:,i,j));
h_ini(i,j)=norm(sum(h_ini_1(:,i,j)))^2;
end
end


for i=1:channel
    for j=1:max_c
        user_index=match_c(i,j);
        Hc_ini2(i,j)=h_ini(user_index,i);
        a(i,j)=real(h_ab(user_index,i));
        b(i,j)=imag(h_ab(user_index,i));
    end
end



for i=1:channel
    for j=1:max_c
       if Hc_ini2(i,j)>=Hc_ini2(i,mod(j,2)+1)
           do(i,j)=2;
       else
           do(i,j)=1;
       end
    end
end
do_ini=do;