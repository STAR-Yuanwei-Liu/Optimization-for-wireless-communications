%%near exhaustive
function[swap_or_not,match_c,p,Hc]=swap_many_3(user,match_c,p,Hc,do,opt_H,first,second,channel)
%swap blocking pair 
% clear all;
% clc;
%user=6;
%channel=4;
delta=1e-11;
max_c=2;
%Hc=zeros(channel,max_c);
%opt_H=zeros(user,channel);
%p=ones(channel,max_c)*0.1;
%p=[0.12,0.18,0.1,0.08,0.14,0.04];
%first=[1,1];
%second=[2,1]; %second>first
%load('channel','Q','h','rate');
%match_c=[1,6;2,5;2,4;3,5]; %user that channel matches


%do=[2,1;2,1;2,1;2,1];
%tr=zeros(channel,max_c); 

c1=first(1);
i1=first(2);
c2=second(1);
i2=second(2);
%% utility before swap
%first user
if  do(c1,i1)==2  
u_fir_u=log2(1+p(c1,i1)*Hc(c1,i1)/delta);
u_fir_c=log2(1+p(c1,i1)*Hc(c1,i1)/delta)+log2(1+p(c1,mod(i1,2)+1)*Hc(c1,mod(i1,2)+1)/(p(c1,i1)*Hc(c1,mod(i1,2)+1)+delta));
else
u_fir_u=log2(1+p(c1,i1)*Hc(c1,i1)/(p(c1,mod(i1,2)+1)*Hc(c1,i1)+delta)); 
u_fir_c= log2(1+p(c1,mod(i1,2)+1)*Hc(c1,mod(i1,2)+1)/delta)+log2(1+p(c1,i1)*Hc(c1,i1)/(p(c1,mod(i1,2)+1)*Hc(c1,i1)+delta)); 
end
%second user
if  do(c2,i2)==2  
u_sec_u=log2(1+p(c2,i2)*Hc(c2,i2)/delta);
u_sec_c=log2(1+p(c2,i2)*Hc(c2,i2)/delta)+log2(1+p(c2,mod(i2,2)+1)*Hc(c2,mod(i2,2)+1)/(p(c2,i2)*Hc(c2,mod(i2,2)+1)+delta));
else
u_sec_u=log2(1+p(c2,i2)*Hc(c2,i2)/(p(c2,mod(i2,2)+1)*Hc(c2,i2)+delta)); 
u_sec_c= log2(1+p(c2,mod(i2,2)+1)*Hc(c2,mod(i2,2)+1)/delta)+log2(1+p(c2,i2)*Hc(c2,i2)/(p(c2,mod(i2,2)+1)*Hc(c2,i2)+delta)); 
end


%% utility after swap 
p_swap=p;
p_temp=p(c1,i1);
p_swap(c1,i1)=p_swap(c2,i2);
p_swap(c2,i2)=p_temp;


Hc_temp1=opt_H(match_c(c1,i1),c2);  
Hc_temp2=opt_H(match_c(c2,i2),c1);
Hc_swap=Hc;
Hc_swap(c1,i1)=Hc_temp2;
Hc_swap(c2,i2)=Hc_temp1;
if Hc_swap(c2,i2)>=Hc_swap(c2,mod(i2,2)+1)
    do_fir=2; %docoding order of first user after swapping
else
    do_fir=1;
end
if Hc_swap(c1,i1)>=Hc_swap(c1,mod(i1,2)+1)
    do_sec=2;
else
    do_sec=1;
end

%% utility after swap  
%first user
if  do_fir==2  
u_fir_u2=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/delta);
u_sec_c2=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/delta)+log2(1+p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,mod(i2,2)+1)/(p_swap(c2,i2)*Hc_swap(c2,mod(i2,2)+1)+delta));
weak1=log2(1+p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,mod(i2,2)+1)/(p_swap(c2,i2)*Hc_swap(c2,mod(i2,2)+1)+delta));
else
u_fir_u2=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/(p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,i2)+delta));
u_sec_c2=log2(1+p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,mod(i2,2)+1)/delta)+log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/(p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,i2)+delta));
weak1=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/(p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,i2)+delta));
end
%second user
if  do_sec==2  
u_sec_u2=log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/delta);
u_fir_c2=log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/delta)+log2(1+p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,mod(i1,2)+1)/(p_swap(c1,i1)*Hc_swap(c1,mod(i1,2)+1)+delta));
weak2=log2(1+p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,mod(i1,2)+1)/(p_swap(c1,i1)*Hc_swap(c1,mod(i1,2)+1)+delta));
else
u_sec_u2= log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/(p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,i1)+delta)); 
u_fir_c2= log2(1+p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,mod(i1,2)+1)/delta)+log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/(p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,i1)+delta)); 
weak2= log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/(p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,i1)+delta)); 
end

% channel utility 原用户交换过去后的utility+原用户旧配对用户的新utility
if do_fir==2 && do_sec==2
u_fir_c2=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/delta)+log2(1+p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,mod(i1,2)+1)/(p_swap(c1,i1)*Hc_swap(c1,mod(i1,2)+1)+delta));
u_sec_c2=log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/delta)+log2(1+p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,mod(i2,2)+1)/(p_swap(c2,i2)*Hc_swap(c2,mod(i2,2)+1)+delta));
elseif do_fir==2 && do_sec==1
u_fir_c2=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/delta)+log2(1+p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,mod(i1,2)+1)/delta);
u_sec_c2=log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/(p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,i1)+delta))+log2(1+p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,mod(i2,2)+1)/(p_swap(c2,i2)*Hc_swap(c2,mod(i2,2)+1)+delta));
elseif do_fir==1 && do_sec==2
u_fir_c2=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/(p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,i2)+delta))+log2(1+p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,mod(i1,2)+1)/(p_swap(c1,i1)*Hc_swap(c1,mod(i1,2)+1)+delta));
u_sec_c2=log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/delta)+log2(1+p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,mod(i2,2)+1)/delta);
else %fir==1,sec==1
u_fir_c2=log2(1+p_swap(c2,i2)*Hc_swap(c2,i2)/(p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,i2)+delta))+ log2(1+p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,mod(i1,2)+1)/delta);
u_sec_c2=log2(1+p_swap(c1,i1)*Hc_swap(c1,i1)/(p_swap(c1,mod(i1,2)+1)*Hc_swap(c1,i1)+delta))+log2(1+p_swap(c2,mod(i2,2)+1)*Hc_swap(c2,mod(i2,2)+1)/delta);
end


%% whether to swap
%if u_fir_u2>=u_fir_u && u_fir_c2>=u_fir_c && u_sec_u2>=u_sec_u && u_sec_c2>=u_sec_c && (u_fir_u2>u_fir_u || u_fir_c2>u_fir_c || u_sec_u2>u_sec_u || u_sec_c2>u_sec_c)
if u_fir_u2>=u_fir_u && u_sec_u2>=u_sec_u  && (u_fir_u2>u_fir_u || u_sec_u2>u_sec_u)  
swap_or_not=1;

% %direct change
% Hc=Hc_swap;
% p=p_swap;
% %decoding order
% do(c2,i2)=do_fir;
% do(c1,i1)=do_sec;


%%%matching state  important!
match_temp=match_c(c1,i1);
match_c(c1,i1)=match_c(c2,i2);
match_c(c2,i2)=match_temp;

% p_temp=p(c1,i1);
% p(c1,i1)=p(c2,i2);
% p(c2,i2)=p_temp;
% alpha_temp=alpha(c1,i1);
% alpha(c1,i1)=alpha(c2,i2);
% alpha(c2,i2)=alpha_temp;
Hc=Hc_swap;
p=p_swap;
% do(c2,i2)=do_fir;
% do(c2,mod(i2,2)+1)=mod(do_fir,2)+1;
% do(c1,i1)=do_sec;
% do(c1,mod(i1,2)+1)=mod(do_sec,2)+1;
% if do_fir==2
%     alpha(c2)=0.01/Hc(c2,mod(i2,2)+1)/1e11;
% else
%     alpha(c2)=0.01/Hc(c2,i2)/1e11;
% end
% 
% if do_sec==2
%     alpha(c1)=0.01/Hc(c1,mod(i1,2)+1)/1e11;
% else
%     alpha(c1)=0.01/Hc(c1,i1)/1e11;
% end

else
    u_total=u_fir_u+u_sec_u;
    u_total_new=u_fir_u2+u_sec_u2;
    random_num= rand(1);
    prob=1/(1+exp(-0.5*((u_total_new-u_total)/u_total)));%probability
    if weak1>=0.01&&weak2>=0.01
    if random_num<prob&&((u_total_new-u_total)>0)
        swap_or_not=1; %swap
        match_temp=match_c(c1,i1);
        match_c(c1,i1)=match_c(c2,i2);
        match_c(c2,i2)=match_temp;
        Hc=Hc_swap;
        p=p_swap;
%         do(c2,i2)=do_fir;
%         do(c2,mod(i2,2)+1)=mod(do_fir,2)+1;
%         do(c1,i1)=do_sec;
%         do(c1,mod(i1,2)+1)=mod(do_sec,2)+1;
%         if do_fir==2
%             alpha(c2)=(2^0.01-1)/Hc(c2,mod(i2,2)+1)/1e11;
%         else
%             alpha(c2)=(2^0.01-1)/Hc(c2,i2)/1e11;
%         end
% 
%         if do_sec==2
%             alpha(c1)=(2^0.01-1)/Hc(c1,mod(i1,2)+1)/1e11;
%         else
%             alpha(c1)=(2^0.01-1)/Hc(c1,i1)/1e11;
%         end
    else
        swap_or_not=0;
    end 
    else
        swap_or_not=0;
    end
end

