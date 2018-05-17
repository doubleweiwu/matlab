# matlab
doublewei的matlab
function Y_out=FOS_convert (total_no,sampling_frequency, frequeny_areas, candidate_no,y_in)
% clear
% clc
  doration=1/sampling_frequency;
  t=doration:doration:total_no*doration;
  N=total_no;
  y=y_in;
  sum=0;
  %初始函数
  M=candidate_no; 
    for i=1:M
        w(i)=i*100;
        p_candidate(2*i-1,:)=sin(w(i)*t);
        p_candidate(2*i,:)=cos(w(i)*t);
    end
    
D(1,1)=1;
C(1)=mean(y);
sumgm=0;%统计GM初值
p(1,1:N)=1;%derectly 
%%%%%%%%%candidate functional expansion
for m_screening=1:M% 一共是100对，进行每一对的计算
    p(2,:)=p_candidate(2*m_screening-1,:);
    p(3,:)=p_candidate(2*m_screening,:);
    for m=2:3
        D(m,1)=mean(p(m,:)) ;
    end
    for m=2:3
        for r=1:m-1
            afa(m,r)=D(m,r)/D(r,r);
            for ii=1:r+1-1
                sum=sum+afa(r+1,ii)*D(m,ii);
            end
            D(m,r+1)=p(m,:)*p(r+1,:)'/N-sum;
            sum=0;
        end
    end
    for m=2:3
         for ii=1:m-1
                sum=sum+afa(m,ii)*C(ii);
         end
        C(m)=y*p(m,:)'/N-sum;
        sum=0;
    end
    for m=1:3
        g(m)=C(m)/D(m,m);  
    end
    Q(m_screening)=g(2)^2*D(2,2)+g(3)^2*D(3,3);
end
plot(Q)
[x_maxnum,y_maxnum]=find(Q==max(Q)); %find the max Q 
%%%%选出来一次 重新算
p(2,:)= p_candidate(2*y_maxnum(1)-1,:);% PROTECT FOR DOUBLE VAR
p(3,:)= p_candidate(2*y_maxnum(1),:);% PROTECT FOR DOUBLE VAR
    for m=2:3
        D(m,1)=mean(p(m,:)) ;
    end
    for m=2:3;
        for r=1:m-1
            afa(m,r)=D(m,r)/D(r,r);
            for ii=1:r+1-1
                sum=sum+afa(r+1,ii)*D(m,ii);
            end
            D(m,r+1)=p(m,:)*p(r+1,:)'/N-sum;
            sum=0;
        end
    end
    for m=2:3
         for ii=1:m-1
                sum=sum+afa(m,ii)*C(ii);
         end
        C(m)=y*p(m,:)'/N-sum;
        sum=0;
    end
    for m=1:3;
        g(m)=C(m)/D(m,m);  
    end    
for m=1:3
    sumgm=sumgm+g(m)^2*D(m,m);
end
mse=y*y'/N-sumgm;
%%%%%%%%%%%选出来的计算完毕    
threthold=0.01*var(y);
if mse>=threthold %大于阈值 循环开始 terminate until mse achieve the threshold
    order=2;%begin from the second order 
%     while(mse>=0.08) 
      for iii=1:10
        fp=[p_candidate(1:2*(y_maxnum-1),:); p_candidate( 2*(y_maxnum)+1: 2*(M-order+2) ,:)];
        p_candidate(1:2*(M-order+1),:)=fp;
        for m_screening=1:M-order% 一共是100对，进行每一对的计算
            p(2*order,:)=p_candidate(2*m_screening-1,:);
            p(2*order+1,:)=p_candidate(2*m_screening,:);
            for m=2*order:2*order+1;%caculate D(m,0)
                D(m,1)=mean(p(m,:));
            end
            for m=2*order:2*order+1;%caculate D(m,r)
                for r=1:m-1
                    afa(m,r)=D(m,r)/D(r,r);
                    for ii=1:r+1-1
                        sum=sum+afa(r+1,ii)*D(m,ii);
                    end
                    D(m,r+1)=p(m,:)*p(r+1,:)'/N-sum;
                    sum=0;
                end
            end
            for m=2*order:2*order+1;%caculate C(m)
                 for ii=1:m-1
                        sum=sum+afa(m,ii)*C(ii);
                 end
                C(m)=y*p(m,:)'/N-sum;
                sum=0;
            end
            for m=2*order:2*order+1;
                g(m)=C(m)/D(m,m);  %caculate G(m)
            end
            Q(m_screening)=g(2*order)^2*D(2*order,2*order)+g(2*order+1)^2*D(2*order+1,2*order+1);%每次order增加 Q会减少一个值
        end
%         figure(1)
%         plot(Q)
        [x_maxnum,y_maxnum]=find(Q==max(Q(1:m_screening))); %find the max Q y_maxnum被更替了
        %%%%选出来一个 重新算
        p(2*order,:)= p_candidate(2*y_maxnum(1)-1,:);% PROTECT FOR DOUBLE VAR
        p(2*order+1,:)= p_candidate(2*y_maxnum(1),:);% PROTECT FOR DOUBLE VAR
        for m=2*order:2*order+1
            D(m,1)=mean(p(m,:)) ;
        end
        for m=2*order:2*order+1;
            for r=1:m-1
                afa(m,r)=D(m,r)/D(r,r);
                for ii=1:r+1-1
                    sum=sum+afa(r+1,ii)*D(m,ii);
                end
                D(m,r+1)=p(m,:)*p(r+1,:)'/N-sum;
                sum=0;
            end
        end
        for m=2*order:2*order+1
             for ii=1:m-1
                    sum=sum+afa(m,ii)*C(ii);
             end
            C(m)=y*p(m,:)'/N-sum;
            sum=0;
        end
        for m=2*order:2*order+1;
            g(m)=C(m)/D(m,m);  
        end
        for m=2*order:2*order+1;
            sumgm=sumgm+g(m)^2*D(m,m);%在原有基础上再次减去新加入的函数
        end
        mse=y*y'/N-sumgm;
        order=order+1;
    end% for while
end % for if

suma=0;
sumv=0;

for m=1:size(p,1);%caculate the am
    for i=m:size(p,1)
        if(m==i)
            v(m)=1;
        else
            for r=m:i-1
              sumv=sumv+afa(i,r)*v(r);
            end
            v(i)=-sumv;
            sumv=0;
        end       
        suma=suma+g(i)*v(i); 
    end    
    a(m)=suma;
    suma=0;
end
yout=0;
for m=1:size(p,1);
   yout=yout+a(m)*p(m,:);
end
figure (2)
plot(y)
hold on 
plot(yout,'-r*')
legend('real','FOS')



%%%%%%%%%%%%%%%
% for m=1:M;
%     for i=m:M
%         if(m==i)
%             v(m)=1;
%         else
%             for r=m:i-1
%               sumv=sumv+afa(i,r)*v(r);
%             end
%             v(i)=-sumv;
%             sumv=0;
%         end       
%         suma=suma+g(i)*v(i); 
%     end
%     a(m)=suma;
% end
% yout=0;
% for m=1:M;
%    yout=yout+a(m)*p(m,:)
% end
% figure
% plot(yout)
% hold on
% plot(y,'-r')
% 
% 
% mse=(y-yout)*(y-yout)'/N;
% 
% sumgm=0;
% for m=1:M
%     sumgm=sumgm+g(m)^2*D(m,m);
% end
% 
% mse1=y*y'/N-sumgm;
% 
% MSE=g(M)^2*D(M,M);
% 
% 
% fs=100;%设定采样频率
% N=128;
% n=0:N-1;
% y=fft(y,N);%进行fft变换
% mag=abs(y);%求幅值
% f=(0:length(y)-1)'*fs/length(y);%进行对应的频率转换
% figure(2);
% plot(f,mag);%做频谱图
% axis([0,100,0,80]);
% xlabel('频率(Hz)');
% ylabel('幅值');
% title('正弦信号y=2*pi*10t幅频谱图N=128');
% grid;
       
fs=20000;%设定采样频率
N=1024;
n=0:N-1;
y=fft(y,N);%进行fft变换
mag=abs(y);%求幅值
f=(0:length(y)-1)'*fs/length(y);%进行对应的频率转换
figure(3);
plot(f,mag);%做频谱图
% axis([0,100,0,80]);
xlabel('频率(Hz)');
ylabel('幅值');
title('正弦信号y=2*pi*10t幅频谱图N=128');
grid;
% for m=2:M
%     D(m,1)=mean(p(m,:)) ;
% end
% 
% for m=2:M;
%     for r=1:m-1
%         afa(m,r)=D(m,r)/D(r,r);
%         for ii=1:r+1-1
%             sum=sum+afa(r+1,ii)*D(m,ii);
%         end
%         D(m,r+1)=p(m,:)*p(r+1,:)'/N-sum;
%         sum=0;
%     end
% end
% 
% C(1)=mean(y);
% 
% for m=2:M
%      for ii=1:m-1
%             sum=sum+afa(m,ii)*C(ii);
%      end
%     C(m)=y*p(m,:)'/N-sum;
%     sum=0;
% end
% 
% for m=1:M;
%     g(m)=C(m)/D(m,m);  
% end
% suma=0;
% sumv=0;
% for m=1:M;
%     for i=m:M
%         if(m==i)
%             v(m)=1;
%         else
%             for r=m:i-1
%               sumv=sumv+afa(i,r)*v(r);
%             end
%             v(i)=-sumv;
%             sumv=0;
%         end       
%         suma=suma+g(i)*v(i); 
%     end
%     a(m)=suma;
% end
% yout=0;
% for m=1:M;
%    yout=yout+a(m)*p(m,:)
% end
% figure
% plot(yout)
% hold on
% plot(y,'-r')
% 
% 
% mse=(y-yout)*(y-yout)'/N;
% 
% sumgm=0;
% for m=1:M
%     sumgm=sumgm+g(m)^2*D(m,m);
% end
% 
% mse1=y*y'/N-sumgm;
% 
% MSE=g(M)^2*D(M,M);
% 
% 
% fs=100;%设定采样频率
% N=128;
% n=0:N-1;
% y=fft(y,N);%进行fft变换
% mag=abs(y);%求幅值
% f=(0:length(y)-1)'*fs/length(y);%进行对应的频率转换
% figure(2);
% plot(f,mag);%做频谱图
% axis([0,100,0,80]);
% xlabel('频率(Hz)');
% ylabel('幅值');
% title('正弦信号y=2*pi*10t幅频谱图N=128');
% grid;
%        
% fs=100;%设定采样频率
% N=128;
% n=0:N-1;
% y=fft(yout,N);%进行fft变换
% mag=abs(y);%求幅值
% f=(0:length(y)-1)'*fs/length(y);%进行对应的频率转换
% figure(3);
% plot(f,mag);%做频谱图
% axis([0,100,0,80]);
% xlabel('频率(Hz)');
% ylabel('幅值');
% title('正弦信号y=2*pi*10t幅频谱图N=128');
% grid;
        
        
