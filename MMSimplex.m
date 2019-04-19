% �����η�������Թ滮���� 
% max c'x
% s.t.
% Ax=b
% x>=0
% ���� A\in R^{m\times n},c,x'\in R^n,b\in R^m,b\geq 0
% ����Ҫ����˹�����ʱ���������׶η������������
% ����
% xstarΪ���Ž�;fxstarΪ����ֵ��iterΪ��������; A0Ϊ���յ����α�
% IBΪ���Ż������±�;INΪ�ǻ������±ꣻSAΪ�ɳڱ����±�
% xSAΪ�ɳ���ȡֵ;sigmaΪ���ռ�����(�ɲ���ʾ)
% A0Ϊ���յ����Ա�(��һ��Ϊb);
% InverseOfB=inv(B)Ϊ������B�������(���������ȷ���)��������Ҫѡ�� 
% By Gongnong Li 2013
function[xstar,fxstar,iter,A0,IB,IN,SA,xSA,InverseOfB,exitflag]=MMSimplex(A,b,c) 
A0=A;
[m,n]=size(A0);
E=eye(m);
IB=zeros(1,m);
SA1=zeros(1,n);
IR1=zeros(1,m);
IR=1:m;
k=0;
%���ԭ����(��׼��ʽ)ϵ���������Ƿ���E(:,i)
tic;
for i=1:m
    for j=1:n
        if A0(:,j)==E(:,i) 
            IB(i)=j;IR1(i)=i;SA1(i)=j; 
        elseif A(:,j)==(-E(:,i))
            SA1(i)=j;
        end
    end
end
s1=find(SA1~=0);
if length(s1)~=0
    for i=1:length(s1)
        SA(i)=SA1(s1(i));
    end
else
    SA=[];
end
IR=find(IR~=IR1);s=find(IR~=0);
for p=1:length(s)
    A0(:,n+p)=E(:,IR(p)); IB(IR(p))=n+p; IR(p)=n+p;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IB��¼��ԭ����ϵ�������ж��ٸ�E(:,i)����m-length(s)��,��Ӧ��x(i)Ϊ��ʼ
% ��������IR���¼��ԭ����ϵ������ȱ�ٵ�E(:,i)�±꣬��i����Щ����Ҫͨ���� 
% ����������Ĺ�length(s)���˹���������Щ����Ҳ�ǳ�ʼ��������SA��¼�ɳ�
 % ����(ʣ�����)�±�IB��¼�˻��������±꣬��IR��¼���˹��������±�(���� 
% length(s)���˹�����)���˳�ʱ����A0����һ����λ�Ӿ��󣬿��ܺ����˹��� 
% ���������˹������������������һ�׶����⣬����õ�ԭ�����һ����ʼ�� 
% ���н� 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A0=[b,A0];flag=0;
while (length(IR)~=0)&(flag==0) %��������˹���������Ҫ����һ�׶�����
    c0=zeros(n+length(s),1); c0(IR)=-ones(length(s),1); %��һ�׶ε���ؾ��������
    N=1:n+length(s); N(IB)=[]; IN=N; IN(find(IN==0))=[];
    %IB��¼�����н���±꣬IN��¼�ǻ����н���±�
    x(IN)=zeros(1,length(IN)); 
    x(IB)=A0(:,1)'; 
    cB=c0(IB);
    %��һ�׶εĳ�ʼ�����н⼰���ֵϵ��
    sigma=c0'-cB'*A0(:,2:n+length(s)+1); %����������һ��������
    t=length(find(sigma>0)); %�������������t��������ļ�����
    while t~=0
        [sigmaJ,jj]=max(sigma);
        %�����jj��sigma�о���ֵ����������У���A0�еĵ�jj+1��(A0�е�һ��Ϊb)����Ӧ�ķǻ�
	    %����x(jj)Ϊ�����������sigmaJ������Ӧ�ļ�����
        tt=find(A0(:,jj+1)>0); kk=length(tt);
        %�������ϵ������A0�е�jj+1��Ԫ���Ƿ��д������Ԫ��
        if kk==0
            disp('ԭ����Ϊ�޽��'); 
            %��A0�ĵ�jj+1��Ԫ��ȫ��С�ڻ������
            xstar=[];fxstar=[];A0=[];IB=[];iter=k;
            flag=1;
        else
            theta=zeros(1,kk);
            for i=1:kk
                theta(i)=A0(tt(i),1)/A0(tt(i),jj+1);
            end
            [thetaI,ii]=min(theta); Temp=tt(ii);
            %��ֵ��С��thetaֵ��ѡ����������TempΪ���������±ꡣ��ʱA0(Temp,jj+1)Ϊ��ת��Ԫ
            for i=1:m
                if i~=Temp
                    A0(i,:)=A0(i,:)-(A0(Temp,:)/A0(Temp,jj+1))*A0(i,jj+1);
                else
                    A0(Temp,:)=A0(Temp,:)/A0(Temp,jj+1);
                end
            end
            TT=IB(Temp); IB(Temp)=jj;
            for i=1:length(IR)
                if IR(i)==TT
                    IR(i)=0;
                end
            end
            d=find(IR==0);IR(d)=[]; %�����¼�����˹������ı仯
            IN(jj)=TT; x(IB)=A0(:,1)'; IN(find(IN==0))=[]; x(IN)=zeros(1,length(IN)); cB=c0(IB); %�µĻ����н⼰��ֵϵ��
            sigma=c0'-cB'*A0(:,2:n+length(s)+1); t=length(find(sigma>0));
            %�ٴμ�����������������������t��������ļ�����
        end
        k=k+1;
    end
    if sum(x(IR))~=0
        disp('ԭ�����޽�');%��ʱû�м�����С���㣬����һ�׶������Ž⣬�Ӷ�ԭ�����޽�
        xstar=[];fxstar=[];A0=[];IB=[];iter=k;
        flag=2;exitflag=flag;
    else
        x=x(1:n);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ�׶����������ϣ��õ�ԭ�����һ�������н� 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag==1)|(flag==2)
    return
else
    IB;N=1:n; N(IB)=[]; IN=N; IN(find(IN==0))=[];x(IN)=zeros(1,length(IN)); cB=c(IB); A0=A0(:,1:n+1); %�ص�ԭ������йؾ��������
    sigma=c'-cB'*A0(:,2:n+1); t=length(find(sigma>0));
    %����ԭ����ļ��������������������t��������ļ�����
    while (t~=0)&(flag==0)
        [sigmaJ,jj]=max(sigma);
        %jj��sigma�о���ֵ����������У���A0�еĵ�jj+1��(A0�е�һ��Ϊb)�����ж�Ӧ %�ķǻ�����x(jj)Ϊ�����������sigmaJ������Ӧ�ļ�����
        tt=find(A0(:,jj+1)>0);kk=length(tt);
        %�������ϵ������A0�е�j+1��Ԫ���Ƿ��д������Ԫ��
        if kk==0
            disp('ԭ����Ϊ�޽��'); 
            xstar=[];fxstar=[];A0=[];IB=[];
            iter=k; flag=1;
        else
            theta=zeros(1,kk);
            for i=1:kk
                theta(i)=A0(tt(i),1)/A0(tt(i),jj+1);
            end
            [thetaI,ii]=min(theta); Temp=tt(ii);
            %��ֵ��С��thetaֵ��ѡ�񻻳���������ʱA0(Temp,jj+1)Ϊ��ת��Ԫ
            for i=1:m
                if i~=Temp
                    A0(i,:)=A0(i,:)-(A0(Temp,:)/A0(Temp,jj+1))*A0(i,jj+1);
                else
                    A0(Temp,:)=A0(Temp,:)/A0(Temp,jj+1);
                end
            end
            TT=IB(Temp);IB(Temp)=jj;IN(jj)=TT; x(IB)=A0(:,1)';
            N=1:n;N(IB)=[];IN=N; IN(find(IN==0))=[];x(IN)=zeros(1,length(IN));
            cB=c(IB);
            %�µĻ����н⼰���ֵϵ��
            sigma=c'-cB'*A0(:,2:n+1);
            t=length(find(sigma>0)); %�ٴμ���������������������t��������
        end
        k=k+1;
    end
end
if flag==1
    xstar=[];fxstar=[];A0=[],IB=[];iter=k;
    disp('ԭ����Ϊ�޽��');exitflag=flag;
else
    xstar=zeros(1,n);xstar(IB)=A0(:,1)';fxstar=xstar(IB)*c(IB);iter=k;
    B=A(:,IB);InverseOfB=inv(B);xSA=x(SA);
    exitflag=flag;
end
toc;





    