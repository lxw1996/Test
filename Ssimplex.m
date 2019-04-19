% Ssimplex.m ���õ����η���MATLAB����
% max c'x
% a.t.
% ax=b
% x>=0
% ����A\in R^{m\times n},c,x'\im R^m,b\geq 0
% �Ҿ���A����һ����λ�Ӿ��󣬲���Ҫ�����˹�����
% By Gongnong Li 2013
function [xstar,fxstar,iter]=Ssimplex(A,b,c)
[m,n]=size(A);E=eye(m);IB=zeros(1,m);k=0;
for i = 1:m
    for j=1:n
        if A(:,j)==E(:,i)
            IB(i)=j;SA(i)=j;%IB��¼�������±꣬SA��¼�ɳڱ����±�
        elseif A(:,j)==(-E(:,i))
            SA(i)=j;%SAҲ��¼ʣ��������ɳڱ������±�
        end
    end
end
A0=[b,A];N=1:n;N(IB)=[];IN=N;x(IB)=A0(:,1)';
x(IN)=zeros(1,length(IN));cB=c(IB);
%INΪ�ǻ������±�
sigma=c'-cB'*A0(:,2:n+1);t=length(find(sigma>0));
%����ԭ����ļ��������������������t��������ļ�����
while t~=0
    [sigmaJ,jj]=max(sigma);
    %�����jj��sigma��ֵ����������У���A0�еĵ�jj+1�У�A0�еĵ�һ��Ϊb�������ж�Ӧ�ķǻ�����x(jj)Ϊ�����������sigmaJ������Ӧ�ļ�����
    tt=find(A0(:,jj+1)>0);kk=length(tt);
    %�������ϵ������A0�е�jj+1��Ԫ���Ƿ��д������Ԫ��
    if kk==0
        disp('ԭ����Ϊ�޽��')
    else
        theta =zeros(1,kk);
        for i=1:kk
            theta(i)=A0(tt(i),1)/A0(tt(i),jj+1);
        end
        [thetaI,ii]=min(theta);Temp=tt(ii);
        %��ֵ��С��thetaֵ�� ѡ�񻻳������� ��ʱA0(Temp,jj+1)Ϊ��ת��Ԫ
        for i=1:m
            if i~=Temp
                A0(i,:)=A0(i,:)-(A0(Temp,:)/A0(Temp,jj+1))*A0(i,jj+1);
            else
                A0(Temp,:)=A0(Temp,:)/A0(Temp,jj+1);
            end
        end
        TT=IB(Temp);IB(Temp)=jj;IN(jj)=TT;x(IB)=A(:,1)';
        N=1:n;N(IB)=[];IN=N;x(IN)=zeros(1,length(IN));cB=c(IB);
        %�µĻ����н⼰���ֵϵ��
        sigma=c'-cB'*A0(:,2:n+1);t=length(find(sigma>0));
        %�ٴμ�����������������������t��������ļ�����
    end
    k=k+1;
end
IB
IN
B=A(:,IB);
InverseOfB=inv(B)
%���ǻ�����B����������������ȷ����������������ȷ���������ע�͵�
xstar=x;fxstar=x(IB)*c(IB);iter=k;

