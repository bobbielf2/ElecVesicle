function [CM,FM]=IntMat(s,si,t)

n=length(s)/2;
m=length(t)/2;

if nargout==2
    FM=MakeFM(s,si,t);
end

h=Findh(s,si);

% Determines target points within 5h
[z,zi]=BoxSort(s,t,h,si);
[z,zi]=MakeUnique(z,zi,si);

[vi,vj,cm,cs]=CloseEval(s,si,t,z,zi);
if isequal(s,t)
    [vi,vj,cm]=CorrectFar(s,si,vi,vj,cm,cs);
end
CM=sparse(vi,vj,cm,2*m,2*n);
end

function FM=MakeFM(s,si,t)
n=length(s)/2;
m=length(t)/2;
M=length(si);
FM=zeros(2*m,2*n);
ls=0;
tf.x=t(1:end/2)+1i*t(end/2+1:end);
for k=1:M
    sf.x=s(ls+1:ls+si(k))+1i*s(end/2+(ls+1:ls+si(k)));
    sf=quadr(sf);
    SLP=SLPmatrix(tf,sf,1);
    FM(:,ls+1:ls+si(k))=SLP(:,1:end/2);
    FM(:,end/2+(ls+1:ls+si(k)))=SLP(:,end/2+1:end);
    ls=ls+si(k);
end
if isequal(s,t)
    ls=0;
    for k=1:M
        for j=1:si(k)
            FM(ls+j,ls+j)=0;
            FM(end/2+ls+j,ls+j)=0;
            FM(end/2+ls+j,end/2+ls+j)=0;
            FM(ls+j,end/2+ls+j)=0;
        end
        ls=ls+si(k);
    end
end
end

function h=Findh(s,si)
M=length(si);
h=zeros(sum(si),1);
ls=0;
for k=1:M
    h(ls+1:ls+si(k))=5*ArcLength(s([ls+si(k),ls+1:ls+si(k)]),s(end/2+[ls+si(k),ls+1:ls+si(k)]),2*pi)/si(k);
    ls=ls+si(k);
end
end

function [z,zi]=MakeUnique(z,zi,si)
M=length(si);
ls=0;
rs=0;
uniz=cell(M,1);
for k=1:M
    szi=sum(zi(rs+1:rs+si(k)));
    uniz{k}=unique(z(ls+1:ls+szi));
    ls=ls+szi;
    rs=rs+si(k);
end
z=cell2mat(uniz(cellfun('length',uniz)>0));
zi=cellfun('length',uniz);
end

function [vi,vj,cm,cs]=CloseEval(s,si,t,z,zi)
n=length(s)/2;
m=length(t)/2;
M=length(si);

if isequal(s,t)
    vi=zeros(4*(zi'*si+si'*si),1);
else
    vi=zeros(4*zi'*si,1);
end
vj=vi;
cm=vi;

ls=0;
rs=0;
cs=0;
for k=1:M
    ri=z(ls+1:ls+zi(k));
    ci=(rs+1:rs+si(k))';
    sz=4*si(k)*zi(k);
    tp.x=t(ri)+1i*t(m+ri);
    sp.x=s(ci)+1i*s(n+ci);
    sp=quadr(sp);
    if isempty(tp.x)
        CSLP=zeros(0,2*length(sp.x));
        SLP=CSLP;
    else
        CSLP=StokesScloseevalF(tp.x,sp,'e');
        %CSLP=StokesScloseeval(tp.x,sp,[],'e');
        SLP=SLPmatrix(tp,sp,1);
    end
    [c,r]=meshgrid([ci;n+ci],[ri;m+ri]);
    vi(cs+1:cs+sz)=reshape(r,[sz,1]);
    vj(cs+1:cs+sz)=reshape(c,[sz,1]);
    cm(cs+1:cs+sz)=reshape(CSLP-SLP,[sz,1]);
%                 scatter(real(sp.x),imag(sp.x),'b')
%                 hold on
%                 scatter(real(tp.x),imag(tp.x),'r.')
    ls=ls+zi(k);
    rs=rs+si(k);
    cs=cs+sz;
end
end

function [vi,vj,cm]=CorrectFar(s,si,vi,vj,cm,cs)
n=length(s)/2;
M=length(si);
rs=0;
for k=1:M
    ci=(rs+1:rs+si(k))';
    sz=4*si(k)*si(k);
    st.x=s(ci)+1i*s(n+ci);
    st=quadr(st);
    SLP=SLPmatrix(st,st,1);
    for j=1:si(k)
        SLP(j,j)=0;
        SLP(si(k)+j,si(k)+j)=0;
        SLP(si(k)+j,j)=0;
        SLP(j,si(k)+j)=0;
    end
    [c,r]=meshgrid([ci;n+ci],[ci;n+ci]);
    vi(cs+1:cs+sz)=reshape(r,[sz,1]);
    vj(cs+1:cs+sz)=reshape(c,[sz,1]);
    cm(cs+1:cs+sz)=reshape(-SLP,[sz,1]);
    rs=rs+si(k);
    cs=cs+sz;
end
end