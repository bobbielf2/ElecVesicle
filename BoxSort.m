function [z,zi]=BoxSort(s,t,d,si)

m=length(s)/2;
n=length(t)/2;

if isequal(size(d),[1 1]);
    dmax=d;
else
    dmax=max(d);
end
if nargin==4
    sil=zeros(sum(si),1);
    ls=0;
    for k=1:length(si)
        sil(ls+1:ls+si(k))=k*ones(si(k),1);
        ls=ls+si(k);
    end
end
scell=cell(m,1);
ind=(1:m+n).';
v=[ind,floor([s(1:m),s(m+1:2*m);t(1:n),t(n+1:2*n)]/dmax)];
[~,I]=sort(v(:,2));
v=v(I,:);
cut=[0;ind(v(1:end-1,2)<v(2:end,2));m+n];
lc=length(cut);
for k=1:lc-1
    if k~=1
        if v(cut(k+1),2)-v(cut(k),2)==1
            ls=cut(k-1)+1;
        else
            ls=cut(k)+1;
        end
    else
        ls=1;
    end
    if k~=lc-1
        if v(cut(k+2),2)-v(cut(k+1),2)==1
            rs=cut(k+2);
        else
            rs=cut(k+1);
        end
    else
        rs=cut(lc);
    end
    vp=v(ls:rs,:);
    [~,I]=sort(vp(:,3));
    vp=vp(I,:);
    indp=(1:size(vp,1))';
    cutp=[0;indp(vp(1:end-1,3)<vp(2:end,3));size(vp,1)];
    lcp=length(cutp);
    for j=1:lcp-1
        if j~=1
            if vp(cutp(j+1),3)-vp(cutp(j),3)==1
                lsp=cutp(j-1)+1;
            else
                lsp=cutp(j)+1;
            end
        else
            lsp=1;
        end
        if j~=lcp-1
            if vp(cutp(j+2),3)-vp(cutp(j+1),3)==1
                rsp=cutp(j+2);
            else
                rsp=cutp(j+1);
            end
        else
            rsp=cutp(lcp);
        end
        vpp=vp(lsp:rsp,:);
        %[v(cut(k+1),2),vp(cutp(j+1),3)];
        
        sppo=vpp(vpp(:,1)<=m,1);
        sppi=vpp(:,1)<=m&vpp(:,2)==v(cut(k+1),2)&vpp(:,3)==vp(cutp(j+1),3);
        spp=vpp(sppi,1);
        if ~isempty(spp)
            tppi=vpp(:,1)>m;
            tpp=vpp(tppi,1);
            dpp=sqrt((t(tpp-m)*ones(1,length(spp))-ones(length(tpp),1)*s(spp)').^2+(t(tpp-m+n)*ones(1,length(spp))-ones(length(tpp),1)*s(spp+m)').^2);
            dpp1=dpp>1e-14;
            if nargin==4
                dppo=sqrt((t(tpp-m)*ones(1,length(sppo))-ones(length(tpp),1)*s(sppo)').^2+(t(tpp-m+n)*ones(1,length(sppo))-ones(length(tpp),1)*s(sppo+m)').^2);
                dppo=dppo>1e-14;
                [msi,mti]=max(~dppo);
                if ~isempty(mti)
                    mti(~dppo(1,:)')=0;
                    mti=mti(mti~=1);
                    mti(mti==0)=1;
                end
                tppn=zeros(length(tpp),1);
                tppn(mti)=sil(sppo(msi));
                nins=tppn*ones(1,length(spp))~=ones(length(tpp),1)*sil(spp)';
                dpp1=dpp1&nins;
            end
            if isequal(d,dmax)
                dpp=dpp<dmax&dpp1;
            else
                dpp=dpp<ones(length(tpp),1)*d(spp)'&dpp1;
            end
            for i=1:size(spp,1)
                scell{spp(i)}=tpp(dpp(:,i))-m;
                
            end
        end
        
    end
end
z=cell2mat(scell);
zi=cellfun('length',scell);