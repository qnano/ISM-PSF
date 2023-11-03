n=30;
I=1e4;
rate=0.9;
Iall=zeros(n);
for i=1:n
    Iall(i)=I*rate^i;
end
Iall=flip(Iall);
bgarr=linspace(5,50,10);
matfile=load('C:/code/ExamplePhaseRetrieval/data/simdata/vector_fitter_amp.mat');
PSF=matfile.PSF;
mat=load('C:/code/ExamplePhaseRetrieval/data/simdata/chiarraybg10fine_amp.mat');
zs=21;
tn=301;
modelfact=zeros(10,n,21);
for bi=1:10
    for ni=1:20
        for k=1:21
            cm=PSF.*Iall(ni)+bgarr(bi);
            modelfact(bi,ni,k)=sum(1./cm,'all');
        end
    end
end


chizstackbg=mat.chiarraybg;
Icstore=zeros(10,zs,tn);
crlb60=zeros(10,zs,tn);
m=1;
expv=31*31;
for bgind=1:10
    for i =1:zs
        cspline_zcintensitytmp=zeros(tn);
        for j =1:tn
            find=0;
            for k =2:n
                if (chizstackbg(bgind,k,i,j)>(expv+m*sqrt(2*expv+(modelfact(bgind,k,i)))) && chizstackbg(bgind,k-1,i,j)>(m*sqrt(modelfact(bgind,k-1,i)+2*expv)) && find==0)
                    Icstore(bgind,i,j)=Iall(k);
                end
            end
        end
    end
end