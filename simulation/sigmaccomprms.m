clear all
mat=load('C:/code/ExamplePhaseRetrieval/data/rmsdata/18rms-18rms_crlb.mat');
Icstore1818=mat.CRLBarr;

mat=load('C:/code/ExamplePhaseRetrieval/data/rmsdata/36rms-18rms_crlb.mat');
Icstore3618=mat.CRLBarr;

mat=load('C:/code/ExamplePhaseRetrieval/data/rmsdata/36rms-36rms_crlb.mat');
Icstore3636=mat.CRLBarr;

if true
    z=linspace(-1,1,11);

    m=mean(Icstore1818,3);
    s=std(Icstore1818,0,3);
    errorbar(z,m(1,:),s(1,:));
    hold on
    m=mean(Icstore3618,3);
    s=std(Icstore3618,0,3);
    errorbar(z,m(1,:),s(1,:));
    hold on
    m=mean(Icstore3636,3);
    s=std(Icstore3636,0,3);
    errorbar(z,m(1,:),s(1,:));
    legend({'18rms data+18rms model','36rms data+18rms model','36rms data+36rms model'},'Location','northwest')
    xlabel("z-offset (um)");
    ylabel("\sigma_{c}(nm)");
end


if false
    bgarr=linspace(3,50,10);
    m=mean(Icstore1818,3);
    s=std(Icstore1818,0,3);
    errorbar(bgarr,m(:,6),s(:,6));
    hold on
    m=mean(Icstore3618,3);
    s=std(Icstore3618,0,3);
    errorbar(bgarr,m(:,6),s(:,6));
    hold on
    m=mean(Icstore3636,3);
    s=std(Icstore3636,0,3);
    errorbar(bgarr,m(:,6),s(:,6));
    legend({'18rms data+18rms model','36rms data+18rms model','36rms data+36rms model'},'Location','northwest')
    xlabel("bg");
    ylabel("\sigma_{c}(nm)");
    %set(gca, 'YScale', 'log')
end