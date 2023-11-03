function OTF_im=imaging_scanning_OTF(OTF_ex,OTF_em,parameters)

    
    Npupil = parameters.Npupil;

    PupilSize = 4.0;
    
    DxyPupil = 2*PupilSize/Npupil;
    XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
    [YPupil,XPupil] = meshgrid(XYPupil,XYPupil);
    
    DxyPupil = PupilSize/Npupil;
    XYPupil = -PupilSize*0.5+DxyPupil/2:DxyPupil:PupilSize*0.5;
    [YPupilm,XPupilm] = meshgrid(XYPupil,XYPupil);
    
    ex = scatteredInterpolant(XPupil(:),YPupil(:),OTF_ex(:));
    
    OTF_ex_qm=ex(XPupilm(:),YPupilm(:));
    
    em = scatteredInterpolant(XPupil(:),YPupil(:),OTF_em(:));
    
    OTF_em_qm=em(XPupilm(:),YPupilm(:));
    
    %OTF_ex_qm = interp2(XPupil,YPupil,OTF_ex,XPupilm,YPupilm);
    %OTF_em_qm = interp2(XPupil,YPupil,OTF_em,XPupilm,YPupilm);
    
    OTF_im= OTF_em_qm.*conj(OTF_ex_qm);

        
    

end