function OTF_im=confocal_pinhole_OTF(OTF_ex,OTF_em,OTF_pinhole,parameters)

    
    Npupil = parameters.Npupil;

    PupilSize = 4.0;
    
    %DxyPupil = 2*PupilSize/Npupil;
    %XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
    %[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

    OTF_ex=conj(OTF_ex);
    %OTF_em=fliplr(flipud(OTF_em));
    OTF_im=zeros(Npupil);
    p4=(Npupil-1)/2;
    for i=1:Npupil

        for j=1:Npupil

            %em_shift=zeros(size(OTF_em));
            %ex_shift=zeros(size(OTF_ex));
            nx=i-p4; %Shift units
            ny=j-p4;
            
            ex = OTF_ex;
            em = imtranslate(OTF_em,[-nx -ny]);
            %em=fliplr(flipud(em));
           
            
            OTF_im(i,j)=sum(sum(OTF_pinhole.*ex.*em));

        end
    end
        
        
end
    

        
    

