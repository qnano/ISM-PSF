function OTF_im=imaging_scanning_pinhole_OTF_swipfactor(OTF_ex,OTF_em,OTF_pinhole,parameters,swipfactor)

    
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

        parfor j=1:Npupil

            %em_shift=zeros(size(OTF_em));
            %ex_shift=zeros(size(OTF_ex));
            nx=i-p4; %Shift units
            ny=j-p4;
            
            ex = imtranslate(OTF_ex,[nx*(swipfactor-1.0)/swipfactor ny*(swipfactor-1.0)/swipfactor]);
            em = imtranslate(OTF_em,[-nx/swipfactor -ny/swipfactor]);
            %em=fliplr(flipud(em));
           
            
            OTF_im(i,j)=sum(sum(OTF_pinhole.*ex.*em));

        end
    end
        
        
end
    

        
    

