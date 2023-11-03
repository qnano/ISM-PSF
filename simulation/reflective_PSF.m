classdef reflective_PSF

properties
        n_pixels;
        n_z_slices;
        NA;
        pixelsize;
        lambda;
        zrange;
        d;
    end
    
    
    methods
        function self=reflective_PSF(n_pixels,n_z_slices,NA,pixelsize,lambda,zrange,d,aberrations)
            mat=load('C:/code/jelmer/defult_parameter.mat');
            parameters=mat.parameters;
            d=d;
            count=1;
            ma = 5;
            for i=1:ma
                mind=-i:2:i;
                for j=1:length(mind)
                    aberrations(count,:)=[i,mind(j),0,0];
                    count=count+1;
                end
            end
            aberrations(3,2)=0;
            aberrations(4,2)=-2;

            aberrations(12,3)=0;
            aberrations(13,3)=50;


            parameters.aberrations=aberrations;

            parameters.NA = NA;

            parameters.xemit = 0;
            parameters.yemit = 0;

            parameters.pixelsize=pixelsize;
            %parameters.samplingdistance = parameters.pixelsize;
            parameters.fitmodel='aberrations';
            %aberrations = [1,-1,0,0; 1,1,0.0,0; 2,0,0,0.05; 2,-2,0,0.0; 2,2,50.0,0.01; 3,-1,0.0,0.0; 3,1,0.0,0; 4,0,0.0,0.1; 3,-3,-0.0,0; 3,3,0.0,0; 4,-2,0.0,0; 4,2,0.0,0; 5,-1,0.0,0; 5,1,0.0,0; 6,0,0.0,0; 4,-4,0.0,0; 4,4,0.0,0;  5,-3,0.0,0; 5,3,0.0,0;  6,-2,0.0,0; 6,2,0.0,0; 7,1,0.0,0; 7,-1,0.0,0; 8,0,0.0,0];
            parameters.beaddiameter = 40;
            parameters.lambda = lambda;
            parameters.lambdacentral=lambda;
            parameters.lambdaspread=[lambda lambda];
            parameters.zrange = zrange;
            parameters.doetype='none';
            parameters.dipoletype='free';
            parameters.Mx = n_pixels ;
            parameters.My = n_pixels ;
            parameters.Mz = n_z_slices;

            parameters.xrange = parameters.pixelsize*parameters.Mx/2;
            parameters.yrange = parameters.pixelsize*parameters.My/2;

            if strcmp(parameters.fitmodel,'aberrationsamp') 
                parameters.aberrations=aberrations;
                numparams_phaseaberr=length(aberrations)-3;
                numparams_ampaberr=length(aberrations);
                parameters.numparams=5+numparams_ampaberr+numparams_phaseaberr;
                parameters.numparams_phaseaberr=numparams_phaseaberr;
                parameters.numparams_ampaberr=numparams_ampaberr;
            elseif strcmp(parameters.fitmodel,'aberrations') 
                aberrations=aberrations(4:end,1:3);
                parameters.aberrations=aberrations;
                parameters.numparams=5+length(aberrations);
            end
            self.parameters=parameters;
        end
        function PSF=PSF(self)
            
            [~,~,wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
          get_pupil_matrix(self.parameters);
        
            parameters.zemit = 0;
            [~,~,FieldMatrix,FieldMatrixDerivatives] = ...
              get_field_matrix_derivatives(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,self.parameters);
             if true
               
                parameters.zemit = 2*d;
                [~,~,FieldMatrix1,FieldMatrixDerivatives1] = ...
                  get_field_matrix_derivatives(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters);
                FieldMatrixadd=FieldMatrix;
                FieldMatrixDerivativesadd=FieldMatrixDerivatives;
                for i=1:2
                    for j=1:3
                        for jz=1:n_z_slices
                            FieldMatrixadd{i,j,jz}=FieldMatrix{i,j,jz}+FieldMatrix1{i,j,jz};
                        end
                    end
                end
                for i=1:2
                    for j=1:3
                        for jz=1:n_z_slices
                            for k=1:20
                                FieldMatrixDerivativesadd{i,j,jz,k}=FieldMatrixDerivatives{i,j,jz,k}+FieldMatrixDerivatives1{i,j,jz,k};
                            end
                        end
                    end
                end
                self.parameters.zemit = 0;
            end
            [PSF,PSFderivatives] = get_psfs_derivatives(FieldMatrixadd,FieldMatrixDerivativesadd,self.parameters);
            
        end
    end

end


    


