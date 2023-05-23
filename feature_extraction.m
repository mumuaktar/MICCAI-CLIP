%%%%%%%%%%%prepare volume function and radiomic feature extraction codes
%%%%%%%%%%%are from this link: https://github.com/mvallieres/radiomics
%%%%A package providing MATLAB programming tools for radiomics analysis.



clc;clear;
NUM=64;
fileformat='.nii.gz';
pre='sub_l_d';
pre1='sub_r_d';
prefix='f_sub_brain';
prefix1='right_hem';
feature=[];

m_l=niftiread('seg_left.nii');
m_r=niftiread('seg_right.nii');

c=niftiread('warped_final_m1.nii.gz');
for num=35:35
    tic
     volume=niftiread(strcat(prefix,num2str(num),fileformat));
     left=volume .* m_l;
     right=volume .* m_r;
     mask=niftiread('ctamask.nii');
     aa1=left;
     bb1=right;
     d1=abs((aa1-mean(aa1(:)))/std(aa1(:)));
     d2=abs((bb1-mean(bb1(:)))/std(bb1(:)));


     [ROIonly,levels] = prepareVolume(d1,mask,'Other',4,3.27,1,5,'Matrix','Uniform',32);
 
    [ROIonly1,levels] = prepareVolume(d2,mask,'Other',4,3.27,1,5,'Matrix','Uniform',32); 
 
     BW1 = edge3(aa1,'Sobel',.01);
     BW2 = edge3(bb1,'Sobel',.01);
     aa1(~BW1)=0;
     bb1(~BW2)=0;
 

    aa1=ROIonly;
    bb1=ROIonly1;
    [GLCM] = getGLCM(aa1(aa1~=0),levels);
    [textures] = getGLCMtextures(GLCM);
    [GLCM] = getGLCM(bb1(bb1~=0),levels);
    [textures1] = getGLCMtextures(GLCM);



    names = fieldnames(textures); % extract names of features
    for i = 1:numel(fieldnames(textures))
    class2(i) = textures.(names{i});
    end

% % % %      
     names = fieldnames(textures1); % extract names of features
    for i = 1:numel(fieldnames(textures1))
    class3(i) = textures1.(names{i});
    end
% % %     
   m1=min(class2,class3);
    m2=max(class2,class3);
    f_glcm=[ abs(m1-m2)];
% % % % % %     
% % % %    glcm=class1;
% % %        
     [GLSZM] = getGLSZM(aa1(aa1~=0),levels);
     [textures] = getGLSZMtextures(GLSZM); 
% % 

   [GLSZM] = getGLSZM(bb1(bb1~=0),levels);
     [textures1] = getGLSZMtextures(GLSZM); 
     
     
     names = fieldnames(textures); % extract names of features
    for i = 1:numel(fieldnames(textures))
    class4(i) = textures.(names{i});
    end

     names = fieldnames(textures1); % extract names of features
    for i = 1:numel(fieldnames(textures1))
    class5(i) = textures1.(names{i});
    end
    
   m1=min(class4,class5);
    m2=max(class4,class5);
    f_glszm=[ abs(m1-m2)];

        [NGTDM,countValid] = getNGTDM(aa1(aa1~=0),levels);
     [textures] = getNGTDMtextures(NGTDM,countValid); 
% % %      

      [NGTDM,countValid] = getNGTDM(bb1(bb1~=0),levels);
     [textures1] = getNGTDMtextures(NGTDM,countValid); 
     names = fieldnames(textures); % extract names of features
    for i = 1:numel(fieldnames(textures))
    class6(i) = textures.(names{i});
    end
% %     
% %    glrlm=class;
% % %   
% %      [GLCM] = getGLCM(ROIonly,levels);
% %      [textures] = getGLCMtextures(GLCM);
% %      
     names = fieldnames(textures1); % extract names of features
    for i = 1:numel(fieldnames(textures1))
    class7(i) = textures1.(names{i});
     end
%     
   m1=min(class6,class7);
    m2=max(class6,class7);
    f_ngtdm=[abs(m1-m2)];
     
% %      names = fieldnames(textures); % extract names of features
% %     for i = 1:numel(fieldnames(textures))
% %     class3(i) = textures.(names{i});
% %     end
% %     
% %    ngtdm=class3;
% % % %      
% %       ROIonly=volume;
 [textures] = getGlobalTextures(aa1(aa1~=0),9);
 
 [textures1] = getGlobalTextures(bb1(bb1~=0),9);
% [textures] = getGlobalTextures(d1,9);
%  
%  [textures1] = getGlobalTextures(d2,9);
%  
  names = fieldnames(textures); % extract namnames = fieldnames(textures); % extract names of features
    for i = 1:numel(fieldnames(textures))
    class8(i) = textures.(names{i});
    end
% % %     
%    glrlm=class;
% %   
%      [GLCM] = getGLCM(ROIonly,levels);
%      [textures] = getGLCMtextures(GLCM);
% %      
     names = fieldnames(textures1); % extract names of features
    for i = 1:numel(fieldnames(textures1))
    class9(i) = textures1.(names{i});
    end
    
   m1=min(class8,class9);
    m2=max(class8,class9);
    f_glo=[ abs(m1-m2)];
    


 f=[f_glcm f_glszm f_glo f_ngtdm];

feature=[feature f'];
toc
 end
%  
  feature=feature';
%  
% %    feature=reshape(feature,64,10);
    feature(isnan(feature))=0;
  save('cotrast1.mat','feature');
  clc;clear;
  load('chk1.mat')
sample = ff;
%rowNames = {'a','b','c'};
colNames = {'d_Energy','d_Contrast','d_Entropy','d_Homogeneity','d_Correlation','d_Variance','d_SumAverage','d_Dissimilarity','d_Autocorrelation','d_SZE','d_LZE','d_GLN','d_ZSN','d_ZP','d_LGZE','d_HGZE','d_SZLGE','d_SZHGE','d_LZLGE','d_LZHGE','d_GLV','d_ZSV','d_Mean','d_G_Energy','d_G_Entropy','d_G_Variance','d_Skewness','d_Kurtosis', 'c_Energy','c_Contrast','c_Entropy','c_Homogeneity','c_Correlation','c_Variance','c_SumAverage','c_Dissimilarity','c_Autocorrelation','c_SZE','c_LZE','c_GLN','c_ZSN','c_ZP','c_LGZE','c_HGZE','c_SZLGE','c_SZHGE','c_LZLGE','c_LZHGE','c_GLV','c_ZSV','c_Mean','c_G_Energy','c_G_Entropy','c_G_Variance','c_Skewness','c_Kurtosis'};

sTable = array2table(sample,'VariableNames',colNames);
writetable(sTable,'feature_set.csv')  
