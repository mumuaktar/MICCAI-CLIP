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
 
% [GLRLM] = getGLRLM(aa1(aa1~=0),levels);
% [textures] = getGLRLMtextures(GLRLM);
% [GLRLM] = getGLRLM(bb1(bb1~=0),levels);
% [textures1] = getGLRLMtextures(GLRLM);
% 
% names = fieldnames(textures); % extract names of features
%     for i = 1:numel(fieldnames(textures))
%     class(i) = textures.(names{i});
%     end
%     names = fieldnames(textures1); % extract names of features
%     for i = 1:numel(fieldnames(textures1))
%     class1(i) = textures1.(names{i});
%     end
%     m1=min(class,class1);
%     m2=max(class,class1);
%     f_glrlm=[(m1 ./ m2) abs(m1-m2)];
    
     
    %  f=[];
%      volume=niftiread(strcat(pre,num2str(num),fileformat));
%      left=volume .* m_l;
%     right=volume .* m_r;
%   mask=niftiread('ctamask.nii');
% %  [ROIonly,levels] = prepareVolume(volume,mask,'Other',4,3.27,1,5,'Matrix','Lloyd',32);
%    [ROIonly,levels] = prepareVolume(left,mask,'Other',4,3.27,1,5,'Matrix','Lloyd',32);
%     [ROIonly1,levels] = prepareVolume(right,mask,'Other',4,3.27,1,5,'Matrix','Lloyd',32); 
%      [r,c,d]=size(volume);
%      if d~=51
%          vv=volume(:,:,65:115);
%          v1=vv;
%      else
%          v1=volume;
%      end
%  
%     p='\\perf-loy-nas.concordia.ca\home\home\m_ktar\brain_registration_ncct\brain_extraction_fsl';
%    outputFileName = fullfile('l_quantized', ['sub_l' num2str(num) '.nii']);
%       outputFileName1 = fullfile('r_quantized', ['sub_r' num2str(num) '.nii']);
% %    outputFileName1 = fullfile('quantized', ['sub' num2str(num) '.nii']);
%      niftiwrite(ROIonly, outputFileName); 
%      niftiwrite(ROIonly1, outputFileName1);
% % save(outputFileName,'levels');
%  end


 
%      
%      volume=imresize3(volume,[221 221 160]);
%      volume(volume>=100)=0;
%        volume=impyramid(volume, 'reduce');
%       max(volume(:))
%       min(volume(:))
%       [row,col,dim]=size(volume);
%        volume=imresize3(volume,[row col 80]);
% %      volume=medfilt3(volume);
% %volume=volume(:,:,75:110);
%     volume(volume<0)=0;
% %      volume=im2double(volume);
% [row,col,dim]=size(volume);
% % volume=imgaussfilt3(volume);
% % v_ax=max(volume,[],3);
% % v_cor=max(volume,[],1);
% % v_cor=squeeze(v_cor);
% % v_sag=max(volume,[],2);
% % v_sag=squeeze(v_sag);
% % 
% % for i=1:3
% % LBP(:,i)=extractLBPFeatures(v_ax,'Upright',false);
% % LBP(:,i)=extractLBPFeatures(v_cor,'Upright',false);
% % LBP(:,i)=extractLBPFeatures(v_sag,'Upright',false);
% % end
% % LBP1=mean(LBP,2);
% 
% 
% 
% %       m_l=imresize3(m_l,[row col dim]);
% %      m_r=imresize3(m_r,[row col dim]);
% % %      
% % % %      
% %     left=volume .* m_l;
% %     right=volume .* m_r;
% %     
% %     volume=left;
%  [rows,cols,dim]=size(volume);
%   volume=im2double(volume); 
% mask=niftiread('ctamask.nii');
%  [ROIonly,levels] = prepareVolume(volume,mask,'Other',4,3.27,1,5,'Matrix','Uniform',32);
% 
aa1=ROIonly;
bb1=ROIonly1;
[GLCM] = getGLCM(aa1(aa1~=0),levels);
    [textures] = getGLCMtextures(GLCM);
    [GLCM] = getGLCM(bb1(bb1~=0),levels);
[textures1] = getGLCMtextures(GLCM);
    
% %      
% % 
% % %  local1=[];
% % % for i = 0: rows/8 - 2
% % %      for j= 0: cols/8 -2
% % %        for k=0: dim/8 -2 
% % %          block = volume(8*i+1 : 8*i+16 , 8*j+1 : 8*j+16, 8*k+1 :8*k+16);
% % %           for x= 0:1
% % %             for y= 0:1
% % %                 for z=0:1
% % %                 cell   =block(8*x+1:8*x+8, 8*y+1:8*y+8, 8*z+1:8*z+8); 
% % %                 local  =zeros(1,5);
% % %                 m=mean(cell(:));
% % %                 v=var(cell(:));
% % %                 s=skewness(cell(:));
% % %                 k=kurtosis(cell(:));
% % %                 st=std(cell(:));
% % %                 md=median(cell(:));
% % %                 mdd=mode(cell(:));
% % %                 
% % %                 e=entropy(cell(:));
% % %                 mn=min(cell(:));
% % %                 mx=max(cell(:));
% % %                 local=[m v s k st md mdd e mn mx];
% % %                 
% % %                 end
% % %             end
% % %           end
% % %        local1=[local1 local'];            
% % % end
% % % end
% % % end
% % %     lc=mean(local1,2);
% % % %     
% % %     volume=right;
% % %     [rows,cols,dim]=size(volume);
% % %  local1=[];
% % % for i = 0: rows/8 - 2
% % %      for j= 0: cols/8 -2
% % %        for k=0: dim/8 -2 
% % %          block = volume(8*i+1 : 8*i+16 , 8*j+1 : 8*j+16, 8*k+1 :8*k+16);
% % %           for x= 0:1
% % %             for y= 0:1
% % %                 for z=0:1
% % %                 cell   =block(8*x+1:8*x+8, 8*y+1:8*y+8, 8*z+1:8*z+8); 
% % %                 local  =zeros(1,5);
% % %                 m=mean(cell(:));
% % %                 v=var(cell(:));
% % %                 s=skewness(cell(:));
% % %                 k=kurtosis(cell(:));
% % %                 s=std(cell(:));
% % % %                 e=entropy(cell(:));
% % % %                 mn=min(cell(:));
% % % %                 mx=max(cell(:));
% % %                 local=[m v s k s];
% % %                 
% % %                 end
% % %             end
% % %           end
% % %        local1=[local1 local'];            
% % % end
% % % end
% % % end
% % %     l2=mean(local1,2);
% %     
% %     
% %      %vv=max(volume,[],3);
% %      %vv=imrotate(vv,90);
% %      %figure,imagesc(vv),colormap(gray)
% % theta_histogram_bins = 9;
% % phi_histogram_bins = 18;
% % cell_size = 16;
% % block_size = 2;
% % step_size = 2;
% %  volume1=rescale(volume);
% % features = hog3d(volume1, cell_size, block_size, theta_histogram_bins, phi_histogram_bins, step_size);
% % % %features = hog3d(Bottle_volume, cell_size, block_size, theta_histogram_bins, phi_histogram_bins, step_size);
% % % %VISUALISE THE FEATURES
% % % plot_hog3d(features, cell_size, theta_histogram_bins, phi_histogram_bins);
% % % view([45,45]);
% % % a1=[];
% % for i=1:phi_histogram_bins
% % a=features(i).Features;
% % aa=mean(a);
% % a1=[a1; aa];
% %  end
% % ff=mean(a1);
% [row,col,dim]=size(aa1);
%   for i=1:dim
%    LBP(i,:)=extractLBPFeatures(aa1(:,:,i),'Upright',false);
% end
% LBP1=mean(LBP); 
% 
% for i=1:dim
%    LBP_r(i,:)=extractLBPFeatures(bb1(:,:,i),'Upright',false);
% end
% LBP1_r=mean(LBP_r);
% 
% m1=min(LBP1,LBP1_r);
%      m2=max(LBP1,LBP1_r);    
%      ratio_lbp=(m1 ./ m2);
%      dif_lbp=abs(m1-m2);


% % %      [row,col,dim]=size(volume);
% % % for i=1:80
% % % lbpFeatures = extractLBPFeatures(volume(:,:,i),'NumNeighbors',26);
% % % lbp(i,:)=lbpFeatures;
% % % numNeighbors = 16;
% % % numBins = numNeighbors*(numNeighbors-1)+3;
% % % lbpCellHists = reshape(lbpFeatures,numBins,[]);
% % % lbpCellHists = bsxfun(@rdivide,lbpCellHists,sum(lbpCellHists));
% % % lbpCellHists=mean(lbpCellHists);
% % % lbpFeatures = reshape(lbpCellHists,1,[]);
% % %  lbp(i,:)=lbpFeatures;
% % % end
% % % LBP1=mean(lbp);
% % % % 
% %         %mask=niftiread('ctamask.nii');
% %  %[ROIonly,levels] = prepareVolume(volume,mask,'Other',4,3.27,1,5,'Matrix','Uniform',32);
% %  [GLRLM] = getGLRLM(ROIonly,levels);
% % [textures] = getGLRLMtextures(GLRLM);
% % % %  
names = fieldnames(textures); % extract names of features
    for i = 1:numel(fieldnames(textures))
    class2(i) = textures.(names{i});
    end
% % % % %     
% % % %    glrlm=class;
% % % % %   
% % % %      [GLCM] = getGLCM(ROIonly,levels);
% % % %      [textures] = getGLCMtextures(GLCM);
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
% % %     
% %    glrlm=class;
% % %   
% %      [GLCM] = getGLCM(ROIonly,levels);
% %      [textures] = getGLCMtextures(GLCM);
% %      
     names = fieldnames(textures1); % extract names of features
    for i = 1:numel(fieldnames(textures1))
    class5(i) = textures1.(names{i});
    end
    
   m1=min(class4,class5);
    m2=max(class4,class5);
    f_glszm=[ abs(m1-m2)];
% %      
% % %      names = fieldnames(textures); % extract names of features
% % %     for i = 1:numel(fieldnames(textures))
% % %     class2(i) = textures.(names{i});
% % %     end
% % %     
% % %    glszm=class2;
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
    
%     for i = 1:numel(fieldnames(textures))
%     class4(i) = textures.(names{i});
%     end
% %     
% %    glo=class4;
% % % % %      
% % % % % %        
% % % % [ROIonly,levels] = prepareVolume(right,mask,'Other',4,3.27,1,5,'Matrix','Uniform',32);
% % % %  [GLRLM] = getGLRLM(ROIonly,levels);
% % % % [textures] = getGLRLMtextures(GLRLM);
% % % %  
% % % % names = fieldnames(textures); % extract names of features
% % % %     for i = 1:numel(fieldnames(textures))
% % % %     class_r(i) = textures.(names{i});
% % % %     end
% % % % %     
% % % %    glrlm_r=class_r;
% % % %   
% % % %      [GLCM] = getGLCM(ROIonly,levels);
% % % %      [textures] = getGLCMtextures(GLCM);
% % % %      
% % % %      names = fieldnames(textures); % extract names of features
% % % %     for i = 1:numel(fieldnames(textures))
% % % %     class1_r(i) = textures.(names{i});
% % % %     end
% % % % %     
% % %    glcm_r=class1_r;
% % %      
% % % 
% % % 
% % % 
% % %    [GLSZM] = getGLSZM(ROIonly,levels);
% % %      [textures] = getGLSZMtextures(GLSZM); 
% % %      names = fieldnames(textures); % extract names of features
% % %     for i = 1:numel(fieldnames(textures))
% % %     class2_r(i) = textures.(names{i});
% % %     end
% % %     
% % %    glszm_r=class2_r;
% % % %      
% % %       [NGTDM,countValid] = getNGTDM(ROIonly,levels);
% % %      [textures] = getNGTDMtextures(NGTDM,countValid); 
% % %      names = fieldnames(textures); % extract names of features
% % %     for i = 1:numel(fieldnames(textures))
% % %     class3_r(i) = textures.(names{i});
% % %     end
% % %     
% % %    ngtdm_r=class3_r;
% % % %      
% % %       ROIonly=volume;
% % %  [textures] = getGlobalTextures(ROIonly,9);
% % %  
% % %  names = fieldnames(textures); % extract names of features
% % %     for i = 1:numel(fieldnames(textures))
% % %     class4_r(i) = textures.(names{i});
% % %     end
% % %     
% % %    glo_r=class4_r;
% % % % var_r=textures.Variance;
% % % % skew_r=textures.Skewness;
% % % % kurt_r=textures.Kurtosis;
% % % % % 
% % % %  for i = 1:9
% % % %     dif(i) = abs(class1(i)-class(i));
% % % %     end
% % %   mean_l=mean(left(:));
% % %   mean_r=mean(right(:));
% % %   dif_intensity=abs(mean_l-mean_r);
% % % % std_l=std(left(:));
% % % % std_r=std(right(:));
% % % % d_std=abs(std_l-std_r);
% % % % 
% % % % m1=max(volume(:));
% % % % m2=min(volume(:));
% % % % aa=max(volume,[],3);
% % % % lbp=extractLBPFeatures(aa,'Upright',false);
% % % % 
% % 
% % % % % var_d=abs(var_l-var_r);
% % % % % skew_d=abs(skew_l-skew_r);
% % % % % kurt_d=abs(kurt_l-kurt_r);
% % % % % 
% %  %  f=[abs(glrlm-glrlm_r) abs(glcm-glcm_r) abs(glszm-glszm_r) abs(ngtdm-ngtdm_r) abs(glo-glo_r)];
% %      f=[glrlm glcm glszm ngtdm glo ff ];
% % % f=[LBP1];
% % %f=[ff];
% %    feature=[feature f'];
% %    toc
% % % % 
% % % %     [NGTDM,countValid] = getNGTDM(ROIonly,levels);
% % % %    [textures] = getNGTDMtextures(NGTDM,countValid);
% % % %     names = fieldnames(textures); % extract names of features
% % % %     for i = 1:5
% % % %     ngt(i) = textures.(names{i});
% % % %     end
% % % %    
% % %  f=[class class1 class2 class3 class4];
% % % feature=[feature f'];
% % % % %      Con_l=textures.Contrast;
% % % % %      Cor_l=textures.Correlation;
% % % % %      SumAverage_l=textures.SumAverage;
% % % % %      auto_l=textures.AutoCorrelation;
% % % % %      Dis_l=textures.Dissimilarity;
% % % % %      En_l=textures.Energy;
% % % % %      var_l=textures.Variance;
% % % %      Hom_l=textures.Homogeneity;
% % % %      Ent_l=textures.Entropy;
% % %       [ROIonly,levels] = prepareVolume(right,mask,'Other',4,3.27,1,5,'Matrix','Uniform',32);
% % % %      [GLCM] = getGLCM(ROIonly,levels);
% % % %      [textures] = getGLCMtextures(GLCM);
% % % %     [GLRLM] = getGLRLM(ROIonly,levels);
% % % %     [textures] = getGLRLMtextures(GLRLM);
% % 
% % 
% % %      Con_r=textures.Contrast;
% % %      Cor_r=textures.Correlation;
% % %      SumAverage_r=textures.SumAverage;
% % %      auto_r=textures.AutoCorrelation;
% % %      Dis_r=textures.Dissimilarity;
% % %      En_r=textures.Energy;
% % %      var_r=textures.Variance;
% % %      Hom_r=textures.Homogeneity;
% % %      Ent_r=textures.Entropy;
% % %      f=[abs(Con_l-Con_r) abs(Cor_l-Cor_r) abs(SumAverage_l-SumAverage_r) abs(Dis_l-Dis_r) abs(En_l-En_r) abs(Ent_l-Ent_r) abs(var_l-var_r) abs(Hom_l-Hom_r) abs(auto_l-auto_r)];
% % 
% % % ROIonly=left;
% % % [textures] = getGlobalTextures(ROIonly,32);
% % % var_l=textures.Variance;
% % % skew_l=textures.Sksewness;
% % % kurt_l=textures.Kurtosis;
% % % 
% % 
% % 
% m11=mean(aa1(:));
% m22=mean(bb1(:));
% m=abs(m11-m22);
% s11=std(aa1(:));
% s22=std(bb1(:));
% s=abs(s11-s22);
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
