function T=asl_analysis(src,subj)


[cvr_default, dims, scales, bpp, endian]= read_avw([src subj '/func/' subj '_cvr_default']);
[cvr_paramA, dims, scales, bpp, endian] = read_avw([src subj '/func/' subj '_cvr_paramA']);
[cvr_paramB, dims, scales, bpp, endian] = read_avw([src subj '/func/' subj '_cvr_paramB']);
gm_cortex=read_avw([src subj '/func/' subj '_GM_MZeroScan']);

figure;
imagesc(cvr_default(:,:,6),[0 2]);
colorbar;

cvr_defaultr=reshape(cvr_default,64*64*12,1);
cvr_paramAr=reshape(cvr_paramA,64*64*12,1);
cvr_paramBr=reshape(cvr_paramB,64*64*12,1);
gm_cortexr=reshape(gm_cortex,64*64*12,1);;

cvr_default_vals=cvr_defaultr(gm_cortexr>0,:);
cvr_paramA_vals=cvr_paramAr(gm_cortexr>0,:);
cvr_paramB_vals=cvr_paramBr(gm_cortexr>0,:);

rows={'default'; 'paramA'; 'paramB'};

cvrmeans(1,:)=mean(cvr_default_vals);
cvrmeansnoZ(1,:)=mean(cvr_default_vals(cvr_default_vals~=0));
cvrmedians(1,:)=median(cvr_default_vals);
cvrmediansnoZ(1,:)=median(cvr_default_vals(cvr_default_vals~=0));
cvrmeans(2,:)=mean(cvr_paramA_vals);
cvrmeansnoZ (2,:)=mean(cvr_paramA_vals(cvr_paramA_vals~=0));
cvrmedians(2,:)=median(cvr_paramA_vals);
cvrmediansnoZ(2,:)=median(cvr_paramA_vals(cvr_paramA_vals~=0));
cvrmeans(3,:)=mean(cvr_paramB_vals);
cvrmeansnoZ(3,:)=mean(cvr_paramB_vals(cvr_paramB_vals~=0));
cvrmedians(3,:)=median(cvr_paramB_vals);
cvrmediansnoZ(3,:)=median(cvr_paramB_vals(cvr_paramB_vals~=0));

T=table(rows,cvrmeans,cvrmeansnoZ,cvrmedians,cvrmediansnoZ);

figure;
hist(cvr_default_vals,1000)
title('CVR default');

figure;
hist(cvr_paramA_vals,1000)
figure;
hist(cvr_paramB_vals,1000)