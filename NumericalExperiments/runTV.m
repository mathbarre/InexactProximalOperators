clear all;
clc;
addpath(genpath('./data/'))
addpath(genpath('./tools/'))

L = 0.5;

lambda_reg = 0.1;

maxiter=3000;
maxProxIter=3000;
xi = "zero";
noise = 0.1*randn(256,256);
blur = ones(5)/25;
dataset = "boat";

relatives = [[0.9;0],[0.0;0.9],[0.5;0.5]];
names = ["distance sigma = "+relatives(1,1) "gradient zeta = "+relatives(2,2) "distance sigma = "+relatives(1,3) + "gradient zeta = "+relatives(2,3)  ];
[~,n] = size(relatives);
infos = cell(1,n);
X = cell(1,n);
minf = Inf;
for i = 1:n
   r = relatives(:,i);
   sigma = r(1);
   zeta = r(2);
   
   [infos{i},X{i}] = TVexp(dataset,blur,noise,maxiter,lambda_reg,maxProxIter,sigma,zeta,"zero",L);
   subplot(1,2,2)
   imshow(X{i})
   title(names(i));
   minf = min(minf,min(infos{i}(:,1)));
   
end    




infosa = cell(1,3);
Xa = cell(1,3);
for i = 1:3
xi = "poly"+(i+1);    
[infosa{i},Xa{i}] = TVexp(dataset,blur,noise,maxiter,lambda_reg,maxProxIter,0,0,xi,L);
subplot(1,2,2)
imshow(Xa{i});
title(xi);
minf = min(minf,min(infosa{i}(:,1)));    
end    




%minf = min(minf,min(infoBello(:,1)));



figure()
for i = 1:n
   semilogy(cumsum(infos{i}(:,2)),infos{i}(:,1)-minf,'Linewidth',2);
   hold on;
end    


for i = 1:3
    semilogy(cumsum(infosa{i}(:,2)),infosa{i}(:,1)-minf,'Linewidth',2);
end    



legend([names "poly2" "poly3" "poly4"]);

resultDir = sprintf('./TV_results/');
if ~exist(resultDir,'dir')
	mkdir(resultDir);
end
path = "TV_results/TV_"+dataset+"_"+L+"_"+lambda_reg+"_"+maxiter+"_";

for i = 1:n
   tot = [cumsum(infos{i}(:,2)) infos{i}(:,1)-minf];  
   r = relatives(:,i);
   sigma = r(1);
   zeta = r(2);
   fid = fopen(path+sigma+"_"+zeta+"_zero.txt",'wt');
   fprintf(fid,'%s\t',["k" "f"]);
   fprintf(fid,'\n');
    for ii = 1:size(tot,1)
        fprintf(fid,'%.16g\t',tot(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end    

for i = 1:n
   tot = [cumsum(infosa{i}(:,2)) infosa{i}(:,1)-minf];  
   fid = fopen(path+"0_0_poly"+(i+1)+".txt",'wt');
   fprintf(fid,'%s\t',["k" "f"]);
   fprintf(fid,'\n');
    for ii = 1:size(tot,1)
        fprintf(fid,'%.16g\t',tot(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end    


