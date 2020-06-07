clear all;
clc;
addpath(genpath('./data/'))
addpath(genpath('./tools/'))

L = 0.5;

lambda_row = 0.01;
lambda_col = 0.01;
maxiter=2000;
maxProxIter=2000;
xi = "zero";
%dataset = "secom";
%dataset = "wdbc";
dataset = "musk";
%dataset = "lena256";
%dataset = "madelon";
relatives = [[0.9;0] [0;0.9],[0.5;0.5]];
names = ["distance sigma = "+relatives(1,1) "gradient zeta = "+relatives(2,2) "distance sigma = "+relatives(1,3) + "gradient zeta = "+relatives(2,3)  ];
[~,n] = size(relatives);
infos = cell(1,n);


minf = Inf;
for i = 1:n
   r = relatives(:,i);
   sigma = r(1);
   zeta = r(2);
   infos{i} = CURexp(dataset,maxiter,lambda_row,lambda_col,maxProxIter,sigma,zeta,"zero",L);
   minf = min(minf,min(infos{i}(:,1)));
   
end    

infosa = cell(1,3);
for i = 1:3
xi = "poly"+(i+1);    
infosa{i} = CURexp(dataset,maxiter,lambda_row,lambda_col,maxProxIter,0,0,xi,L);   
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

resultDir = sprintf('./CUR_results/');
if ~exist(resultDir,'dir')
	mkdir(resultDir);
end
path = "CUR_results/CUR_"+dataset+"_"+L+"_"+lambda_row+"_"+maxiter+"_";

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

% 
% 
% 
% sigma = 0.01;
% zeta = 0.01;
% 
% info2 = CURexp("wdbc",maxiter,lambda_row,lambda_col,maxProxIter,sigma,zeta,xi,L,mu);
% 
% hold on;
% 
% semilogy(cumsum(info2(:,2)),info2(:,1)-min(info2(:,1)),'Linewidth',2);
% 
% 
% sigma = 0.01;
% zeta = 0.00;
% 
% info3 = CURexp("wdbc",maxiter,lambda_row,lambda_col,maxProxIter,sigma,zeta,xi,L,mu);
% 
% hold on;
% 
% semilogy(cumsum(info3(:,2)),info3(:,1)-min(info3(:,1)),'Linewidth',2);