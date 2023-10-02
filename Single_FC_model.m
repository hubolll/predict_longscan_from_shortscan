clc
clear

%% load information and data
load('D:\HCP_project\SA.txt');
path = 'D:\HCP_project\HCP_PTN1200\3T_HCP1200_MSMAll_d15_ts2';
list = dir([path,filesep,'*.txt']);
for i = 1:length(list)
    name = [path,filesep,list(i).name];
    TS(:,:,i) = load(name);
end

%% randomization
rd = load('D:\HCP_project\newtest\randomization.txt');
SA_test = SA(rd,2:3);
SA_train = SA(:,2:3);
SA_train(rd,:) = [];

%% calculation of FC
TP = []%time points
longscan = TS;
shortscan = longscan(1:TP,:,:);


for j = 1:1003
        FC_longscan = corrcoef(longscan(:,:,j));
        FC_longscan2(:,:,j) =FC_longscan;
        FC_shortscan = corrcoef(shortscan(:,:,j));
        FC_shortscan2(:,:,j) = FC_shortscan;
end

test_longscan = FC_longscan2(:,:,rd);
train_longscan = FC_longscan2;
train_longscan(:,:,rd) = [];

test_shortscan = FC_shortscan2(:,:,rd);
train_shortscan = FC_shortscan2;
train_shortscan(:,:,rd) = [];

%%model training
for i = 1:15
for j = (i+1):15
        
    train_longscan2 = train_longscan(i,j,:);
    train_longscan3 = reshape(train_longscan2,603,1);
    train_shortscan2 = train_shortscan(i,j,:);
    train_shortscan3 = reshape(train_shortscan2,603,1);
    a = regress(train_longscan3,[train_shortscan3,SA_train,ones(603,1)]);
    test_shortscan2 = test_shortscan(i,j,:);
    test_shortscan3 = reshape(test_shortscan2,400,1);
    test_longscan2 = test_longscan(i,j,:);
    test_longscan3 = reshape(test_longscan2,400,1);
    prediction = test_shortscan3*a(1)+SA_test(:,1)*a(2)+SA_test(:,2)*a(3)+a(4);
          r1= -abs(test_longscan3-prediction); %% the negative sign was used to make sure shorter absolute distance will derive a higher effect size
          r2= -abs(test_longscan3-test_shortscan3);
 es(i,j) = computeCohen_d(r1, r2, 'paired');
 
 predic2(i,j,:) = prediction;
end
end

%%calculation of mean effect size of single FC
 es2 = [es;zeros(1,15)];   
 ess = triu(es2);
 esss = ess(:);
 esss(esss==0) = [];
 essss = mean(esss);
 
%%calculation of mean effect size of summed FC
for f= 1:400
   sub = predic2(:,:,f);
   sum_pos(f) = sum(sum(sub.*(sub>0)));
   sum_neg(f) = sum(sum(sub.*(sub<0)));
   
   test_short2 = test_shortscan(:,:,f);
   test_short = triu(test_short2,1);
   test_short_pos(f) = sum(sum(test_short.*(test_short>0)));
   test_short_neg(f) = sum(sum(test_short.*(test_short<0)));
   
   test_long_real2 = test_longscan(:,:,f);
   test_long_real = triu(test_long_real2,1);
   test_long_pos(f) = sum(sum(test_long_real.*(test_long_real>0)));
   test_long_neg(f) = sum(sum(test_long_real.*(test_long_real<0)));
   

end
       pos1= -abs(test_long_pos-sum_pos);
       pos2= -abs(test_long_pos-test_short_pos);
       es_pos= computeCohen_d(pos1,pos2, 'paired');
          
       neg1= -abs(test_long_neg-sum_neg);
       neg2= -abs(test_long_neg-test_short_neg);
       es_neg= computeCohen_d(neg1,neg2, 'paired');
