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

TP = []%time point
longscan = TS;
shortscan = longscan(1:TP,:,:);

longscan_test = TS(:,:,rd);
longscan_train = TS;
longscan_train(:,:,rd) = [];

shortscan_test = shortscan(:,:,rd);
shortscan_train = shortscan;
shortscan_train(:,:,rd) = [];
%%calculation of dynamic FC
width = 40;
leth = 3;
k = 0;
for j = 1:603
    for i = width:leth:600
        DFC_short_train(j,:,:,(k+1)) = corr(shortscan_train(((k*leth)+1):i,:,j));
        k = k+1;
    end
    k=0;
end

k = 0;
for j = 1:400
    for i = width:leth:600
        DFC_short_test(j,:,:,(k+1)) = corr(shortscan_test(((k*leth)+1):i,:,j));
        k = k+1;
    end
    k=0;
end

%%calculation of static FC
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

size1 = size(DFC_short_train);
size2 = size(DFC_short_test);

%%model training
for i = 1:15
for j = (i+1):15
        
    train_longscan2 = train_longscan(i,j,:);
    train_longscan3 = reshape(train_longscan2,603,1);
    train_dynamic2 = DFC_short_train(:,i,j,:);
    train_dynamic = reshape(train_dynamic2,603,size1(4));
    a = regress(train_longscan3,[train_dynamic,SA_train,ones(603,1)]);
    test_shortscan2 = test_shortscan(i,j,:);
    test_shortscan3 = reshape(test_shortscan2,400,1);
    test_longscan2 = FC_long_test(i,j,:);
    test_longscan3 = reshape(test_longscan2,400,1);
    test_dynamic2 = DFC_short_test(:,i,j,:);
    test_dynamic = reshape(test_dynamic2,400,size2(4));
    prediction = sum([test_dynamic,SA_test,ones(400,1)].*a',2);
    r= -abs(test_longscan3-prediction);
    r2= -abs(test_longscan3-test_shortscan3);
    es(i,j) = computeCohen_d(r, r2, 'paired');
 
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

