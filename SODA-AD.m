clc;clear;

n=5;
g=4;
TAR=0.2;
RDD=0.6;
load (strcat('problems\prob_.mat'));


int_dim=n*g+n*g;
con_dim=n*g+n*(g-1)+n*(g-1)+n*(g-1);
Dim=int_dim+con_dim;

Data.number_startpoints=2*(Dim+1);
Data.dim=Dim;

Data.xlow=[ones(1,int_dim),zeros(1,con_dim)];
Data.xup=[3, 5, 5, 4];

Data.xup=[repmat(Data.xup,1,n),l*ones(1,n*g),ones(1,con_dim) ];
Data.integer=[1:int_dim]; %indices of integer variables
Data.continuous=[int_dim+1:Dim]; %indices of continuous variables
InitialPoints = slhd(Data);
xlatin=repmat(Data.xlow, Data.number_startpoints,1) + repmat(Data.xup-Data.xlow,Data.number_startpoints,1).*InitialPoints; 
Data.S=xlatin;
fixrank = false;  
if rank([Data.S,ones(size(Data.S,1),1)]) < Data.dim + 1
    fixrank = true;
end 
while fixrank %rank is too small to fit RBF, add random points to initial design to satisfy rank condition
    n_new = Data.dim+1-rank([Data.S,ones(size(Data.S,1),1)]); %minimum number of new points needed
    randpoint = repmat(Data.xlow,n_new,1) + repmat((Data.xup-Data.xlow), n_new,1).*rand(n_new,Data.dim);
    temp=rank([[Data.S;randpoint], ones(size(Data.S,1)+n_new,1)]);
    if temp == Data.dim + 1
        Data.S = [Data.S; randpoint];
        fixrank = false;
    end 
end
Data.S(:,Data.integer)=round(Data.S(:,Data.integer));


save (strcat('problems\Data1.mat'));


[xbest, fbest] = miso('datainput_fun',1000, 'rbf_c', [], 'slhd', 'cp1',[],Data); %SODA-AD

Data.save_name=strcat('sodaad','_ans\prob_.mat');
