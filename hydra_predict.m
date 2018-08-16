%  HYDRA PREDICT
%  Version 1.0.0 --- August 2018
%  Section of Biomedical Image Analysis
%  Department of Radiology
%  University of Pennsylvania
%  Richard Building
%  3700 Hamilton Walk, 7th Floor
%  Philadelphia, PA 19104
%
%  Web:   https://www.med.upenn.edu/sbia/
%  Email: sbia-software at uphs.upenn.edu
%
%  Copyright (c) 2018 University of Pennsylvania. All rights reserved.
%  See https://www.med.upenn.edu/sbia/software-agreement.html or COPYING file.
%
%  Author:
%  Erdem Varol
%  software@cbica.upenn.edu


function CIDX = hydra_predict(varargin)

if nargin==0
    printhelp()
    return
end

if( strcmp(varargin{1},'--help') || isempty(varargin))
    printhelp()
    return;
end

if( strcmp(varargin{1},'-h') || isempty(varargin) )
    printhelp()
    return
end

if( strcmp(varargin{1},'--version') || isempty(varargin) )
    fprintf('Version 1.0.\n')
    return
end

if( strcmp(varargin{1},'-v') || isempty(varargin) )
    fprintf('Version 1.0.\n')
    return
end

if( strcmp(varargin{1},'-u') || isempty(varargin) )
    fprintf(' EXAMPLE USE (in matlab) \n');
    fprintf(' hydra_predict(''-i'',''test.csv'',''-o'',''.'',''-m'',''hydra_model.mat'') \n');
    fprintf(' EXAMPLE USE (in command line) \n');
    fprintf(' hydra_predict -i test.csv -o . -m hydra_model.mat \n');
    return
end

if( strcmp(varargin{1},'--usage') || isempty(varargin) )
    fprintf(' EXAMPLE USE (in matlab) \n');
    fprintf(' hydra_predict(''-i'',''test.csv'',''-o'',''.'',''-m'',''hydra_model.mat'') \n');
    fprintf(' EXAMPLE USE (in command line) \n');
    fprintf(' hydra_predict -i test.csv -o . -m hydra_model.mat \n');
    return
end

if( sum(or(strcmpi(varargin,'--input'),strcmpi(varargin,'-i')))==1)
    featureCSV=varargin{find(or(strcmpi(varargin,'--input'),strcmp(varargin,'-i')))+1};
else
    error('hydra:argChk','Please specify input csv file!');
end


if( sum(or(strcmpi(varargin,'--outputDir'),strcmpi(varargin,'-o')))==1)
    outputDir=varargin{find(or(strcmp(varargin,'--outputDir'),strcmp(varargin,'-o')))+1};
else
    error('hydra:argChk','Please specify output directory!');
end

if( sum(or(strcmpi(varargin,'--model'),strcmpi(varargin,'-m')))==1)
    modelFile=varargin{find(or(strcmp(varargin,'--model'),strcmp(varargin,'-m')))+1};
else
    error('hydra:argChk','Please specify model file!');
end

% create output directory
if (~exist(outputDir,'dir'))
    [status,~,~] = mkdir(outputDir);
    if (status == 0)
        error('hydra:argChk','Cannot create output directory!');
    end
end

% csv with features
fname=featureCSV;
if (~exist(fname,'file'))
    error('hydra:argChk','Input feature .csv file does not exist');
end

% input data
% assumption is that the first column contains IDs, and the last contains
% labels
disp('Loading features...');
input=readtable(fname);
ID=input{:,1};
XK=input{:,2:end};
HYDRA_model=load(modelFile);

CIDX = zeros(size(XK,1),length(HYDRA_model.fullmodel));
for k=1:length(HYDRA_model.fullmodel)
    clear S
    S = zeros(size(XK,1),k);
for j=1:HYDRA_model.fullmodel{k}.params.k
    S(:,j)=w_svmpredict(XK,HYDRA_model.fullmodel{k}.mdl{j},HYDRA_model.fullmodel{k}.params.kernel);
end
[~,CIDX(:,k)] = max(S,[],2);
end


disp('Saving results...')
save([outputDir '/HYDRA_predict_results.mat'],'CIDX','ID');
disp('Done')

end


function S=w_svmpredict(X,mdl,dual)
%Function makes svm prediction using model
if dual==0
    X=X;
elseif dual==1
    [U,S,~]=svd(X);
    X=U*sqrt(S);
    
end

S=X*mdl.w+mdl.b;

end


function printhelp()
fprintf(' function returns estimated subgroups of out of sample subjects \n')
fprintf(' through the clustering assignments obtained by HYDRA\n')
fprintf('\n')
fprintf(' INPUT\n')
fprintf('\n')
fprintf(' REQUIRED\n')
fprintf(' [--input, -i] : .csv file containing the input features. (REQUIRED)\n')
fprintf('              every column of the file contains values for a feature, with\n')
fprintf('              the exception of the first column. We assume that\n')
fprintf('              the first column contains subject identifying information\n')
fprintf('              First line of the file should contain header information.\n')
fprintf('              VERY IMPORTANT: Make sure the input is normalized in the same fashion\n')
fprintf('              as data used to train model.\n')
fprintf('\n')
fprintf(' [--outputDir, -o] : directory where the output clustering assignments \n')
fprintf('will be saved. (REQUIRED)\n')
fprintf('\n')
fprintf(' [--model, -m] : HYDRA model that is previously trained. (REQUIRED)\n')
fprintf('\n')
fprintf(' OUTPUT:\n')
fprintf(' CIDX: sub-clustering assignments of the inputs.\n')
fprintf('\n')
fprintf(' NOTE: to compile this function do\n')
fprintf(' mcc -m  hydra_predict.m\n')
fprintf('\n')
fprintf('\n')
fprintf(' EXAMPLE USE (in matlab) \n');
fprintf(' hydra_predict(''-i'',''test.csv'',''-o'',''.'',''-m'',''hydra_model.mat'') \n');
fprintf(' EXAMPLE USE (in command line) \n');
fprintf(' hydra_predict -i test.csv -o . -m hydra_model.mat \n');
fprintf('======================================================================\n');
fprintf('Contact: software@cbica.upenn.edu\n');
fprintf('\n');
fprintf('Copyright (c) 2018 University of Pennsylvania. All rights reserved.\n');
fprintf('See COPYING file or http://www.med.upenn.edu/sbia/software/license.html\n');
fprintf('======================================================================\n');
end