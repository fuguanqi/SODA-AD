% function Data= datainput_cplex
% Data.objfunction=@(x)cplex(x); %objective function handle
% end %function


function y=datainput_fun(x) %objective function
xx=x(:)'; % make sure vector is row vector
load (strcat('problems\Data1.mat'));

...


y=obj;


% sampledata = [sampledata; x ,y, t]; %collect sample data (point x, value y, time t)
end %
