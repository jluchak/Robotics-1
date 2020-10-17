%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File name Jacobian2.m
% This function is to calculate the Jacobian matrix, and will be called
%  in the script Tutorial_2_question_2
function J=Jacobian2(q_l,q_r)
l=0.2;
J=l*[-sin(q_l) -sin(q_r); cos(q_l) cos(q_r)];
end
% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
