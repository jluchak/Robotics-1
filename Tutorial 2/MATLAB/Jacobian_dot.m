%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File name Jacobian_dot.m
% This function is to calculate the Jacobian_dot matrix, and will be called
%  in the script Tutorial_2_question_2
function J_dot=Jacobian_dot(q_l,q_r,q_l_dot,q_r_dot)
l=0.2;
J_dot = [-l*cos(q_l)*q_l_dot, -l*cos(q_r)*q_r_dot;
         -l*sin(q_l)*q_l_dot, -l*sin(q_r)*q_r_dot];
end
% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
