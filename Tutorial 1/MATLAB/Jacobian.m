% File name Jacobian.m
% This function is to calculate the Jacobian matrix, and will be called
%  in the script Tutorial_1_question_2
function J=Jacobian(ql,qr,theta)
a=0.3;
b=0.3;
c=0.2;
J=[-a*sin(ql) -b*sin(qr) -c*sin(theta); a*cos(ql) b*cos(qr) c*cos(theta)];
end
% End of function