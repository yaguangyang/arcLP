function infeasible_method(algorithm,A,b,c);
% infeasible_method('arcLP',A,b,c) calls arc-search algorithm
% infeasible_method('mehrotra',A,b,c) calls Mehrotra algorithm
if strcmp(algorithm,'arcLP')
        [x,obj,kk,infe]=arcLP(A,b,c,1);
else
        [x,obj,kk,infe]=mehrotra(A,b,c,1);
end
disp(['iter: ', num2str(kk), ',   obj: ', num2str(obj), ...
    ',   infe: ', num2str(infe)]);
