function [s_emp_recession] = s_emp_recession(s_bar, s)

t = length(s);
s_emp_recession = nan(t, 1);

% Creating recession dummey
for i = 1:t
   if s(i) < s_bar
       s_emp_recession(i) = 1;
   else 
       s_emp_recession(i) = 0;
   end
end

end