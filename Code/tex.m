disp("This is a simple test")

% uncomment to erase log
% by default diary appends new data
%delete("matlabDiary.txt")

diary("matlabDiary.txt")
diary on

x = 5;
y = 3;

disp(["x + y = ",num2str(x+y)])

diary off


disp(["x - y = ",num2str(x-y)])

diary on


disp(["x * y = ",num2str(x*y)])

diary off