function error_table(h,k,E)
%
% Print out table of errors, ratios, and observed order of accuracy.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

ntest = length(h);
ratio = nan(size(h));   % initialize storage
order = nan(size(h));   % initialize storage

ntestk = length(k);
ratiok = nan(size(k));   % initialize storage
orderk = nan(size(k));

for j=2:ntest
   ratio(j) = E(j-1)/E(j);
   ratiok(j) = E(j-1)/E(j);
   order(j) = log(abs(ratio(j))) / log(abs(h(j-1)/h(j)));
   orderk(j) = log(abs(ratiok(j))) / log(abs(k(j-1)/k(j)));
end


% print out table:

disp(' ')
disp('      h        error       ratio       observed order')
for j=1:ntest
   %fileID = fopen('Error.txt','w');
   fprintf(' %9.5f  %12.5e %9.5f %15.5f\n',h(j),E(j),ratio(j),order(j));
end
%fclose(fileID);
disp(' ')

disp(' ')
disp('      k        error       ratio       observed order')
for i=1:ntestk
   fprintf(' %9.5f  %12.5e %9.5f %15.5f\n',k(i),E(i),ratiok(i),orderk(i));
end
disp(' ')

