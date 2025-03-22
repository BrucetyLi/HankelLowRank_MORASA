function value = mat_compare(new,old)
value = norm(new-old,"fro")/norm(new,'fro');
end