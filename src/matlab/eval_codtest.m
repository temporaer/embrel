try
	rv = evalc('codtest');
	fprintf('-----SUCCESS\n%s', rv);
catch
	fprintf('-----ERROR\n%s\n%s', lasterr);
end
quit;
