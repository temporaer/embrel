function print_fig(fname)

eval(sprintf('print -dpng %s.png',fname));
eval(sprintf('print -depsc %s.eps',fname));
eval(sprintf('hgsave  %s',fname));