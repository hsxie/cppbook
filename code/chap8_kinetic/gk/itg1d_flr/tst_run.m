close all; clear; clc;

for prun=1:3
    run gkpic1d_itg_flr;
    print(gcf,'-dpng',[str,'_',num2str(prun),'_history.png']);
end