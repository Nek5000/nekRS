
function myprint(format, fileName, option)

fig = gcf;
set(fig,'PaperPositionMode','auto');
fig_pos = get(fig,'PaperPosition');
set(fig,'PaperSize',[fig_pos(3) fig_pos(4)]);
 
print(format, fileName);
