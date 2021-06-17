function figurep(base, profile, name)

print('-dpng','-r300',sprintf('figures/%s_profile_%d_%s.png',base,profile,name));
print('-depsc',sprintf('figures/%s_profile_%d_%s.eps',base,profile,name));
print('-dsvg',sprintf('figures/%s_profile_%d_%s.svg',base,profile,name));
fprintf('Wrote figures/%s_profile_%d_%s.{png,eps,svg}\n',base,profile,name)
