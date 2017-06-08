function ax=PlotExSpec(avg1,chamchi,chameps,iwind)
%%
%ax=figure(1);clf
ax=gca;

h1=loglog(avg1.ks(iwind,:),avg1.kspec(iwind,:),'ko-');
hold on
hchifit=loglog(avg1.kks(iwind,:),avg1.kkspec(iwind,:),'k--');

% compute kraichnan spec for actual value measured by chameleon also
clear nu tdif k chi eps kb spec_vals
nu=avg1.nu(iwind);
tdif=avg1.tdif(iwind);
qq=7;
k=avg1.kks(iwind,:);
chi=chamchi(iwind);
eps=chameps(iwind);
kb = (((eps./(nu.^3)).^.25 )/2/pi).*sqrt(nu./tdif);
[spec_vals]=kraichnan(nu,k,kb,tdif,chi,qq);
hcham=loglog(k,spec_vals,'m','linewidth',2);
%
ylim([1e-12 1e0])
grid on
xlabel('wavenumber [cpm]','fontsize',16)
ylabel('\Phi_{dT/dz}[K^2m^{-2}cpm^{-1}]','fontsize',16)
xlim([0.1 1.5e3])

% mark kmax
freqline(avg1.Params.fmax./avg1.fspd(iwind),'r--');
text(avg1.Params.fmax./avg1.fspd(iwind),1e-2,'k_{max}','fontsize',15);

% mark batchelor wavenumber k_b
hkb=freqline(kb,'m--');
set(hkb,'linewidth',3);
text(kb,1e-2,'k_b','fontsize',15);

title(['P=' sprintf('%.2f',avg1.P(iwind)) 'db, Chameleon \epsilon=' sprintf('%.2e',chameps(iwind)) ])
%legend([h1 h2 h3],['fc=' num2str(avg1.Params.fc)],['fc=' num2str(avg2.Params.fc)],['fc=' num2str(avg3.Params.fc)])
% legend([h1 h2 hcham],'\chi pod corr.','\chi pod no corr.','chameleon')
legend([h1 hchifit hcham],'obs.','\chi pod fit ','Chameleon','location','northwest')

freqline(avg1.kstart(iwind),'b--');
%freqline(avg1.kstop(iwind),'b--')

fs=14;
text(1e-1,1e-8,['\chi_{\chi} =' sprintf('%.2e',avg1.chi1(iwind))],'fontweight','bold','fontsize',fs)
text(1e-1,1e-9,['\chi_{\epsilon} =' sprintf('%.2e',chamchi(iwind))],'color','m','fontweight','bold','fontsize',fs)
text(1e-1,1e-10,['\chi_{\chi}/\chi_{\epsilon} =' sprintf('%.2f',avg1.chi1(iwind)/chamchi(iwind))],'color','b','fontweight','bold','fontsize',fs)

text(1e0,1e-8,['\epsilon_{\chi} =' sprintf('%.2e',avg1.eps1(iwind))],'fontweight','bold','fontsize',fs)
text(1e0,1e-9,['\epsilon_{\epsilon} =' sprintf('%.2e',chameps(iwind))],'color','m','fontweight','bold','fontsize',fs)
text(1e0,1e-10,['\epsilon_{\chi}/\epsilon_{\epsilon} =' sprintf('%.2f',avg1.eps1(iwind)/chameps(iwind))],'color','b','fontweight','bold','fontsize',fs)

%%