x = copy(:,29);
[phat,pci] = mle(x,'pdf',@(x,a,b,c,d) mlPDFfitting(x,a,b,c,d),'Start',[0,0,0,5])

