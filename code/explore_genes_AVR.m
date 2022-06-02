%explore paralogs in base and in c. intermedia models

cint = load('..\models\candida_intermedia\cintGEM_oxido_AVR.mat');
base = load('..\models\yeastGEM.mat');
model_leloir = load('..\models\candida_intermedia\cint_leloir_AVR.mat');

%starting with the GAL cluster

%Just find Gal1 (YBR020W)
cint.model.genes(find(contains(cint.model.genes,'YBR020W'))) % found two
base.model.genes(find(contains(base.model.genes,'YBR020W'))) % found one
model_leloir.model_leloir.genes(find(contains(model_leloir.model_leloir.genes,'YBR020W'))) % found two

%so just find Gal3 (YDR009W)
cint.model.genes(find(contains(cint.model.genes,'YDR009W'))) %didn't find it
base.model.genes(find(contains(base.model.genes,'YDR009W'))) %found it
model_leloir.model_leloir.genes(find(contains(model_leloir.model_leloir.genes,'YDR009W'))) %didn't find it, so Gal3 got clumped up with Gal1

%Just find Gal7 (YBR018C) 
cint.model.genes(find(contains(cint.model.genes,'YBR018C')))%found it
base.model.genes(find(contains(base.model.genes,'YBR018C')))%found it
model_leloir.model_leloir.genes(find(contains(model_leloir.model_leloir.genes,'YBR018C')))%found it

%Just find Gal10 (YBR018C) %found 1
model.genes(find(contains(model.genes,'YBR018C')))
