#LOAD DATA FROM GBIF

#Anopheles gambiae (agam)
agam.gbif<-gbif(genus="Anopheles", species="gambiae*", geo=T, removeZeros=T, concept=T)
save(agam.gbif, file="agam.RData")
agam.use<-agam.gbif

#Apis mellifera --> domesticated?
amel.gbif<-gbif(genus="Apis", species="mellifera*", geo=T, removeZeros=T, concept=T)
save(amel.gbif, file="amel.RData")
amel.use<-amel.gbif

#Arabidopsis thaliana
atha.gbif<-gbif(genus="Arabidopsis", species="thaliana*", geo=T, removeZeros=T, concept=T)
save(atha.gbif, file="atha.RData")
atha.use<-atha.gbif

#Brachypodium distachyon
bdis.gbif<-gbif(genus="Brachypodium", species="distachyon*", geo=T, removeZeros=T, concept=T)
save(bdis.gbif, file="bdis.RData")
bdis.use<-bdis.gbif

#Bombyx mandarina
bman.gbif<-gbif(genus="Bombyx", species="mandarina*", geo=T, removeZeros=T, concept=T)
save(bman.gbif, file="bman.RData")
bman.use<-bman.gbif

#Bos taurus -- domesticated not sure what to use
btau.gbif<-gbif(genus="Bos", species="taurus*", geo=T, removeZeros=T, concept=T)
save(btau.gbif, file="btau.RData")
btau.use<-btau.gbif

#Caenorhabditis briggsae
cbri.gbif<-gbif(genus="Caenorhabditis", species="briggsae*", geo=T, removeZeros=T, concept=T)
save(cbri.gbif, file="cbri.RData")
cbri.use<-cbri.gbif

#Caenorhabditis elegans
cele.gbif<-gbif(genus="Caenorhabditis", species="elegans*", geo=T, removeZeros=T, concept=T)
save(cele.gbif, file="cele.RData")
cele.use<-cele.gbif

#Citrus reticulata
cret.gbif<-gbif(genus="Citrus", species="reticulata*", geo=T, removeZeros=T, concept=T)
save(cret.gbif, file="cret.RData")
cret.use<-cret.gbif

#Citrullus lanatus
clan.gbif<-gbif(genus="Citrullus", species="lanatus*", geo=T, removeZeros=T, concept=T)
save(clan.gbif, file="clan.RData")
clan.use<-clan.gbif

#Canis lupus
clup.gbif<-gbif(genus="Canis", species="lupus*", geo=T, removeZeros=T, concept=T)
save(clup.gbif, file="clup.RData")
clup.use<-clup.gbif

#Capsella rubella
crub.gbif<-gbif(genus="Capsella", species="rubella*", geo=T, removeZeros=T, concept=T)
save(crub.gbif, file="crub.RData")
crub.use<-crub.gbif

#Cucumis sativus
csat.gbif<-gbif(genus="Cucumis", species="sativus*", geo=T, removeZeros=T, concept=T)
save(csat.gbif, file="csat.RData")
csat.use<-csat.gbif

#Drosophila melanogaster
dmel.gbif<-gbif(genus="Drosophila", species="melanogaster*", geo=T, removeZeros=T, concept=T)
save(dmel.gbif, file="dmel.RData")
dmel.use<-dmel.gbif

#Drosophila pseudoobscura
dpse.gbif<-gbif(genus="Drosophila", species="pseudoobscura*", geo=T, removeZeros=T, concept=T)
save(dpse.gbif, file="dpse.RData")
dpse.use<-dpse.gbif

#Danio rerio
drer.gbif<-gbif(genus="Danio", species="rerio*", geo=T, removeZeros=T, concept=T)
save(drer.gbif, file="drer.RData")
drer.use<-drer.gbif

#Equus ferus przewalskii
efer.gbif<-gbif(genus="Equus", species="ferus*", geo=T, removeZeros=T, concept=T)
save(efer.gbif, file="efer.RData")
efer.use<-efer.gbif

#Gasterosteus aculeatus
gacu.gbif<-gbif(genus="Gasterosteus", species="aculeatus*", geo=T, removeZeros=T, concept=T)
save(gacu.gbif, file="gacu.RData")
gacu.use<-gacu.gbif

#Gallus gallus
ggal.gbif<-gbif(genus="Gallus", species="gallus*", geo=T, removeZeros=T, concept=T)
save(ggal.gbif, file="ggal.RData")
ggal.use<-ggal.gbif

#Glycine max
gmax.gbif<-gbif(genus="Glycine", species="max*", geo=T, removeZeros=T, concept=T)
save(gmax.gbif, file="gmax.RData")
gmax.use<-gmax.gbif

#Gossypium raimondii
grai.gbif<-gbif(genus="Gossypium", species="raimondii*", geo=F, removeZeros=F, concept=T)
save(grai.gbif, file="grai.RData")
grai.use<-grai.gbif

#Heliconius melpomene
hmel.gbif<-gbif(genus="Heliconius", species="melpomene*", geo=T, removeZeros=T, concept=T)
save(hmel.gbif, file="hmel.RData")
hmel.use<-hmel.gbif

#Homo sapiens
hsap.gbif<-gbif(genus="Homo", species="sapiens*", geo=T, removeZeros=T, concept=T)
save(hsap.gbif, file="hsap.RData")
hsap.use<-hsap.gbif

#Lepisosteus oculatus
locu.gbif<-gbif(genus="Lepisosteus", species="oculatus*", geo=T, removeZeros=T, concept=T)
save(locu.gbif, file="locu.RData")
locu.use<-locu.gbif

#Meleagris gallopavo
mgal.gbif<-gbif(genus="Meleagris", species="gallopavo*", geo=T, removeZeros=T, concept=T)
save(mgal.gbif, file="mgal.RData")
mgal.use<-mgal.gbif

#Macaca mulatta
mmul.gbif<-gbif(genus="Macaca", species="mulatta*", geo=T, removeZeros=T, concept=T)
save(mmul.gbif, file="mmul.RData")
mmul.use<-mmul.gbif

#Mus musculus castaneus
mcas.gbif<-gbif(genus="Mus", species="musculus*", geo=T, removeZeros=T, concept=T)
save(mcas.gbif, file="mcas.RData")
mcas.use<-mcas.gbif

#Medicago truncatula
mtru.gbif<-gbif(genus="Medicago", species="truncatula*", geo=T, removeZeros=T, concept=T)
save(mtru.gbif, file="mtru.RData")
mtru.use<-mtru.gbif

#Ovis canadensis
ocan.gbif<-gbif(genus="Ovis", species="canadensis*", geo=T, removeZeros=T, concept=T)
save(ocan.gbif, file="ocan.RData")
ocan.use<-ocan.gbif

#Oryzias latipes
olat.gbif<-gbif(genus="Oryzias", species="latipes*", geo=T, removeZeros=T, concept=T)
save(olat.gbif, file="olat.RData")
olat.use<-olat.gbif

#Oryza rufipogon
oruf.gbif<-gbif(genus="Oryza", species="rufipogon*", geo=T, removeZeros=T, concept=T)
save(oruf.gbif, file="oruf.RData")
oruf.use<-oruf.gbif

#Papio anubis
panu.gbif<-gbif(genus="Papio", species="anubis*", geo=T, removeZeros=T, concept=T)
save(panu.gbif, file="panu.RData")
panu.use<-panu.gbif

#Prunus davidiana
pdav.gbif<-gbif(genus="Prunus", species="davidiana*", geo=F, removeZeros=T, concept=T)
save(pdav.gbif, file="pdav.RData")
pdav.use<-pdav.gbif

#Populus trichocarpa
ptri.gbif<-gbif(genus="Populus", species="trichocarpa*", geo=T, removeZeros=T, concept=T)
save(ptri.gbif, file="ptri.RData")
ptri.use<-ptri.gbif

#Sorghum bicolor -> domesticated
sbic.gbif<-gbif(genus="Sorghum", species="bicolor*", geo=T, removeZeros=T, concept=T)
save(sbic.gbif, file="sbic.RData")
sbic.use<-sbic.gbif

#Setaria italica -> domesticated
sita.gbif<-gbif(genus="Setaria", species="italica*", geo=T, removeZeros=T, concept=T)
save(sita.gbif, file="sita.RData")
sita.use<-sita.gbif

#Sus scrofa
sscr.gbif<-gbif(genus="Sus", species="scrofa*", geo=T, removeZeros=T, concept=T)
save(sscr.gbif, file="sscr.RData")
sscr.use<-sscr.gbif

#Zea mays ssp parviglumis
zmay.gbif<-gbif(genus="Zea", species="mays*", geo=T, removeZeros=T, concept=T)
save(zmay.gbif, file="zmay.RData")
zmay.use<-zmay.gbif

#Ficedula albicollis
falb.gbif<-gbif(genus="Ficedula", species="albicollis*", geo=T, removeZeros=T, concept=T)
save(falb.gbif, file="falb.RData")
falb.use<-falb.gbif

#Cynoglossus semilaevis
csem.gbif<-gbif(genus="Cynoglossus", species="semilaevis*", geo=F, removeZeros=F, concept=T)
save(csem.gbif, file="csem.RData")
csem.use<-csem.gbif

