# Universal stuff that needs to be run first
methods::setOldClass("UERE")
new.UERE <- methods::setClass("UERE",contains="list",representation=methods::representation(info="list"),
                              prototype=methods::prototype(list(UERE=cbind(0),DOF=cbind(0),AICc=Inf,Zsq=Inf,VAR.Zsq=Inf,N=0),info=list()))

methods::setOldClass("variogram")
new.variogram <- methods::setClass("variogram",representation=methods::representation("data.frame",info="list",UERE="UERE"),
                                   prototype=methods::prototype(data.frame(),info=list(),UERE=new.UERE()) )
