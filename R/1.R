# Universal stuff that needs to be run first
methods::setOldClass("variogram")
new.variogram <- methods::setClass("variogram",representation=methods::representation("data.frame",info="list",UERE="UERE"),
                                   prototype=methods::prototype(data.frame(),info=list(),UERE=new.UERE()) )
