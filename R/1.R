# Universal stuff that needs to be run first
methods::setOldClass("variogram")

new.variogram <- methods::setClass("variogram",representation=methods::representation("data.frame",info="list"),
                                   prototype=methods::prototype(data.frame(),info=list()) )
