# Universal stuff that needs to be run first

methods::setOldClass("variogram")
new.variogram <- methods::setClass("variogram",representation(Distance="numeric",
                                                              Gamma="numeric",
                                                              CI_min="numeric",
                                                              CI_max="numeric"))
