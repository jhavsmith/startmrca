
cprobback.func <- function(newhap.vec, hapdata.M, hapfreqs.vec,rhobar, lambda.vec, pos.vec, prob, prob.vec, copiedh.vec, seed) {
#        library.dynam("cprobback",package="startmrca")
      	return(.C("cprobback",
      	newhap    = as.integer(newhap.vec),
      	hapdata   = as.integer(hapdata.M),
      	nloci     = as.integer(dim(hapdata.M)[2]),
      	nhaps     = as.integer(dim(hapdata.M)[1]),
      	hapfreqs  = as.integer(hapfreqs.vec ),
      	nchroms   = as.integer(sum(hapfreqs.vec)),
      	meanrho   = as.double(rhobar),
      	lambda    = as.double(lambda.vec), # length nloci - 1
      	positions = as.double(pos.vec),
      	pihat     = as.double(prob),
      	pihat.vec = as.double(prob.vec),
      	copiedh   = as.integer(copiedh.vec),
      	seed      = as.integer(seed))) 
}



