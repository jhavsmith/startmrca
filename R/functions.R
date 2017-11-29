#This runs the full sequence of functions below.
run.startmrca = function(vcf.file, rec.file, mut.rate, rec.rate = NULL, nsel = NULL, nanc = NULL, chain.length = 15000, proposal.sd = 20, nanc.post = 100, sample.ids, refsample.ids = NULL, pos, sel.allele = 1, bed.file = NULL, upper.t.limit = 2000) {
        params = make.params(vcf.file, rec.file, mut.rate, rec.rate, nsel, nanc, chain.length, proposal.sd, nanc.post, sample.ids, refsample.ids, pos, sel.allele, bed.file, upper.t.limit)
        input.list = get.vcfdata.func(params)
        prep.list = prep.func(input.list,params)
        anc.list = anc.func(prep.list)
        new.recmap = rec.func(params,prep.list)
        mcmc.list = one.list.func(prep.list,anc.list,new.recmap)
        mcmc.list = mcmc.init.func(mcmc.list)
        mcmc.output = startchain.func(mcmc.list,params)
        save(mcmc.output,file = paste(params$file.name,"_mcmc_list.RDATA",sep=''))
}

#This puts the parameter values into list to pass through the following functions.
make.params = function(vcf.file, rec.file, mut.rate, rec.rate, nsel, nanc, chain.length, proposal.sd, nanc.post, sample.ids, refsample.ids, pos, sel.allele, bed.file, upper.t.limit) {
    file.name = strsplit(vcf.file, ".vcf.gz")
    params = list("vcf.file"      = vcf.file,     "rec.file"      = rec.file,
                  "rec.rate"      = rec.rate,     "mut.rate"      = mut.rate,
                  "nsel"          = nsel,         "nanc"          = nanc,
                  "chain.length"  = chain.length, "proposal.sd"   = proposal.sd,
                  "nanc.post"     = nanc.post,    "sample.ids"    = sample.ids,
                  "refsample.ids" = refsample.ids,"pos"           = pos,
                  "sel.allele"    = sel.allele,   "file.name"     = file.name,
                  "bed.file"      = bed.file,     "upper.t.limit" = upper.t.limit)
    if (chain.length < nanc.post) { stop("nanc.post must be less than chain.length")}
    return(params)
}

#This puts the genotype data from the the vcf into an R matrix format.
get.vcfdata.func = function(params) {
    print("Getting data from the vcf.")
    vcf.file = params$vcf.file
    file.name = params$file.name
    sel.pos = params$pos
    # pulls the sample labels from the vcf
    system(paste('zcat < ',vcf.file,' | grep CHROM -m 1 > ', paste(file.name,'_fields.txt',sep='')))
	fields = scan(paste(file.name,'_fields.txt',sep=''),what='character')
	system(paste("rm", paste(file.name,'_fields.txt',sep='')))
	sample.fields = fields[-c(1:9)]
	ref.fields = "NA"
	if (!is.character(params$sample.ids)) {
	    # If no carriers are specified it uses the full sample.
	    print("No selected allele carrier IDs are specified.")
	    print("Using carriers from the full vcf.")
	    selid = sample.fields
	    sel.fields = c(1:length(selid))
	}
	if (is.character(params$sample.ids)) {
	    # If carriers are specified it identifies which columns to use in the vcf.
		if (!file.exists(params$sample.ids)) {
		    # Checks to see that the sample.id file was specified correctly.
	       print(paste(params$sample.ids,"doesn't exist.",sep=' '))
	    } else {
	       print("Selected allele carrier IDs are specified.")
	       selid = scan(params$sample.ids,what="character")
	       sel.fields = which(sample.fields%in%selid==TRUE)
	       if (length(sel.fields)==0) {
	          # The ids in the sample.id file may not be represented in the vcf.
              print("Selected allele carrier IDs are not in this vcf.")
              print("Using any carriers in the full vcf instead.")
           }
        }
	}
    if (!is.character(params$refsample.ids)) {
        # if no reference individuals are specified, it uses the full sample.
        print("No reference panel IDs are specified.")
	    print("Using non-carriers from the full vcf.")
	    refid = sample.fields
	    ref.fields = which(sample.fields%in%sample.fields[-sel.fields]==TRUE)
	}
	if (is.character(params$refsample.ids)) {
	    if (!file.exists(params$refsample.ids)) {
	       print(paste(params$refsample.ids,"doesn't exist.",sep=' '))
	       break
	    } else {
	       print("Reference panel IDs are specified.")
	       # This pulls out the reference individuals from the full sample.
	       if (length(which(params$refsample.ids%in%params$sample.ids==TRUE))==length(params$sample.ids)) {
	          print("Reference panel is the same as the carrier panel.")
	          refid = selid
	          ref.fields = sel.fields
	       } else {
              refid = scan(params$refsample.ids,what="character")
              ref.fields = which(sample.fields%in%refid==TRUE)
           }
           if (length(ref.fields)==0) {
              print("Reference panel IDs are not in this vcf.")
              print("Using any non-carriers in the full vcf instead.")
           }
        }
    }
    # If there's redundancy between the carrier and reference panels it removes duplicates.
	pop.fields        = c(c(1:9),c(sel.fields,ref.fields) + 9)
	if (length(which(sel.fields%in%ref.fields==TRUE))==length(sel.fields)) {
        pop.fields = c(c(1:9),c(sel.fields) + 9)
	}
	# This trims the vcf file down to the sample IDs we want using the columns from pop.fields.
    system(paste('zcat < ',vcf.file, '| cut -f',paste(pop.fields,collapse=","),'> ',paste(file.name,".txt",sep='')))
	vcf.sample        = read.table(paste(file.name,".txt",sep=''))
	system(paste("rm", paste(file.name,".txt",sep='')))
	chrom             = as.numeric(paste(unlist(vcf.sample[1,1])))
	# We want to remove any sites that aren't biallelic.
    alt.allele.number = nchar(paste(unlist(vcf.sample[,5])))
    ref.allele.number = nchar(paste(unlist(vcf.sample[,4])))
    alt.remove        = which(alt.allele.number!=1)
    ref.remove        = which(ref.allele.number!=1)
    remove.alleles    = sort(unique(c(alt.remove,ref.remove)))
    if (!identical(remove.alleles,integer(0))) {
        vcf.sample    = vcf.sample[-remove.alleles,]
    }
    options(scipen=0)
    # This is where we make our matrix of genotypes from the vcf.
    genotype.matrix   = matrix(nrow=(ncol(vcf.sample)-9)*2,ncol=nrow(vcf.sample))
    odds              = seq(1,nrow(genotype.matrix),2)
    pb                = txtProgressBar(1,length(odds),1,style=3)
    for (i in 1:length(odds)) {
        index = odds[i]
        genotype.matrix[index,]     = as.numeric(substr(vcf.sample[,-c(1:9)][,i],1,1))
        genotype.matrix[(index+1),] = as.numeric(substr(vcf.sample[,-c(1:9)][,i],3,3))
        setTxtProgressBar(pb, i)
    }
    close(pb)
    whole.sample      = genotype.matrix
    # Getting allele frequencies.
    allele.freqs      = c(1:ncol(whole.sample))
    for (i in 1:ncol(whole.sample)) {
        allele.freqs[i] = sum(whole.sample[,i])/nrow(whole.sample)
    }
    # We want to remove any fixed sites in our sample.
    invariant.sites   = c(which(allele.freqs==1),which(allele.freqs==0))
    if (!identical(invariant.sites,integer(0))) {
        vcf.sample    = vcf.sample[-invariant.sites,]
        allele.freqs  = allele.freqs[-invariant.sites]
        whole.sample  = whole.sample[,-invariant.sites]
    }
    positions         = as.numeric(paste(unlist(vcf.sample[,2])))
    fields.list       = list("sel.fields"=sel.fields,"ref.fields"=ref.fields)
    input.list = list("whole.sample" = whole.sample, "positions"  = positions,
                      "allele.freqs" = allele.freqs, "const.mult" = params$const.mult,
                      "fields.list"  = fields.list,  "chr"        = chrom)
    return(input.list)
}

# This is another data processing function before running the MCMC.
prep.func = function(input.list,params) {
    # Here we define the variables we want from the lists 'input.list' and 'params'.
    psel          = params$pos
    nanc          = params$nanc
    nsel          = params$nsel
    mut.rate      = params$mut.rate
    sel.allele    = params$sel.allele
    chain.length  = params$chain.length
    chr           = input.list$chr
    whole.sample  = input.list$whole.sample
    positions     = input.list$positions
    all.positions = cumsum(c(0,diff(positions)))+1
    allele.freqs  = input.list$allele.freqs
    sel.site      = which(positions==psel)
    bed.file      = params$bed.file
    # Checking to make sure the selected site is in the vcf.
    if (length(sel.site)==0) {
        print("The selected site is not a position in the VCF file.")
    }
    fields.list   = input.list$fields.list
    sel.fields    = fields.list$sel.fields
    ref.fields    = fields.list$ref.fields
    sel.fields    = rbind(sel.fields,sel.fields)[1:(length(sel.fields)*2)]
    sorted.fields = sel.fields
    ref.sample    = "NA"
    # If the reference panel was specified, then we put those haplotypes into their own matrix.
    if (!is.na(ref.fields[1])) {
        # In case the reference panel columns are the same as the carrier panel.
        if (length(which(fields.list$sel.fields%in%ref.fields==TRUE))==length(ref.fields)) {
           ref.cols = which(sorted.fields%in%sel.fields)
           ref.sample    = whole.sample[ref.cols,]
        # Otherwise we put them in ref.sample.
        } else {
           ref.fields    = rbind(ref.fields,ref.fields)[1:(length(ref.fields)*2)]
           sorted.fields = sort(c(sel.fields,ref.fields))
           ref.cols      = which(sorted.fields%in%ref.fields)
           ref.sample    = whole.sample[ref.cols,]
        }
    }
    # Now we put the carrier panel into its own matrix.
    sel.cols      = which(sorted.fields%in%sel.fields)
    sel.sample    = whole.sample[sel.cols,]
    # Here we figure out which individuals have the selected allele in the carrier panel.
    if (nsel=="NA") {
        # If the sample size for the selected panel wasn't specified.
        sel.haps = which(sel.sample[,sel.site]==sel.allele)
        sample   = sel.sample[sel.haps,]
    }
    if (nsel!="NA") {
        # If the sample size for the selected panel was specified, then we randomly select that many.
        sel.haps = which(sel.sample[,sel.site]==sel.allele)
        if (length(sel.haps)>=nsel) {
            sample = sel.sample[sample(sel.haps,nsel),]
        }
        # In case there are fewer haplotypes than specified by nsel.
        if (length(sel.haps)<nsel) {
            sample = sel.sample[sel.haps,]
        }
    }
    # Same deal as above, but for the reference panel.
    if (ref.sample[1]=="NA") {
        if (nanc=="NA") {
            ancestral.haps = which(sel.sample[,sel.site]==abs(sel.allele-1))
            cont.sample    = sel.sample[ancestral.haps,]
        }
        if (nanc!="NA") {
            ancestral.haps = which(sel.sample[,sel.site]==abs(sel.allele-1))
            if (length(ancestral.haps)>=nanc) {
                cont.sample = sel.sample[c(1:nanc),]
            }
            if (length(ancestral.haps)<nanc) {
                cont.sample = sel.sample[ancestral.haps,]
            }
        }
    }
    if (ref.sample[1]!="NA") {
        if (nanc=="NA") {
            ancestral.haps = which(ref.sample[,sel.site]==abs(sel.allele-1))
            cont.sample    = ref.sample[ancestral.haps,]
        }
        if (nanc!="NA") {
            ancestral.haps = which(ref.sample[,sel.site]==abs(sel.allele-1))
            if (length(ancestral.haps)>=nanc) {
                cont.sample = ref.sample[sample(ancestral.haps,nanc),]
            }
            if (length(ancestral.haps)<nanc) {
                cont.sample = ref.sample[ancestral.haps,]
            }
        }
    }
    # Checks to make sure there are individuals in the reference panel without the selected allele.
    if (length(ancestral.haps)==0) {
        print("There are no individuals for the reference panel in this VCF.")
    }
    # Now we divide the data up into halves. The model runs from left to right starting
    # from the selected site, so we have to flip the left side around.
    left.sample  = sample[,sel.site:1]
    left.cont    = cont.sample[,sel.site:1]
    left.pos     = abs(positions[sel.site:1]-positions[sel.site])+1
    right.sample = sample[,sel.site:ncol(sample)]
    right.cont   = cont.sample[,sel.site:ncol(cont.sample)]
    right.pos    = (positions[sel.site:ncol(sample)]-positions[sel.site])+1
    # The bedfile can tell us which (if any) positions were not actually sequenced in the vcf.
    # We need to know this because invariant sites may instead just be unobserved (or unsequenced) sites.
    if (typeof(bed.file)=="character") {
        bed.list = bed.file.func(params$bed.file,psel,left.pos,right.pos)
    }
    if (typeof(bed.file)!="character") {
        left.gaps = NULL
        right.gaps = NULL
        bed.list = list("left.gaps"=left.gaps,"right.gaps"=right.gaps)
    }
    # this list gets passed to the next function(s).
    prep.list = list("left.sample"   = left.sample,  "right.sample"  = right.sample,
                     "left.cont"     = left.cont,    "right.cont"    = right.cont,
                     "left.pos"      = left.pos,     "right.pos"     = right.pos,
                     "sample"        = sample,       "cont.sample"   = cont.sample,
                     "all.positions" = all.positions,"positions"     = positions,
                     "sel.site"      = sel.site,     "mut.rate"      = mut.rate,
                     "chr"           = chr,          "chain.length"  = chain.length,
                     "bed.list"      = bed.list,     "upper.t.limit" = params$upper.t.limit)
   return(prep.list)
}

# This function uses a bed file to correct for any gaps that may be in the sequence data.
bed.file.func = function(bed.file,psel,left.pos,right.pos) {
    new.bed.file = no.list.func(read.table(bed.file),"num")
    chroms = length(unique(new.bed.file[,1]))
    if (chroms!=1) {
       print("The bed file has too many chromosomes")
       break
    }
    print("Bed file is specified.")
    segs = new.bed.file[-1,2]-new.bed.file[-nrow(new.bed.file),3]
    seg.gaps.b = c(0,which(segs!=1))
    seg.gaps.a = seg.gaps.b+1
    seg.gaps.b = c(seg.gaps.b[-1],(length(segs)+1))
    compressed.bed = matrix(nrow=length(seg.gaps.a),ncol=2)
    for (i in 1:length(seg.gaps.a)) {
        compressed.bed[i,] = c(new.bed.file[seg.gaps.a[i],2],new.bed.file[seg.gaps.b[i],3])
    }
    left.segs  = compressed.bed[which(compressed.bed[,1]<psel),]
    left.segs = cbind(rev(psel-left.segs[,2]),rev(psel-left.segs[,1]))
    right.segs = compressed.bed[which(compressed.bed[,2]>psel),]
    right.segs = cbind(right.segs[,1]-psel,right.segs[,2]-psel)
    left.segs[1,1] = 1
    right.segs[1,1] = 1
    left.gaps = cbind(left.segs[-nrow(left.segs),2]+1,left.segs[-1,1]-1)
    right.gaps = cbind(right.segs[-nrow(right.segs),2]+1,right.segs[-1,1]-1)
    bed.list = list("left.gaps"=left.gaps,"right.gaps"=right.gaps)
    return(bed.list)
}

# This initializes the ancestral haplotype using an algorithm described in the Appendix of the paper.
anc.func = function(prep.list) {
    print("Estimating the ancestral haplotype.")
    # The left and right sides of the selected allele are done seperately.
    for (k in 1:2) {
        # Defines some relevant objects.
        sample               = prep.list[[k]]
        cont.sample          = prep.list[[2+k]]
        anc.table            = sample
        pos                  = prep.list[[4+k]]
        anc.zeros            = which(anc.table==0)
        anc.table[anc.zeros] = 1
        hap.init             = rep(2,ncol(anc.table))
        sample.freqs         = as.numeric(apply(sample,2,function(x) sum(x)/length(x)))
        cont.sample.freqs    = as.numeric(apply(cont.sample,2,function(x) sum(x)/length(x)))
        half.freq.sites      = which(sample.freqs==0.5)
        singletons           = which(sample.freqs==(1/nrow(sample)))
        removed.sites        = c(half.freq.sites,singletons)
        # Removes singletons and sites with frequency 0.5.
        if (length(removed.sites)!=0) {
            temp.sample          = sample[,-removed.sites]
            temp.freqs           = sample.freqs[-removed.sites]
            temp.cont.freqs      = cont.sample.freqs[-removed.sites]
            temp.anc.table       = anc.table[,-removed.sites]
            temp.hap.init        = hap.init[-removed.sites]
            temp.pos             = pos[-removed.sites]
        }
        if (length(removed.sites)==0) {
            temp.sample          = sample
            temp.freqs           = sample.freqs
            temp.cont.freqs      = cont.sample.freqs
            temp.anc.table       = anc.table
            temp.hap.init        = hap.init
            temp.pos             = pos
        }
        # identifies the major and minor allele.
        for (i in 2:ncol(temp.sample)) {
            major.allele = round(mean(temp.sample[which(temp.anc.table[,(i-1)]==1),i]))
            minor.allele = abs(1-major.allele)
            if (length(which(temp.sample[,i]==minor.allele))==0) {
               next
            }
            if (length(which(temp.sample[,i]==minor.allele))!=0) {
               minor.inds = which(temp.sample[,i]==minor.allele)
               temp.anc.table[minor.inds,c(i:ncol(temp.sample))] = 0
            }
        }
        order.anc.table     = rbind(temp.pos,temp.anc.table)
        # Adds the removed sites back in from before
        if (length(removed.sites)!=0) {
            if (length(removed.sites)==1) {
                add.anc.table.sites = c(pos[removed.sites],anc.table[,removed.sites])
                order.anc.table = cbind(order.anc.table,add.anc.table.sites)
            }
            if (length(removed.sites)!=1) {
                add.anc.table.sites = rbind(pos[removed.sites],anc.table[,removed.sites])
                order.anc.table = cbind(order.anc.table,add.anc.table.sites)
            }
            new.anc.table = order.anc.table[-1,order(order.anc.table[1,])]
        }
        if (length(removed.sites)==0) {
            add.anc.table.sites = rbind(pos,anc.table)
            order.anc.table     = cbind(order.anc.table,add.anc.table.sites)
            new.anc.table       = order.anc.table[-1,order(order.anc.table[1,])]
        }
        # identifies the first "non-ancestral" site in each chromosome.
        for (j in 1:nrow(new.anc.table)) {
            if (length(which(new.anc.table[j,]==0))==0) {
               next
            }
            if (min(which(new.anc.table[j,]==0))>ncol(anc.table)) {
               next
            }
            else {
               anc.table[j,min(which(new.anc.table[j,]==0)):ncol(anc.table)] = 0
            }
        }
        # Uses the major allele among "ancestral" sites to define the ancestral haplotype.
        for (i in 1:ncol(sample)) {
            if (length(which(anc.table[,i]==1))!=0) {
              hap.init[i] = round(mean(sample[which(anc.table[,i]==1),i]))
            } else {
               hap.init[i] = rbinom(1,1,sample.freqs[i])
            }
        }
        # Re-formats for the output.
        if (k==1) {
           left.hap.init  = rev(hap.init)
           left.anc.table = apply(anc.table,1,rev)
        }
        if (k==2) {
           full.hap.init  = c(left.hap.init,hap.init[-1])
           full.anc.table = cbind(t(left.anc.table),anc.table[,-1])
        }
    }
    anc.list = list("full.hap.init"=full.hap.init,"full.anc.table"=full.anc.table)
    return(anc.list)
}

# Takes a recombination map and calculates recombination rates between each SNP in the VCF.
rec.func = function(params,prep.list) {
    print("Getting the recombination map.")
    # Defines a few objects.
    recmap.file = params$rec.file
    r           = params$rec.rate
    chrom       = prep.list$chr
    positions   = prep.list$positions
    # If no recombination map file is specified, it uses a uniform recombination rate defined by rec.rate.
    if (!is.character(recmap.file)) {
	    windows      = cbind(positions[-length(positions)],positions[-1])
        r.vec        = rep(r,nrow(windows))
    }
    if (is.character(recmap.file)) {
        # Defines the new recombination intervals from the VCF positions.
        rec.map         = no.list.func(read.table(recmap.file),"num")
        recmap.chr      = rec.map[which(rec.map[,1]==chrom),]
        recmap.windows  = cbind((recmap.chr[-nrow(recmap.chr),2]),(recmap.chr[-1,2]))
	    windows         = cbind(positions[-length(positions)],positions[-1])
        new.recs        = rep(0,nrow(windows))
	    window.size     = recmap.windows[,2]-recmap.windows[,1]
	    pb              = txtProgressBar(1,nrow(windows),1,style=3)
	    # Gets the recombination rate for each SNP window.
	    for (j in 1:nrow(windows)) {
	        left.side  = max(which(recmap.windows[,1]<=windows[j,1]))
	        if (length(which(recmap.windows[,2]>=windows[j,2]))==0) {
	            print("The recombination map is shorter than the VCF.")
	            break
	        } else {
	            right.side = min(which(recmap.windows[,2]>=windows[j,2]))
	            if (left.side==right.side) {
	                new.recs[j] = recmap.chr[left.side,3]
	            }
	            if (right.side==left.side+1) {
	                left.fraction  = (recmap.windows[left.side,2]-windows[j,1])/as.numeric(window.size[left.side])
		            right.fraction = 1-((recmap.windows[right.side,2]-windows[j,2])/as.numeric(window.size[right.side]))
		            left.recs      = left.fraction*recmap.chr[left.side,3]
		            right.recs     = right.fraction*recmap.chr[right.side,3]
		            new.recs[j]    = left.recs+right.recs
	            }
	        }
            setTxtProgressBar(pb, j)
        }
	    close(pb)
	    r.vec = new.recs
    }
	return(r.vec)
}

# This does some reorganizing. Gets things ready for the MCMC.
one.list.func = function(prep.list,anc.list,r.vec) {
    sample       = prep.list$sample
    cont.sample  = prep.list$cont.sample
    sel.site     = prep.list$sel.site
    chain.length = prep.list$chain.length
    full.anc     = anc.list[[1]]
    full.haps    = anc.list[[2]]
    mut.rate     = prep.list$mut.rate
    r            = mean(r.vec)
    left.r       = r.vec[sel.site:1]
    right.r      = r.vec[sel.site:(ncol(cont.sample))]
    t.chain      = c(0,0,0,0)
    anchap.post  = 0
    mcmc.list    = list("left.sample"  = prep.list[[1]], "right.sample" = prep.list[[2]],
                       "left.cont"     = prep.list[[3]], "right.cont"   = prep.list[[4]],
                       "left.pos"      = prep.list[[5]], "right.pos"    = prep.list[[6]],
                       "full.anc"      = full.anc,       "full.haps"    = full.haps,
                       "sample"        = prep.list[[7]], "cont.sample"  = prep.list[[8]],
                       "all.positions" = prep.list[[9]], "positions"    = prep.list[[10]],
                       "left.r"        = left.r,         "right.r"      = right.r,
                       "t.chain"       = t.chain,        "anchap.post"  = anchap.post,
                       "sel.site"      = sel.site,        "mu"           = mut.rate,
                       "r"             = r,               "bed.list"     = prep.list$bed.list,
                       "tmrca"         = prep.list$tmrca, "chain.length" = chain.length,
                       "upper.t.limit" = prep.list$upper.t.limit)
    return(mcmc.list)
}

# Computes the transition probabilities.
transition.probs.func = function(mcmc.list,t) {
    left.sample      = mcmc.list$left.sample
    left.pos         = mcmc.list$left.pos
    left.r           = mcmc.list$left.r
    left.lambda      = left.r*t
    right.sample     = mcmc.list$right.sample
    right.pos        = mcmc.list$right.pos
    right.r          = mcmc.list$right.r
    r.len            = length(right.r)
    right.r[r.len]   = right.r[(r.len-1)]
    right.lambda     = right.r*t
    temp.left.pos    = left.pos
    zl.on            = -left.lambda*(temp.left.pos)
    zl.on[1]         = 0
    zl.off           = log((1-exp(-left.lambda[-1]*(diff(temp.left.pos)))))
    # If there are any zeros in the recombination map, this changes them to the lowest non-zero value.
    if (length(which(zl.off==-Inf))!=0) {
        if (length(which(zl.off==-Inf))!=length(zl.off)) {
            min.value =  min(zl.off[which(zl.off!=-Inf)])-5
            zl.off[which(zl.off==-Inf)] = min.value
        }
    }
    zl               = zl.on+c(zl.off,0)
    temp.right.pos   = right.pos
    zr.on            = -right.lambda*(temp.right.pos)
    zr.on[1]         = 0
    zr.off           = log((1-exp(-right.lambda[-1]*(diff(temp.right.pos)))))
    # Same as above, but for the right side.
    if (length(which(zr.off==-Inf))!=0) {
        if (length(which(zr.off==-Inf))!=length(zr.off)) {
            min.value = min(zr.off[which(zr.off!=-Inf)])-5
            zr.off[which(zr.off==-Inf)] = min.value
        }
    }
    zr               = zr.on+c(zr.off,0)
    if (length(which(zl==-Inf))!=0) {
        zl[which(zl==-Inf)] = min(zr)
    }
    if (length(which(zr==-Inf))!=0) {
        zr[which(zr==-Inf)] = min(zl)
    }
    if (length(which(c(zr,zl)==-Inf))!=0) {
        print("The recombination map has no non-zero values - transition probabilities can't be calculated.")
        break
    }
    transition.probs = list("zl"=zl,"zr"=zr)
    return(transition.probs)
}

# Computes the emission probabilities.
emission.probs.func = function(mcmc.list,transition.probs,t.m) {
	mu = mcmc.list$mu
    r  = mcmc.list$r
    Ne = 10000
    const.mult = mcmc.list$const.mult
    bed.list = mcmc.list$bed.list
    sel.site = mcmc.list$sel.site
    # Does each side seperately.
	for (k in 1:2) {
		if (k==1) {
			sample      = mcmc.list$left.sample
			cont.sample = mcmc.list$left.cont
			anc         = mcmc.list$full.anc[sel.site:1]
			pos         = mcmc.list$left.pos
			prior       = transition.probs[[1]]
			gaps        = bed.list$left.gaps
		}
		if (k==2) {
			sample      = mcmc.list$right.sample
			cont.sample = mcmc.list$right.cont
			anc         = mcmc.list$full.anc[sel.site:length(mcmc.list$full.anc)]
			pos         = mcmc.list$right.pos
			prior       = transition.probs[[2]]
			gaps        = bed.list$right.gaps
		}
		hapcounts    = rep(1,nrow(cont.sample))
		prob.matrix  = matrix(nrow=nrow(sample),ncol=(ncol(sample)))
		mismatch.vec = 0
		mismatches   = 0
		# If a bed file was specified, this figures out how many sites to call "invariant".
		if (length(gaps)!=0) {
		    gaps = gaps[-which(gaps[,1]>pos[length(pos)]),]
	        inv.sites = diff(pos)-1
	        gap.index = c(1:nrow(gaps))
	        new.inv.sites = inv.sites
		    for (i in 1:nrow(gaps)) {
	            new.inv.sites[min(which(pos>gaps[i,2]))-1] = inv.sites[min(which(pos>gaps[i,2]))-1] - length(c(gaps[i,1]:gaps[i,2]))
		    }
	        invariant.matches = c(0,cumsum(new.inv.sites))
	    }
	    if (length(gaps)==0) {
	        inv.sites = diff(pos)-1
	        invariant.matches = c(0,cumsum(inv.sites))
	    }
	    # For each chromosome in the sample.
		for (i in 1:nrow(sample)) {
			ind                     = sample[i,]
			# Finds differences from the ancestral haplotype.
			differences             = ind-anc
			# Finds matches.
			variant.matches         = c(cumsum(differences==0))
			cumulative.matches      = c(variant.matches+invariant.matches)
			cumulative.mismatches   = cumsum(differences!=0)
			# Uses C code to perform the forwards algorithm on the reference panel.
            a.vec                   = rev(c(ind))
			H.M                     = cont.sample
			for (w in 1:nrow(cont.sample)) {
				H.M[w,] = rev(H.M[w,])
			}
			rho                     = 4*r*Ne
			Lambda.vec              <- rep(1,length(a.vec)-1)
			seqpos.vec              <- pos[length(pos)]-rev(c(pos,pos[length(pos)]))
			seqpos.vec              = seqpos.vec[-1]
			seqpos.vec[1]           = 1
			dummy.prob              <- -9.0
			copiedh.vec             = rep(0,length(a.vec))
			dummy.prob.vec          = rep(0,length(a.vec))
			result                  <- cprobback.func(a.vec, H.M, hapcounts, rho, Lambda.vec, seqpos.vec, dummy.prob, dummy.prob.vec, copiedh.vec, 1222)
			pihat.vec               = rev(result$pihat.vec)
			alpha.matches           = -t.m*mu*cumulative.matches
			alpha.mismatches        = log(1-exp(-t.m*mu))*cumulative.mismatches
			zeros                   = cumulative.mismatches==0
			alpha.mismatches[zeros] = 0
			alpha                   = alpha.matches+alpha.mismatches
			alpha[1]                = 0
			beta                    = pihat.vec
			beta[length(beta)]      = 0
			posterior               = alpha+beta+prior
			prob.matrix[i,]         = posterior
		}
		if (k==1) {
		    l.prob.table  = prob.matrix;
		    l.prob.matrix = prob.matrix;
		    l.prob.vec    = c(1:nrow(prob.matrix))
		}
		if (k==2) {
		    r.prob.table  = prob.matrix;
		    r.prob.matrix = prob.matrix;
		    r.prob.vec    = c(1:nrow(prob.matrix))
	    }
	}
	# This takes the product of likelihoods across individuals.
	for (i in 1:nrow(l.prob.table)) {
	    a = max(l.prob.table[i,])
		l.prob.matrix[i,] = exp(l.prob.table[i,] - a)
		l.prob.vec[i] = a + log(sum(l.prob.matrix[i,]))
        a = max(r.prob.table[i,])
		r.prob.matrix[i,] = exp(r.prob.table[i,] - a)
	    r.prob.vec[i] = a + log(sum(r.prob.matrix[i,]))
	}
    ind.prod.l     = sum(l.prob.vec)
    ind.prod.r     = sum(r.prob.vec)
	prob.total     = ind.prod.l+ind.prod.r
	emission.probs = list("l.prob.table" = l.prob.table,
	                      "r.prob.table" = r.prob.table,
	                      "prob.total"   = prob.total)
	return(emission.probs)
}

# Computes the first iteration of the MCMC
mcmc.init.func = function(mcmc.list) {
    print("Initializing the MCMC.")
	hap.init         = mcmc.list$full.anc
	chain.length     = mcmc.list$chain.length
        t.b              = runif(1,100,mcmc.list$upper.t.limit)
	transition.probs = transition.probs.func(mcmc.list,t.b)
	emission.probs.b = emission.probs.func(mcmc.list,transition.probs,t.b)[[3]]
	t.chain          = matrix(nrow=(chain.length+1),ncol=4)
	t.chain[1,]      = c(t.b, emission.probs.b,0,0)
	mcmc.list$t.chain  = t.chain
	return(mcmc.list)
}


startchain.func = function(mcmc.list,params) {
    print("Starting the MCMC.")
    chain.length    = params$chain.length
    nanc.post       = params$nanc.post
    proposal.sd     = params$proposal.sd
    t.chain         = mcmc.list$t.chain
    last.prop       = max(which(!is.na(t.chain[,1])))
	q.vec           = rep(0,(chain.length+1))
	evens           = seq(2,(chain.length+1),2)
	q.vec[evens]    = 1
	q.index         = 1
	anchap.post     = matrix(nrow=nanc.post,ncol=length(mcmc.list$full.anc))
	anchap.post[1,] = mcmc.list$full.anc
	anchap.index    = 1
	pb              = txtProgressBar(last.prop+1,chain.length+1,last.prop+1,style=3)
	for (y in (last.prop+1):(chain.length+1)) {
	    t.b              = t.chain[(y-1),1]
	    A.b              = mcmc.list$full.anc
	    emission.probs.b = t.chain[y-1,2]
	    A.a              = A.b
	    # propose new t
        if (q.vec[q.index]==0) {
            t.a = 0
            while(t.a<=0) {
                  t.a = t.b+rnorm(1,0,proposal.sd)
             }
        }
        # propose site to switch
        if (q.vec[q.index]==1) {
           site.to.flip       = sample(length(A.b),1)
		   A.a[site.to.flip]  = abs(A.a[site.to.flip]-1)
	       mcmc.list$full.anc = A.a
	       t.a                = t.b
	    }
	    # compute acceptance ratio
		transition.probs = transition.probs.func(mcmc.list,t.a)
		emission.probs.a = emission.probs.func(mcmc.list,transition.probs,t.a)[[3]]
		acceptance.ratio = emission.probs.a - emission.probs.b
		# accept the proposal
		if (acceptance.ratio>=0) {
		   t.b                = t.a
	       A.b                = A.a
           emission.probs.b   = emission.probs.a
           t.chain[y,]        = c(t.b,emission.probs.b,acceptance.ratio,1)
	       mcmc.list$full.anc = A.b
	    }
	    # reject the proposal
	    if (acceptance.ratio<0) {
		   accept = rbinom(1,1,exp(acceptance.ratio))
		   if (accept==1) {
		      t.b                = t.a
		      A.b                = A.a
		      emission.probs.b   = emission.probs.a
              t.chain[y,]        = c(t.b,emission.probs.b,acceptance.ratio,1)
	          mcmc.list$full.anc = A.b
	       } else {
              t.chain[y,]        = c(t.b,emission.probs.b,acceptance.ratio,0)
	          mcmc.list$full.anc = A.b
	       }
		}
		if (y%in%c(((chain.length-nanc.post)+1):chain.length)) {
		    anchap.post[anchap.index,] = A.b
		    anchap.index               = anchap.index+1
		}
	    q.index = q.index+1
        setTxtProgressBar(pb, y)
	}
	close(pb)
    mcmc.output    = list("t.chain"= t.chain,"anchap.post"=anchap.post)
    return(mcmc.output)
}

# This switches data from a list into a matrix.
no.list.func = function(data,char.or.num) {
	new.matrix = matrix(nrow=nrow(data),ncol=ncol(data))
	if (char.or.num=="num") {
	    for (i in 1:ncol(data)) {
		    new.matrix[,i] = as.numeric(paste(unlist(data[,i])))
	    }
	}
	if (char.or.num=="char") {
		for (i in 1:ncol(data)) {
		    new.matrix[,i] = as.character(paste(unlist(data[,i])))
	    }
	}
	return(new.matrix)
}
