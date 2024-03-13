#! /usr/bin/env -S awk --traditional -f
## Copyright (c) 2019 alberto Otero de la Roza <aoterodelaroza@gmail.com>
## Minor modifications by Rich Townsend <rhtownsend@me.com>
## This file is free software; distributed under GNU/GPL version 3.

function dirname(file){
    ## function dirname by Aleksey Cheusov
    ## https://github.com/cheusov/runawk/blob/master/modules/dirname.awk
    if (!sub(/\/[^\/]*\/?$/,"",file))
	return "."
    else if (file != "")
	return file
    else
	return "/"
}

FNR==1{
    if ((FILENAME in include) && include[FILENAME])
	file = include[FILENAME]
    else{
	file = FILENAME
	sub(/.(f|F|fpp|FPP|for|FOR|ftn|FTN|f90|F90|f95|F95|f03|F03|f08|F08)$/,"",file)
    }
}
tolower($1) == "module" && tolower($0) !~ /^[^!]+(subroutine|function|procedure)[[:blank:]]+[^!]/{
    name = tolower($2)
    sub(/!.*$/,"",name)
    mod[name]=file
}
tolower($1) == "submodule"{
    gsub(/[[:blank:]]+/,"",$0)
    gsub(/!.*$/,"",$0)
    n = split(tolower($0),arr,/[):(]/)
    name = arr[2]"@"arr[n]
    smod[name]=file
    ancestor[name] = arr[2]
    isancestor[ancestor[name]] = 1
    if (n >= 4){
	parent[name] = arr[2]"@"arr[3]
	isparent[parent[name]] = 1
    }
}
tolower($1) == "include"{
    incfile = tolower($0)
    sub(/^[[:blank:]]*include[[:blank:]]*.[[:blank:]]*/,"",incfile)
    sub(/[[:blank:]]*.[[:blank:]]*(!.*)?$/,"",incfile)
    idx = index(tolower($0),incfile)
    incfile = substr($0,idx,length(incfile))
    incfile = dirname(file)"/"incfile
    include[incfile] = file
    ARGV[ARGC++] = incfile
}
tolower($1) == "use"{
    name = tolower($0)
    sub(/^[[:blank:]]*use[[:blank:]]*/,"",name)
    sub(/^(.*::)?[[:blank:]]*/,"",name)
    sub(/[[:blank:]]*((,|!).*)?$/,"",name)
    usedmod[name]++
    fileuse[usedmod[name],name] = file
}
END{
    for (i in mod){
	## Rule 1: the anchor of a source file depends on the mod files of its modules
	printf("$(DEPEND_DIR)/%s.anc : %s.mod\n",mod[i],i)
	printf("%s.mod :\n",i)
	if ((i in isancestor) && isancestor[i]){
	    ## Rule 2: the anchor of a source file depends on the smod file of those of its modules that are ancestors of a submodule
	    printf("$(DEPEND_DIR)/%s.anc : %s.smod\n",mod[i],i)
	    printf("%s.smod :\n",i)
	}
    }

    for (i in smod){
	## Rule 3: the anchor of a source file depends on the smod file of those of its submodules that are parents of a submodule
	if ((i in isparent) && isparent[i]){
	    printf("$(DEPEND_DIR)/%s.anc : %s.smod\n",smod[i],i)
	    printf("%s.smod :\n",i)
	}

	## Rule 6: submodule anchor files depend on their ancestor's anchor files
	if ((ancestor[i] in mod) && mod[ancestor[i]] && (smod[i] != mod[ancestor[i]]))
	    printf("$(DEPEND_DIR)/%s.anc : $(DEPEND_DIR)/%s.anc\n",smod[i],mod[ancestor[i]]);

	## Rule 7: submodule anchor files depend on their parent's anchor files
	if ((i in parent) && parent[i] && smod[i] != smod[parent[i]])
	    printf("$(DEPEND_DIR)/%s.anc : $(DEPEND_DIR)/%s.anc\n",smod[i],smod[parent[i]]);
    }

    ## Rule 5: the anchor of a source file depends on the anchor of
    ## all the non-intrinsic modules it uses
    split("", filuniq)
    for (i in usedmod){
	if ((i in mod) && mod[i]){
	    for (j=1;j<=usedmod[i];j++){
		if (!filuniq[fileuse[j,i],mod[i]] && fileuse[j,i] != mod[i]){
		    filuniq[fileuse[j,i],mod[i]] = 1
		    printf("$(DEPEND_DIR)/%s.anc : $(DEPEND_DIR)/%s.anc\n",fileuse[j,i],mod[i])
		}
	    }
	}
    }

    ## Rule 4: the anchor of a source file depends on all its included files and their contents
    for (i in include)
	printf("$(DEPEND_DIR)/%s.anc : %s\n",include[i],i)
}
