#!/usr/bin/env Rscript

library(data.table)
library(biomaRt)

setGeneric(
  name="merge_enrich_terms_0.05",
  def=function(
    Input,
    cutoff=0.05,
    envir=.GlobalEnv
  ){
    standardGeneric("merge_enrich_terms_0.05")
  }
)

setMethod(
  "merge_enrich_terms_0.05",
  signature(
    Input="list"
  ),
  definition=function(Input,cutoff,envir){
    
    ## cutoff
    
    # check cutoff length according Input
    if(length(cutoff)>1 & length(cutoff)!=length(Input)){
      stop("cutoff must be a single value to repet for each list element, or the same length than list")
    }
    
    # repet cutoff value according Input length
    if(length(cutoff)==1){
      cutoff<-rep(cutoff,length(Input))
    }
    
    ## check ontology
    
    # same ontology ?
    check.onto=unique(
      unlist(
        lapply(seq_along(Input),function(x){
          
          # extract quering objects names
          x=Input[[x]]
          
          # check existence
          values<-ls(envir=envir)
          
          # check if available
          available<-x%in%values
          
          # stop if not found
          if(!all(available)){
            stop(paste("objects not found:",paste(x[!available],collapse=", "),sep="\n"))
          }
          
          # get objects
          x=mget(x,envir=envir)
          
          # objects type
          obj.type=vapply(x,class,"")
          
          # extract ontology type
          vapply(seq_along(x),function(y){
            
            # extract ontology slot
            if(obj.type[y]%in%c("topGOdata","fgsea")){
              
              # for topGOdata
              slot(x[[y]],"ontology")
              
            }else{
              
              # for topGOresult
              sub("^.+\nOntology: ","",slot(x[[y]],"description"))
            }
          },"")
        })
      )
    )
    
    # check ontology
    if(length(check.onto)>1){
      stop("Only one ontology supported")
    }
    
    # extract method used (topGO or fgsea)
    check.method=unique(
      unlist(
        lapply(seq_along(Input),function(x){
          
          # extract quering objects names
          x=Input[[x]]
          
          # get objects
          x=mget(x,envir=envir)
          
          # objects type
          vapply(x,class,"")
        })
      )
    )
    
    # same method used
    if(length(check.method)==3){
      stop("topGO enrichment and fgsea results can't be merged")
    }
    
    # if two object --> topGO
    if(length(check.method)==2){check.method<-"topGO"}
    
    ## gene symbols add functions
    genes_symbols_add=function(db,genes){
      
      # if db match to bioconductor
      if(db[1]=="Bioconductor"){
        
        # load GeneID and symbols
        annot<-data.table:::data.table(
          select(
            get(db[2]),
            keys=keys(get(db[2])),
            columns =c("ENTREZID","SYMBOL")
          )
        )
        
        # if found symbols
        if(nrow(annot[!is.na(annot$ENTREZID)])>0){
          
          # load GeneID and symbols
          genes<-merge(
            genes,
            annot,
            by.x="Significant_genes",
            by.y="ENTREZID",
            all.x=TRUE
          )
          
        }else{
          
          # add empty name columns
          genes[,Name:=NA]
        }
      }
      
      # if db match to enstrezGene
      if(db[1]=="EntrezGene"){
        
        # function for generate data packets defined by the "by" argument
        pos=function(Data,by=""){
          if(by>length(Data)){
            by=length(Data)
          }
          if(length(Data)<=by){
            data.table:::data.table(
              start=1,
              end=length(Data)
            )
          }else{
            data.table:::data.table(
              start=seq(1,length(Data),by=by),
              end=unique(c(seq(by,length(Data),by=by),length(Data)))
            )
          }
        }
        
        # pattern.extract
        pattern.extract=function(query,m){
          query=lapply(seq_along(query),function(i){
            
            if(length(na.omit(m[[i]][1]))>0){
              a=attr(m[[i]],"capture.start")
              t(substring(query[i],a,a+attr(m[[i]],"capture.length")-1))
            }else{
              NA
            }
          })
          query=do.call("rbind",query)
          query[query%in%""]<-NA
          query[query%in%"\t"]<-NA
          query
        }
        
        # esummary
        esummary<-function(...){
          
          # submitted Data
          Data<-na.omit(unique(unique(...)))
          
          # batch size
          wpos=pos(Data,by=400)
          
          # core address
          core="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?version=2.0&db="
          
          # Data submission an retrieve
          query=lapply(seq_len(nrow(wpos)),function(x){
            
            # build block of id
            Dat<-paste(
              Data[wpos[x,start]:wpos[x,end]],
              collapse=","
            )
            
            # submit block
            query <- paste(core,"gene","&id=",Dat,sep = "")
            
            # read and return reults
            readLines(query)
          })
          query=paste(
            unlist(query),
            collapse=""
          )
          
          # parse results
          query<-substring(
            query,
            unlist(gregexpr("<DocumentSummary ",query)),
            unlist(gregexpr("</DocumentSummary>",query))
          )
          
          # extraction pattern
          pattern="<DocumentSummary uid=\"(?<Id>[[:digit:]]*)\">\t<Name>(?<Name>.*)</Name>"
          
          # locate elements
          m1=gregexpr(
            paste(pattern,collapse=""),
            query,
            perl=TRUE
          )
          
          # extract elements
          res=data.table:::data.table(
            pattern.extract(query,m1)
          )
          
          # add an header
          names(res)=attr(m1[[1]],"capture.name")
          
          # return
          return(res)
        }
        
        # get esummary from NCBI
        annot<-esummary(genes$Significant_genes)
        
        # if found symbols
        if(nrow(annot[!is.na(annot$Id)])>0){
          
          # load GeneID and symbols
          genes<-merge(
            genes,
            annot,
            by.x="Significant_genes",
            by.y="Id",
            all.x=TRUE
          )
          
        }else{
          
          # add empty name columns
          genes[,Name:=NA]
        }
      }
      
      # if db  match to Ensembl
      if(db[1]=="Ensembl"){
        
        if(db[3]!="http://grch37.ensembl.org"){
          
          # connect to Ensembl
          mart<-useEnsembl(
            biomart="genes",
            host=db[3],
            version=db[6]
          )
          
        }else{
          # connect to Ensembl
          mart<-useEnsembl(
            biomart="genes",
            GRCh =37,
            version=db[6]
          )
        }
        
        # connect to ensembl specified dataset
        myspecies<-useDataset(
          db[2],
          mart
        )
        
        # load Ensembl genes, with GO annotations
        annot<-data.table:::data.table(
          getBM(
            attributes =c("ensembl_gene_id","external_gene_name"),
            value=TRUE,
            mart =myspecies
          )
        )
        
        # if found symbols
        if(nrow(annot[!"",on="external_gene_name"])>0){
          
          # load GeneID and symbols
          genes<-merge(
            genes,
            annot,
            by.x="Significant_genes",
            by.y="ensembl_gene_id",
            all.x=TRUE
          )
          
        }else{
          
          # add empty name columns
          genes[,Name:=NA]
        }
      }
      
      # if db  match to Custom
      if(db[1]=="Uniprot-GOA"){
        
        # temp file
        temp<-paste(
          tempfile(),
          "gz",
          sep="."
        )
        
        # load the file
        download.file(
          paste(
            'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/',
            toupper(db[2]),
            '/goa_',
            db[2],
            '.gaf.gz',
            sep=""
          ),
          destfile =temp,
          quiet=TRUE,
          method="internal"
        )
        
        # unzip
        gunzip(temp)
        
        # read file
        annot<-unique(
          fread(
            sub("\\.gz","",temp),
            skip=12,
            select=c(2,3),
            col.names=c("gene_id","gene_symbol")
          )
        )
        
        # merge gene_symbol
        genes<-merge(
          genes,
          annot,
          by.x="Significant_genes",
          by.y="gene_id",
          all.x=TRUE
        )
      }
      
      # if db  match to Custom
      if(db[1]=="Custom"){
        
        # load GeneID and symbols
        annot<-unique(
          fread(
            db[3],
            select=c("gene_id","gene_symbol")
          )
        )
        
        # merge gene_symbol
        genes<-merge(
          genes,
          annot,
          by.x="Significant_genes",
          by.y="gene_id",
          all.x=TRUE
        )
      }
      
      # load GeneID and symbols
      names(genes)[3]<-"Significant_genes_symbol"
      
      # reorder columns
      genes<-genes[,c("GO.ID","Significant_genes","Significant_genes_symbol"),with=FALSE]
      
      # collapse results
      genes<-genes[,lapply(.SD,function(x){paste(x,collapse=";")}),.SDcols=2:3,by="GO.ID"]
      
      # replace all blank cells by NA
      genes[genes==""]<-NA
      
      # return genes
      return(genes)
    }
    
    ## apply
    
    # topGO method
    if(check.method=="topGO"){
      
      ## topGO summary informations
      
      # build topGO summary
      obj_summary=lapply(seq_along(Input),function(x){
        
        # keep pos
        pos=x
        
        # extract quering objects names
        x=Input[[x]]
        
        # keep names
        x_names=x
        
        # get objects
        x=mget(x,envir=envir)
        
        # objects type
        obj.type=vapply(x,class,"")
        
        # extract topGO objects summary
        topGO<-lapply(seq_along(x),function(y){
          
          # extract ontology slot
          if(obj.type[y]=="topGOdata"){
            
            # for topGOdata
            list(
              
              # description
              description=slot(x[[y]],"description"),
              
              # availables genes
              available_genes=length(slot(x[[y]],"allGenes")),
              
              # availables genes significant
              available_genes_significant=table(slot(x[[y]],"allScores"))[2],
              
              # feasibles genes
              feasible_genes=table(slot(x[[y]],"feasible"))[2],
              
              # feasibles genes significant
              feasible_genes_significant=table(slot(x[[y]],"allScores")==1 & slot(x[[y]],"feasible")==TRUE)[2],
              
              # nodes with at least  x genes
              genes_nodeSize=slot(x[[y]],"nodeSize"),
              
              # number of nodes
              nodes_number=length(slot(slot(x[[y]],"graph"),"nodes")),
              
              # number of edges
              edges_number=length(slot(slot(slot(x[[y]],"graph"),"edgeData"),"data"))
            )
            
          }else{
            
            # for topGO result
            list(
              
              # description
              description=slot(x[[y]],"description"),
              
              # test name
              test_name=sub(": ","p<",slot(x[[y]],"testName")),
              
              # algorithm name
              algorithm_name=slot(x[[y]],"algorithm"),
              
              # scored GOs
              GO_scored=length(slot(x[[y]],"score")),
              
              # significant GOs according cutOff
              GO_significant=table(slot(x[[y]],"score")<cutoff[pos])[2],
              
              # feasibles genes
              feasible_genes=slot(x[[y]],"geneData")[1],
              
              # feasibles genes significant
              feasible_genes_significant=slot(x[[y]],"geneData")[2],
              
              # nodes with at least  x genes
              genes_nodeSize=slot(x[[y]],"geneData")[3],
              
              # nodes with at least  x genes
              Nontrivial_nodes=slot(x[[y]],"geneData")[4]
            )
          }
        })
        
        # extract topGO objects summary
        names(topGO)<-x_names
        
        # return topGO
        topGO
      })
      
      # add names to topGO summary
      names(obj_summary)<-names(Input)
      
      # find enrich GOs in a least one comparison
      GOs<-lapply(seq_along(Input),function(x){
        
        # objects type
        Data=mget(Input[[x]],envir=envir)
        
        ## checking step
        
        # objects type
        obj.type=vapply(Data,class,"")
        
        # objects type
        if(sum(obj.type%in%"topGOdata")>1){
          
          # stop if more than one godata by list
          stop("Only one topGOdata object is supported by list")
        }
        
        # objects type
        if(sum(obj.type%in%"topGOresult")>1){
          
          # stop if more than one godata by list
          stop("Only one topGOresult object is supported by list")
        }
        
        ## find and extract pvalues
        
        # find topGOresult
        pos=which(obj.type=="topGOresult")
        
        # extract scores
        pvalues<-topGO::score(Data[[pos]])
        
        # extract names of enrich terms
        as.vector(
          names(
            pvalues[pvalues<cutoff[x]]
          )
        )
      })
      
      # remove redondancy and convert to vector
      GOs<-as.vector(unique(unlist(GOs)))
      
      ## check genes background
      
      # extract genes background
      allgenes<-lapply(seq_along(Input),function(x){
        
        # objects type
        Data=mget(Input[[x]],envir=envir)
        
        # objects type
        obj.type=vapply(Data,class,"")
        
        # load GOdata
        pos=which(obj.type=="topGOdata")
        
        slot(Data[[pos]],"allGenes")
      })
      
      # check if same gene background
      if(length(Input)>1){
        same_genes_background=all(
          vapply(2:length(allgenes),function(x){
            identical(sort(allgenes[[1]]),sort(allgenes[[x]]))
          },TRUE)
        )
      }else{
        same_genes_background=TRUE
      }
      
      # stop if no enrich GO terms
      if(length(GOs)==0){
        stop("No enrich GO terms available in at least one condition")
      }
      
      # initialize input
      input=list()
      
      # combine results
      allResults<-lapply(seq_along(Input),function(x){
        
        ## extract Data
        
        # objects type
        Data=mget(Input[[x]],envir=envir)
        
        # objects type
        obj.type=vapply(Data,class,"")
        
        # load GOdata
        pos=which(obj.type=="topGOdata")
        
        # load GOdata
        GOdata=Data[[pos]]
        
        # tested algorithm
        algorithm<-Data[-pos]
        
        ## extract some statistics from initial GOdata object (before enrichment test)
        
        # get the GOdatatermStat(GOdata)
        Stats<-termStat(
          GOdata,
          whichGO = GOs
        )
        
        # convert to data.table
        Stats<-data.table:::data.table(
          GO.ID=row.names(Stats),
          genes_frequency=paste(
            round(Stats$Significant/Stats$Annotated*100,digits=3),
            "% (",Stats$Significant,"/",Stats$Annotated,")",
            sep=""
          )
        )
        
        ## extract genes identifiants
        
        # extract counts genes by term according GeneList
        genes<-scoresInTerm(
          GOdata,GOs,
          use.names = TRUE
        )
        
        # extract  genes Ids according GeneList
        genes=lapply(names(genes),function(x){
          
          # extract significant terms
          val=attr(genes[[x]][genes[[x]]==2],"names")
          
          # build a table
          data.table:::data.table(
            GO.ID=x,
            Significant_genes=if(length(val)>0){val}else{NA}
          )
        })
        
        # convert to data.table
        genes<-rbindlist(genes)
        
        ## add Genes symbols
        
        # get db
        db=strsplit(slot(GOdata,"description")," ")[[1]]
        
        # add genes symbols
        genes<-genes_symbols_add(db,genes)
        
        ## extract pvalue according the algorithm results
        
        # extract all pvalues from topGOresult
        pvalue<-topGO::score(algorithm[[1]])
        
        # select pvalues from topGOresult
        pvalue<-pvalue[GOs]
        
        # extract pvalue in data.table
        pvalues<-data.table:::data.table(
          GO.ID=names(pvalue),
          pvalue=as.numeric(format(pvalue,scientific = T)),
          `-log10_pvalue`=round(-log10(pvalue),digits=2)
        )
        
        # algorithm
        algo=slot(algorithm[[1]],"algorithm")
        
        # return input params
        assign(
          "input",
          c(input,list(algo)),
          inherits=TRUE
        )
        
        ## combine results
        
        # all results in list
        Results<-list(
          data.table:::data.table(GO.ID=GOs),
          Stats,
          pvalues,
          genes
        )
        
        # merge all
        Results<-Reduce(
          function(...){
            merge(
              ...,
              by ="GO.ID",
              sort=FALSE,
              all=TRUE
            )
          },
          Results
        )
        
        # remove NA in GO.Id column
        Results<-Results[!is.na(Results$GO.ID)]
        
        # Remove gene ID and symbol if GO term not significant
        Results[Results$pvalue>=cutoff[x],`:=`(Significant_genes=NA,Significant_genes_symbol=NA)]
        
        if(!is.null(names(Input))){
          
          # add GOdata name in the header
          names(Results)[-1]<-paste(
            names(Input)[x],
            names(Results)[-1],
            sep="."
          )
        }
        
        # return Results
        Results
      })
    }
    
    # fgsea method
    if(check.method=="fgsea"){
      
      # get fgsea parameters
      obj_summary=lapply(seq_along(Input),function(x){
        
        # extract quering objects names
        x=Input[[x]]
        
        # get objects
        x=get(x,envir=envir)
        
        # extract parameters
        c(
          method=slot(x,"method"),
          slot(x,"params")
        )
      })
      
      # add names to fgsea summary
      names(obj_summary)<-names(Input)
      
      # find enrich GOs in a least one comparison
      GOs<-lapply(seq_along(Input),function(x){
        
        # objects type
        Data=slot(get(Input[[x]],envir=envir),"data")[[1]]
        
        # return enrich terms
        Data[pval<cutoff[x],pathway]
      })
      
      # remove redondancy and convert to vector
      GOs<-as.vector(unique(unlist(GOs)))
      
      ## check genes background
      
      # extract genes background
      allgenes<-lapply(seq_along(Input),function(x){
        
        # load input fgsea data input
        Data=slot(get(Input[[x]],envir=envir),"input")[[1]]
        
        # return gene identifiants
        Data[,Id]
      })
      
      # check if same gene background
      if(length(Input)>1){
        same_genes_background=all(
          vapply(2:length(allgenes),function(x){
            identical(sort(allgenes[[1]]),sort(allgenes[[x]]))
          },TRUE)
        )
      }else{
        same_genes_background=TRUE
      }
      
      # stop if no enrich GO terms
      if(length(GOs)==0){
        stop("No enrich GO terms available in at least one condition")
      }
      
      # initialize input
      input=list()
      
      # combine results
      allResults<-lapply(seq_along(Input),function(x){
        
        # get Data
        Data=get(Input[[x]],envir=envir)
        
        # if fgseaSimple
        if(slot(Data,"method")=="fgseaSimple"){
          
          # build stats table
          Stats<-slot(Data,"data")[[1]][,.(GO.ID=pathway,padj,nMoreExtreme,ES,NES,size,genes_frequency)]
        }
        
        if(slot(Data,"method")=="fgseaMultilevel"){
          
          # build stats table
          Stats<-slot(Data,"data")[[1]][,.(GO.ID=pathway,padj,log2err,ES,NES,size,genes_frequency)]
        }
        
        # build genes table
        genes<-slot(Data,"data")[[1]][,.(GO.ID=pathway,Significant_genes=leadingEdge)]
        
        # unlist Significant_genes
        genes<-genes[,.(Significant_genes=unlist(Significant_genes)),by=GO.ID]
        
        # get db
        db=strsplit(slot(Data,"description")," ")[[1]]
        
        # add genes symbols
        genes<-genes_symbols_add(db,genes)
        
        ## extract pvalue according the algorithm results
        pvalues<-slot(Data,"data")[[1]][,.(GO.ID=pathway,pvalue=as.numeric(format(pval,scientific = T)),`-log10_pvalue`=round(-log10(pval),digits=2))]
        
        # algorithm
        algo=list(slot(Data,"params"))
        
        names(algo)<-names(Input)[x]
        
        # return input params
        assign(
          "input",
          c(input,algo),
          inherits=TRUE
        )
        
        ## combine results
        
        # all results in list
        Results<-list(
          data.table:::data.table(GO.ID=GOs),
          Stats,
          pvalues,
          genes
        )
        
        # merge all
        Results<-Reduce(
          function(...){
            merge(
              ...,
              by ="GO.ID",
              sort=FALSE,
              all.x=TRUE
            )
          },
          Results
        )
        
        # remove NA in GO.Id column
        Results<-Results[!is.na(Results$GO.ID)]
        
        # Remove gene ID and symbol if GO term not significant
        Results[as.numeric(Results$pvalue)>=cutoff[x],`:=`(Significant_genes=NA,Significant_genes_symbol=NA)]
        
        if(!is.null(names(Input))){
          
          # add GOdata name in the header
          names(Results)[-1]<-paste(
            names(Input)[x],
            names(Results)[-1],
            sep="."
          )
        }
        
        # return Results
        Results
      })
      
      
      
    }
    
    # number of Godata
    nb=length(allResults)
    
    # merge results if more than GOdata input
    if(nb>1){
      
      # merge all results
      allResults=Reduce(
        function(...){
          merge(...,
                by ="GO.ID",
                sort=FALSE,
                all=TRUE
          )
        },
        allResults
      )
      
    }else{
      
      # if only one GOdata unlist table
      allResults<-allResults[[1]]
    }
    
    # add Input names
    names(input)<-names(Input)
    
    # extract GO terms
    GO<-data.table:::data.table(
      select(
        GO.db,
        keys=allResults$GO.ID,
        columns = c("GOID","TERM","DEFINITION")
      )
    )
    
    # add GO term description and definition to sResults
    allResults<-data.table:::data.table(
      GO,
      allResults[,"GO.ID":=NULL]
    )
    
    # rename the first 3 columns
    names(allResults)[seq_len(3)]<-c("GO.ID","term","definition")
    
    # significant results in at least one condition
    new(
      "enrich_GO_terms",
      same_genes_background=same_genes_background,
      ont=check.onto,
      data=allResults
    )
  }
)
