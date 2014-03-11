rm(list=ls())#wipe the workspace

#First set the parameters of the model.
threshold = c(0.9,0.91) # The minimum size required to reproduce for the two strategies with a "minimum"/small/big difference between them.
lambda = 1/threshold # Set the parameter of the exponential distribution
popsize = 1000 # preferably large to minimize drift, but also to indicate competition of a finite resource in a closed space
k=5 # seeds per unit plant size above threshold
totalGen=50 # Maximum number of generations for one simulation to run 
totalSim=1000 #System time for this is 10 mins on this machine, reduce to get the gist


#Set up the arrays for results output
Output=matrix(0,totalGen,7) #the seven outputs just for a generationloop are defined below
Outputfinal=matrix(0,totalSim,11) #eleven outputs defined below

###LOOP THE LOOP THE LOOP####
###Set random seed from a uniform distribution of values, and keep track of seeds#######
for(isim in 1:totalSim)
{
  #For deterministic results see below for pseudo random seeds, otherwise these 
  #next two lines can be commented out.
  seeder<-round(runif(min=2, max = 80E4, n=1),2)
  set.seed(seeder)
  
  # Set up a population with half one strategy, half the other (coded as "TRUE" and "FALSE")
  x = c(rep(FALSE,popsize/2),rep(TRUE,popsize/2))
  x = sort(x) # False is sorted ahead of true
  p = length(x[x == TRUE])/length(x)
  
  
  ###Loop loop###		
  ### To loop over multiple generations or just one simulation, start the loop here.##########	
  for (igen in 1:totalGen)   
  {
    Output[igen,1]<- igen     # Store Generation number in 1st column of output matrix
    #Calculate the proportion of each type for every generation and store in the Output
    Output[igen,2] <- p     # proportion of larger morph
    Output[igen,3] <- 1-p   #proportion of smaller morph
    
    # Set up a vector of sizes drawn from an exponential distribution
    # The bigger lambda is, the smaller the average size is
    size = c(rexp(length(x[x==FALSE]),lambda[1]),rexp(length(x[x==TRUE]),lambda[2]))
    Output[igen, 4]<- mean(size)
    #Before getting into the reproduction, set up arrays to keep track of fecundity, seed type and size of parents
    fecundity = rep(0,popsize)
    seeds = array()
    parentsize = array()
    
    # seeds and parentsize will start with an initial entry of NA which needs to be stripped away below
    
    # Loop over the population creating an array of seeds
    for (i in 1:length(size)) {
      if (x[i])	{	# for the "TRUE" morph (larger)...
        if (size[i] > threshold[2])		{	# If it exceeds the minimum size then... 
          fecundity[i] = floor(k*(size[i]-threshold[2]))	# it produces k seeds for every unit mass over the threshold
          # note that this produces a geometric distribution of fecundity (the discontinuous equivalent of exponential)
          # with the same lambda parameter as the size distribution
          seeds = c(seeds,rep(x[i],fecundity[i]))	# We record its seeds in the vector of progeny 
          parentsize = c(parentsize,rep(size[i],fecundity[i])) 	# and we keep track of the seed output of plants of different sizes
        }
      }
      else 	{	# If the plant is the "FALSE" morph (reproductive economy) 
        if (size[i]>threshold[1])		{	# as above, but using hte lower threshold
          fecundity[i] = floor(k*(size[i]-threshold[1]))
          seeds = c(seeds,rep(x[i],fecundity[i]))
          parentsize = c(parentsize,rep(size[i],fecundity[i]))
        }
      }
      
    }
    ###Close parent loop#####
    
    # Since we added all the seeds that started with an NA entry, we strip off the first entry and just use entries 2-n.
    seeds = seeds[2:length(seeds)]
    parentsize = parentsize[2:length(parentsize)]
    
    # Take popsize seeds at random to start the next generation
    x= sample(seeds, popsize, replace=FALSE) #If I remove popsize, this will be more than 1000 individuals and then the loop hangs itself by the seventh or 8th generation- I think
    x = sort(x)	#Sort them (False comes ahead of true)
    p = length(x[x == TRUE])/length(x)	#And calculate the frequency of the "TRUE" (larger) morph
    
    Output[igen,5]<-mean(fecundity)
    Output[igen,6]<-max(fecundity)
    Output[igen,7]<-max(size)
    
    ###Pulling other "useful" summary stats out, using logical conditions - make dataframe
    Outputdf<-as.data.frame(Output)
    maxsize<-Outputdf[Outputdf$V4==max(Outputdf$V4),] #pull out the entire row of max size
    genmaxsize<-maxsize[1,1] #pull out the Gen no in which max size occurs
    maxfec<-Outputdf[Outputdf$V5==max(Outputdf$V5),] #pull out the entire row of max fecundity
    genmaxfec<-maxfec[1,1] #pull out the Gen no in which max fecundity occurs
    fixationall<-Outputdf[Outputdf$V2>=1.0,] #get all the rows where p=1
    fixationfirstgen<-fixationall[1,1] #isolate only the first instance (generation) where p=1
    
  }
  
  #########################closes Igen loop##################################################################################################################
  
  
  
  ###The output with comments but no real names
  Outputfinal[isim,1] <-isim
  Outputfinal[isim,2] <-seeder #If deterministic only, otherwise no need to keep track
  Outputfinal[isim,3] <- mean(Output[igen,2])     # proportion of larger morph p
  Outputfinal[isim,4] <- mean(Output[igen,3])		#proportion of smaller morph q
  Outputfinal[isim,5]	<- mean(Output[igen,4])		#mean size of that simulation
  Outputfinal[isim,6] <- mean(Output[igen,5])		#mean fecundity that simulation
  Outputfinal[isim,7] <- Output[igen,6]			#max fecundity of that simulation
  Outputfinal[isim,8]	<- Output[igen,7]			#max size of that simulation
  Outputfinal[isim,9]	<- genmaxsize		#in which generation of that simulation is maximum size achieved? 
  Outputfinal[isim,10] <- genmaxfec			
  Outputfinal[isim,11] <- fixationfirstgen
  
  
}


####closes Isim loop - run 1000 - takes too much time on this machine####

#Output 
Outputfinal ##prints [1000,11] matrix with no names

#Renaming the Outputfinal columns 
Sim_No<-Outputfinal[,1]
Seed_No<-Outputfinal[,2] 
final_p_freq<-Outputfinal[,3]
final_q_freq<-Outputfinal[,4]
mean_size_sim<-Outputfinal[,5]
mean_fec_sim<-Outputfinal[,6]
max_fec_sim<-Outputfinal[,7]
max_size_sim<-Outputfinal[,8]
gen_max_size_out<-Outputfinal[,9]
gen_max_fec_out<-Outputfinal[,10]
fixation_p_1stgen<-Outputfinal[,11] 

#### Make a dataframe of these names
results.df<-data.frame(cbind(Sim_No, Seed_No,final_p_freq,
                             final_q_freq,mean_size_sim,mean_fec_sim,max_fec_sim,max_size_sim,
                             gen_max_size_out,gen_max_fec_out,
                             fixation_p_1stgen))

####Print the dataframe
results.df
summary(results.df)

####Write it out
write.csv(results.df, file="results_aarssen_very_little_difference.csv")

