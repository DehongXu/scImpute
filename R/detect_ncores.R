detect_ncores <- function()
{
  #number of cpus
  cpu_number_info=try(system("cat /proc/cpuinfo| grep \"physical id\"| sort| uniq| wc -l",intern=TRUE))
  num_of_cpu=as.numeric(cpu_number_info[[1]][1])
  print(paste("number of cpus: ",num_of_cpu,sep = ""))
  
  #number of cores in each cpu
  core_number_info=try(system("cat /proc/cpuinfo| grep \"cpu cores\"| uniq",intern=TRUE))
  list_str_core=strsplit(core_number_info[[1]],":")
  num_of_core=as.numeric(list_str_core[[1]][2])
  print(paste("number of cores in each cpu: ",num_of_core,sep = ""))
  
  #cpu idle rate
  hardware_info=try(system("top -bn 1 -i -c",intern=TRUE))
  cpu_info=strsplit(hardware_info[[3]],",")
  str_idle_rate=strsplit(cpu_info[[1]][4]," ")
  cpu_idle_rate=as.double(str_idle_rate[[1]][2])
  print(paste("cpu idle rate: ",cpu_idle_rate,"%",sep = ""))
  
  #cpu running threads
  tasks_info=strsplit(hardware_info[[2]],",")
  str_threads=strsplit(tasks_info[[1]][2],"running")
  num_of_threads=as.numeric(str_threads)
  print(paste("number of running threads: ",num_of_threads,sep = ""))
  
  num_of_core_idle=floor(cpu_idle_rate/100*num_of_cpu*num_of_core)
  num_of_threads_idle=num_of_cpu*num_of_core - num_of_threads
  min=which.min(c(num_of_core_idle,num_of_threads_idle))
  ncores=c(num_of_core_idle,num_of_threads_idle)[[min]]
  print(paste("ncores: ",ncores,sep=""))
  
  return(ncores)
}