write.landtnt <- function(x, names, k=NULL, filename="data_to_tnt.tnt") { 
  dimdata<-dim(x)
  species<-names
  
  if(length(dimdata)==2){
    if(is.null(k)) {
      stop("Please set k", call. = FALSE)
    } else {
      species_coor<-x
      n <- nrow(species_coor)
      p <- ncol(species_coor)/k
      rownames(species_coor)<-species
    }
  }
  
  if(length(dimdata)==3){
    n <- dimdata[1]
    k <- dimdata[2]
    p <- dimdata[3]
    species_coor<-two.d.array(x) 
  }
  
  complete_data <- NULL
  if(k==2){
    cor_x <- seq(1, by = 2, len = p)
    cor_y <- seq(2, by = 2, len = p)
    for (i in 1:n) {
      row_sp<-as.character(species_coor[i,])
      sal <-paste(row_sp[(cor_x)], "," ,row_sp[(cor_y)], sep="", collapse = " ")
      complete_data <- rbind(complete_data, sal)
    }
  }
  if(k==3){
    cor_x <- seq(1, by = 3, len = p)
    cor_y <- seq(2, by = 3, len = p)
    cor_z <- seq(3, by = 3, len = p)
    for (i in 1:n) {
      row_sp<-as.character(species_coor[i,])
      sal <-paste(row_sp[(cor_x)], "," ,row_sp[(cor_y)],",",row_sp[(cor_z)], sep="", collapse = " ")
      complete_data <- rbind(complete_data, sal)
    }
  }
  
  rownames(complete_data) <- species
  xano <- paste("1", n)
  fileConn<-file(filename)
  writeLines(c("", "xread", xano, paste0("& [landmark ",k,"d]")), fileConn)
  close(fileConn)
  for (i in 1:n) {
    line <- complete_data[i,]
    line <- paste(species[i], line)
    write(line,file=filename,append=TRUE)  
  }
  final_line <- ";"
  write(final_line,file=filename,append=TRUE) 
}